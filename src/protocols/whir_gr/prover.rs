use crate::{
    algebra::galois_ring::{Domain, GrElem, GrError, Result},
    hash::Hash,
    protocols::whir_gr::{
        common::{
            WhirGrCommitment, WhirGrOpening, WhirGrProof, WhirGrPublicParameters, WhirGrRoundProof,
        },
        constraint::{honest_sumcheck_polynomial, sumcheck_next_sigma, WhirConstraint},
        folding::{
            evaluate_repeated_ternary_fold_from_values, repeated_ternary_fold_table,
            virtual_fold_query_indices,
        },
        merkle::{build_oracle_tree, ByteMerkleTree},
        multiquadratic::{pow3_checked, pow_m, MultiQuadraticPolynomial, MultilinearPolynomial},
        serialization::{
            serialize_public_parameters, serialize_ring_vector, serialize_sumcheck_polynomial,
        },
        transcript::Transcript,
    },
};

#[derive(Clone, Debug)]
pub struct WhirGrCommitmentState {
    public_params: WhirGrPublicParameters,
    polynomial: MultiQuadraticPolynomial,
    initial_tree: ByteMerkleTree,
    initial_oracle: Vec<GrElem>,
    oracle_root: Hash,
}

#[derive(Clone, Debug)]
pub struct WhirGrProver<'a> {
    public_params: &'a WhirGrPublicParameters,
}

struct ProverRoundState<'a> {
    polynomial: &'a mut MultiQuadraticPolynomial,
    domain: &'a mut Domain,
    oracle: &'a mut Vec<GrElem>,
    tree: &'a mut ByteMerkleTree,
    constraint: &'a mut WhirConstraint,
    sigma: &'a mut GrElem,
}

struct ShiftQueryData {
    parent_indices: Vec<u64>,
    shift_points: Vec<Vec<GrElem>>,
    shift_values: Vec<GrElem>,
}

impl<'a> WhirGrProver<'a> {
    pub const fn new(public_params: &'a WhirGrPublicParameters) -> Self {
        Self { public_params }
    }

    pub fn commit(
        &self,
        polynomial: &MultiQuadraticPolynomial,
    ) -> Result<(WhirGrCommitment, WhirGrCommitmentState)> {
        validate_public_parameters(self.public_params)?;
        if polynomial.variable_count() != self.public_params.variable_count {
            return Err(GrError::InvalidPolynomial(
                "WHIR_GR commit polynomial variable count mismatch",
            ));
        }

        let initial_oracle = encode_oracle(
            self.public_params,
            &self.public_params.initial_domain,
            polynomial,
        )?;
        let initial_tree = build_oracle_tree(&self.public_params.ctx, &initial_oracle)?;
        let oracle_root = initial_tree.root();
        let commitment = WhirGrCommitment { oracle_root };
        let state = WhirGrCommitmentState {
            public_params: self.public_params.clone(),
            polynomial: polynomial.clone(),
            initial_tree,
            initial_oracle,
            oracle_root,
        };
        Ok((commitment, state))
    }

    pub fn commit_multilinear(
        &self,
        polynomial: &MultilinearPolynomial,
    ) -> Result<(WhirGrCommitment, WhirGrCommitmentState)> {
        self.commit(&polynomial.to_multi_quadratic(&self.public_params.ctx)?)
    }

    pub fn open(
        &self,
        commitment: &WhirGrCommitment,
        state: &WhirGrCommitmentState,
        point: &[GrElem],
    ) -> Result<WhirGrOpening> {
        validate_open_inputs(self.public_params, commitment, state, point)?;

        let ctx = &self.public_params.ctx;
        let mut current_polynomial = state.polynomial.clone();
        let mut current_domain = self.public_params.initial_domain.clone();
        let mut current_oracle = state.initial_oracle.clone();
        let mut current_tree = state.initial_tree.clone();
        let mut constraint = WhirConstraint::new(self.public_params.ternary_grid.clone());
        constraint.add_shift_term(ctx.one(), point.to_vec())?;
        let mut sigma = current_polynomial.evaluate(ctx, point)?;

        let mut opening = WhirGrOpening {
            value: sigma.clone(),
            proof: WhirGrProof {
                rounds: Vec::with_capacity(self.public_params.layer_widths.len()),
                final_constant: ctx.zero(),
                final_openings: current_tree.open(&[0])?,
            },
        };

        let mut transcript = opening_transcript(self.public_params, commitment, point, &sigma);
        for (layer, &width) in self.public_params.layer_widths.iter().enumerate() {
            let mut round_state = ProverRoundState {
                polynomial: &mut current_polynomial,
                domain: &mut current_domain,
                oracle: &mut current_oracle,
                tree: &mut current_tree,
                constraint: &mut constraint,
                sigma: &mut sigma,
            };
            let round = prove_round(
                self.public_params,
                &mut transcript,
                layer as u64,
                width,
                &mut round_state,
            )?;
            opening.proof.rounds.push(round);
        }

        opening.proof.final_constant = current_polynomial.evaluate(ctx, &[])?;
        transcript.absorb_ring_element(ctx, b"whir.final.constant", &opening.proof.final_constant);
        let final_positions = transcript.derive_unique_positions(
            b"whir.final.query",
            current_domain.size(),
            self.public_params.final_repetitions,
        )?;
        opening.proof.final_openings = current_tree.open(&positions_to_sorted_usize(
            final_positions,
            current_domain.size(),
        )?)?;
        Ok(opening)
    }
}

fn prove_round(
    params: &WhirGrPublicParameters,
    transcript: &mut Transcript,
    layer: u64,
    width: u64,
    state: &mut ProverRoundState<'_>,
) -> Result<WhirGrRoundProof> {
    let ctx = &params.ctx;
    let mut round = WhirGrRoundProof {
        sumcheck_polynomials: Vec::with_capacity(width as usize),
        g_root: Hash::default(),
        virtual_fold_openings: state.tree.open(&[0])?,
    };

    let mut alphas = Vec::with_capacity(width as usize);
    for j in 0..width {
        let h = honest_sumcheck_polynomial(ctx, state.polynomial, state.constraint, &alphas)?;
        transcript.absorb_labeled_bytes(
            &indexed_label(b"whir.sumcheck.poly", layer, Some(j)),
            &serialize_sumcheck_polynomial(ctx, &h),
        );
        let alpha =
            transcript.challenge_teichmuller(ctx, &indexed_label(b"whir.alpha", layer, Some(j)))?;
        *state.sigma = sumcheck_next_sigma(ctx, &h, &alpha);
        alphas.push(alpha);
        round.sumcheck_polynomials.push(h);
    }

    let next_polynomial = state.polynomial.restrict_prefix(ctx, &alphas)?;
    let next_domain = state.domain.pow_map(3)?;
    let next_oracle = encode_oracle(params, &next_domain, &next_polynomial)?;
    let next_tree = build_oracle_tree(ctx, &next_oracle)?;
    round.g_root = next_tree.root();
    transcript.absorb_labeled_bytes(&indexed_label(b"whir.g_root", layer, None), &round.g_root.0);

    let shift_domain_size = state.domain.size() / pow3_checked(width)?;
    let shift_positions = transcript.derive_unique_positions(
        &indexed_label(b"whir.shift", layer, None),
        shift_domain_size,
        params.shift_repetitions[layer as usize],
    )?;
    let shift_data = shift_query_data(
        params,
        state.domain,
        state.oracle,
        &next_polynomial,
        &alphas,
        &shift_positions,
        width,
    )?;
    round.virtual_fold_openings = state.tree.open(&positions_to_sorted_usize(
        shift_data.parent_indices,
        state.domain.size(),
    )?)?;

    let gamma =
        transcript.challenge_teichmuller(ctx, &indexed_label(b"whir.gamma", layer, None))?;
    let mut next_constraint = state.constraint.restrict_prefix(ctx, &alphas)?;
    let mut gamma_power = gamma.clone();
    for (point, value) in shift_data.shift_points.iter().zip(&shift_data.shift_values) {
        next_constraint.add_shift_term(gamma_power.clone(), point.clone())?;
        *state.sigma = ctx.add(state.sigma, &ctx.mul(&gamma_power, value));
        gamma_power = ctx.mul(&gamma_power, &gamma);
    }

    *state.polynomial = next_polynomial;
    *state.domain = next_domain;
    *state.oracle = next_oracle;
    *state.tree = next_tree;
    *state.constraint = next_constraint;
    Ok(round)
}

fn shift_query_data(
    params: &WhirGrPublicParameters,
    current_domain: &Domain,
    current_oracle: &[GrElem],
    next_polynomial: &MultiQuadraticPolynomial,
    alphas: &[GrElem],
    shift_positions: &[u64],
    width: u64,
) -> Result<ShiftQueryData> {
    let ctx = &params.ctx;
    let fold_width = pow3_checked(width)?;
    let shift_domain = current_domain.pow_map(fold_width)?;
    let dense_shift_queries =
        shift_positions.len() as u64 >= current_domain.size() / fold_width / 2;
    let folded_for_queries = if dense_shift_queries {
        repeated_ternary_fold_table(current_domain, current_oracle, alphas)?
    } else {
        Vec::new()
    };

    let mut parent_indices = Vec::new();
    let mut shift_points = Vec::with_capacity(shift_positions.len());
    let mut shift_values = Vec::with_capacity(shift_positions.len());
    for &shift_index in shift_positions {
        let indices = virtual_fold_query_indices(current_domain.size(), width, shift_index)?;
        parent_indices.extend_from_slice(&indices);
        let value = if dense_shift_queries {
            folded_for_queries[shift_index as usize].clone()
        } else {
            evaluate_virtual_fold_query_from_oracle(
                current_domain,
                current_oracle,
                &indices,
                alphas,
            )?
        };
        shift_values.push(value);
        shift_points.push(pow_m(
            ctx,
            &shift_domain.element(shift_index)?,
            next_polynomial.variable_count(),
        )?);
    }

    Ok(ShiftQueryData {
        parent_indices,
        shift_points,
        shift_values,
    })
}

fn evaluate_virtual_fold_query_from_oracle(
    domain: &Domain,
    oracle: &[GrElem],
    indices: &[u64],
    alphas: &[GrElem],
) -> Result<GrElem> {
    if oracle.len() != domain.size() as usize {
        return Err(GrError::InvalidDomain(
            "WHIR virtual fold query requires oracle size == domain size",
        ));
    }

    let mut points = Vec::with_capacity(indices.len());
    let mut values = Vec::with_capacity(indices.len());
    for &index in indices {
        points.push(domain.element(index)?);
        values.push(oracle[index as usize].clone());
    }
    evaluate_repeated_ternary_fold_from_values(domain.context(), &points, &values, alphas)
}

pub(crate) fn encode_oracle(
    params: &WhirGrPublicParameters,
    domain: &Domain,
    polynomial: &MultiQuadraticPolynomial,
) -> Result<Vec<GrElem>> {
    if domain.context().config() != params.ctx.config() {
        return Err(GrError::DifferentRings);
    }
    (0..domain.size())
        .map(|index| Ok(polynomial.evaluate_pow(&params.ctx, &domain.element(index)?)))
        .collect()
}

pub(crate) fn opening_transcript(
    params: &WhirGrPublicParameters,
    commitment: &WhirGrCommitment,
    point: &[GrElem],
    value: &GrElem,
) -> Transcript {
    let ctx = &params.ctx;
    let mut transcript = Transcript::new(b"whir-gr.opening.v1");
    transcript.absorb_labeled_bytes(b"whir.pp", &serialize_public_parameters(params));
    transcript.absorb_labeled_bytes(b"whir.commitment", &commitment.oracle_root.0);
    transcript.absorb_labeled_bytes(b"whir.open.point", &serialize_ring_vector(ctx, point));
    transcript.absorb_ring_element(ctx, b"whir.open.value", value);
    transcript
}

pub(crate) fn indexed_label(prefix: &[u8], round: u64, index: Option<u64>) -> Vec<u8> {
    let mut label = Vec::from(prefix);
    label.push(b':');
    label.extend_from_slice(round.to_string().as_bytes());
    if let Some(index) = index {
        label.push(b':');
        label.extend_from_slice(index.to_string().as_bytes());
    }
    label
}

pub(crate) fn positions_to_sorted_usize(mut indices: Vec<u64>, size: u64) -> Result<Vec<usize>> {
    indices.sort_unstable();
    indices.dedup();
    let mut out = Vec::with_capacity(indices.len());
    for index in indices {
        if index >= size {
            return Err(GrError::IndexOutOfRange { index, size });
        }
        out.push(index as usize);
    }
    Ok(out)
}

pub(crate) fn validate_public_parameters(params: &WhirGrPublicParameters) -> Result<()> {
    if params.ctx.config().p != 2
        || params.variable_count == 0
        || params.lambda_target == 0
        || params.layer_widths.is_empty()
    {
        return Err(GrError::InvalidDomain("invalid WHIR_GR public parameters"));
    }
    if params.initial_domain.context().config() != params.ctx.config()
        || !params.initial_domain.is_teichmuller_subset()
    {
        return Err(GrError::InvalidDomain(
            "WHIR_GR initial domain must belong to the parameter ring",
        ));
    }
    if params.shift_repetitions.len() != params.layer_widths.len()
        || params.degree_bounds.len() != params.layer_widths.len()
    {
        return Err(GrError::InvalidDomain(
            "WHIR_GR layer metadata length mismatch",
        ));
    }
    validate_layer_schedule(params)?;
    validate_ternary_grid(params)
}

fn validate_layer_schedule(params: &WhirGrPublicParameters) -> Result<()> {
    let mut summed_width = 0u64;
    let mut live_variables = params.variable_count;
    let mut domain_size = params.initial_domain.size();
    for (layer, &width) in params.layer_widths.iter().enumerate() {
        if width == 0
            || width > live_variables
            || params.degree_bounds[layer] == 0
            || params.shift_repetitions[layer] == 0
        {
            return Err(GrError::InvalidDomain("invalid WHIR_GR layer schedule"));
        }
        let width_divisor = pow3_checked(width)?;
        let live_code_size = pow3_checked(live_variables)?;
        if !domain_size.is_multiple_of(width_divisor) || live_code_size >= domain_size {
            return Err(GrError::InvalidDomain(
                "invalid WHIR_GR layer domain/rate schedule",
            ));
        }
        summed_width += width;
        live_variables -= width;
        if !domain_size.is_multiple_of(3) {
            return Err(GrError::InvalidDomain(
                "WHIR_GR domain must support ternary folding",
            ));
        }
        domain_size /= 3;
    }
    if summed_width != params.variable_count || live_variables != 0 {
        return Err(GrError::InvalidDomain(
            "WHIR_GR layer widths must consume all variables",
        ));
    }
    Ok(())
}

fn validate_ternary_grid(params: &WhirGrPublicParameters) -> Result<()> {
    let ctx = &params.ctx;
    if !ctx.is_unit(&params.omega)
        || params.omega == ctx.one()
        || ctx.pow(&params.omega, 3) != ctx.one()
        || params.ternary_grid[0] != ctx.one()
        || params.ternary_grid[1] != params.omega
        || params.ternary_grid[2] != ctx.square(&params.omega)
    {
        return Err(GrError::InvalidDomain("invalid WHIR_GR ternary grid"));
    }
    Ok(())
}

fn validate_open_inputs(
    params: &WhirGrPublicParameters,
    commitment: &WhirGrCommitment,
    state: &WhirGrCommitmentState,
    point: &[GrElem],
) -> Result<()> {
    validate_public_parameters(params)?;
    if point.len() != params.variable_count as usize {
        return Err(GrError::InvalidPolynomial(
            "WHIR_GR open point length mismatch",
        ));
    }
    if state.public_params.ctx.config() != params.ctx.config()
        || state.public_params.variable_count != params.variable_count
        || state.oracle_root != commitment.oracle_root
        || state.initial_tree.root() != commitment.oracle_root
        || state.initial_tree.leaf_count() != params.initial_domain.size() as usize
        || state.initial_oracle.len() != params.initial_domain.size() as usize
    {
        return Err(GrError::InvalidDomain(
            "WHIR_GR open state does not match commitment",
        ));
    }
    if params.final_repetitions == 0 {
        return Err(GrError::InvalidDomain(
            "WHIR_GR open requires final_repetitions > 0",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::{
        algebra::galois_ring::{Domain, GrConfig, GrContext, GrElem},
        protocols::whir_gr::{
            common::WhirGrPublicParameters, constraint::ternary_grid,
            multiquadratic::MultiQuadraticPolynomial, prover::WhirGrProver,
            verifier::WhirGrVerifier,
        },
    };

    fn context_for_domain_size(domain_size: u64) -> Arc<GrContext> {
        let r = match domain_size {
            9 => 6,
            27 => 18,
            81 => 54,
            _ => panic!("unsupported test domain size"),
        };
        Arc::new(GrContext::new(GrConfig { p: 2, k_exp: 16, r }).unwrap())
    }

    fn layer_widths(variable_count: u64, max_layer_width: u64) -> Vec<u64> {
        let mut widths = Vec::new();
        let mut remaining = variable_count;
        while remaining != 0 {
            let width = max_layer_width.min(remaining);
            widths.push(width);
            remaining -= width;
        }
        widths
    }

    fn public_parameters(variable_count: u64, max_layer_width: u64) -> WhirGrPublicParameters {
        let domain_size = 3u64.pow((variable_count + 1) as u32);
        let ctx = context_for_domain_size(domain_size);
        let domain = Domain::teichmuller_subgroup(Arc::clone(&ctx), domain_size).unwrap();
        let omega = ctx.pow(domain.root(), u128::from(domain_size / 3));
        let grid = ternary_grid(&ctx, &omega).unwrap();
        let widths = layer_widths(variable_count, max_layer_width);
        let mut params =
            WhirGrPublicParameters::new(Arc::clone(&ctx), domain, variable_count, omega, grid);
        params.layer_widths = widths.clone();
        params.shift_repetitions = vec![1; widths.len()];
        params.final_repetitions = 1;
        params.degree_bounds = vec![4; widths.len()];
        params.lambda_target = 32;
        params
    }

    fn polynomial(ctx: &GrContext, variable_count: u64) -> MultiQuadraticPolynomial {
        let coefficient_count = 3u64.pow(variable_count as u32);
        let coefficients = (0..coefficient_count)
            .map(|index| ctx.from_u64((11 * index + 5) % 23))
            .collect();
        MultiQuadraticPolynomial::new(variable_count, coefficients).unwrap()
    }

    fn open_point(ctx: &GrContext, variable_count: u64) -> Vec<GrElem> {
        (0..variable_count)
            .map(|index| ctx.from_u64(7 + 3 * index))
            .collect()
    }

    fn run_roundtrip(variable_count: u64, max_layer_width: u64) {
        let params = public_parameters(variable_count, max_layer_width);
        let poly = polynomial(&params.ctx, variable_count);
        let point = open_point(&params.ctx, variable_count);
        let prover = WhirGrProver::new(&params);
        let verifier = WhirGrVerifier::new(&params);

        let (commitment, state) = prover.commit(&poly).unwrap();
        let opening = prover.open(&commitment, &state, &point).unwrap();

        assert!(verifier.verify(&commitment, &point, &opening).unwrap());
        assert_eq!(opening.value, poly.evaluate(&params.ctx, &point).unwrap());
        assert_eq!(opening.proof.rounds.len(), params.layer_widths.len());
        for (round, &width) in opening.proof.rounds.iter().zip(&params.layer_widths) {
            assert_eq!(round.sumcheck_polynomials.len(), width as usize);
        }
    }

    #[test]
    fn honest_roundtrips_should_verify() {
        run_roundtrip(1, 1);
        run_roundtrip(2, 1);
        run_roundtrip(3, 1);
        run_roundtrip(2, 2);
    }

    #[test]
    fn verifier_should_reject_tampering() {
        let params = public_parameters(2, 1);
        let poly = polynomial(&params.ctx, 2);
        let point = open_point(&params.ctx, 2);
        let prover = WhirGrProver::new(&params);
        let verifier = WhirGrVerifier::new(&params);
        let (commitment, state) = prover.commit(&poly).unwrap();
        let opening = prover.open(&commitment, &state, &point).unwrap();

        let mut bad_commitment = commitment.clone();
        bad_commitment.oracle_root.0[0] ^= 1;
        assert!(!verifier.verify(&bad_commitment, &point, &opening).unwrap());

        let mut bad_opening = opening.clone();
        bad_opening.value = params.ctx.add(&bad_opening.value, &params.ctx.one());
        assert!(!verifier.verify(&commitment, &point, &bad_opening).unwrap());

        let mut bad_final = opening;
        bad_final.proof.final_constant = params
            .ctx
            .add(&bad_final.proof.final_constant, &params.ctx.one());
        assert!(!verifier.verify(&commitment, &point, &bad_final).unwrap());
    }
}
