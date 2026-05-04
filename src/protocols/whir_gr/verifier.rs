use crate::{
    algebra::galois_ring::{Domain, GrElem, GrError, Result},
    hash::Hash,
    protocols::whir_gr::{
        common::{WhirGrCommitment, WhirGrOpening, WhirGrPublicParameters, WhirGrRoundProof},
        constraint::{check_sumcheck_identity, sumcheck_next_sigma, WhirConstraint},
        folding::{evaluate_repeated_ternary_fold_from_values, virtual_fold_query_indices},
        merkle::{ByteMerkleTree, MerkleProof},
        multiquadratic::{pow3_checked, pow_m},
        prover::{
            indexed_label, opening_transcript, positions_to_sorted_usize,
            validate_public_parameters,
        },
        serialization::{serialize_sumcheck_polynomial, WhirGrOpeningPayload},
        transcript::Transcript,
    },
    transcript::{
        DuplexSpongeInterface, ProverMessage, VerificationError, VerificationResult, VerifierState,
    },
    verify,
};

#[derive(Clone, Debug)]
pub struct WhirGrVerifier<'a> {
    public_params: &'a WhirGrPublicParameters,
}

impl<'a> WhirGrVerifier<'a> {
    pub const fn new(public_params: &'a WhirGrPublicParameters) -> Self {
        Self { public_params }
    }

    pub fn verify(
        &self,
        commitment: &WhirGrCommitment,
        point: &[GrElem],
        opening: &WhirGrOpening,
    ) -> Result<bool> {
        if !self.shape_valid(commitment, point, opening)? {
            return Ok(false);
        }

        let ctx = &self.public_params.ctx;
        let mut transcript =
            opening_transcript(self.public_params, commitment, point, &opening.value);
        let mut constraint = WhirConstraint::new(self.public_params.ternary_grid.clone());
        constraint.add_shift_term(ctx.one(), point.to_vec())?;
        let mut sigma = opening.value.clone();
        let mut current_domain = self.public_params.initial_domain.clone();
        let mut current_root = commitment.oracle_root;
        let mut live_variables = self.public_params.variable_count;

        for (layer, &width) in self.public_params.layer_widths.iter().enumerate() {
            let round = &opening.proof.rounds[layer];
            let input = VerifierRoundInput {
                layer: layer as u64,
                width,
                live_variables,
                current_domain: &current_domain,
                current_root,
                round,
            };
            let mut round_state = VerifierRoundState {
                constraint: &mut constraint,
                sigma: &mut sigma,
            };
            let Some(next_state) = verify_round(
                self.public_params,
                &mut transcript,
                &input,
                &mut round_state,
            )?
            else {
                return Ok(false);
            };
            current_domain = next_state.domain;
            current_root = next_state.root;
            live_variables -= width;
        }

        if constraint.evaluate_w(ctx, &opening.proof.final_constant, &[])? != sigma {
            return Ok(false);
        }

        transcript.absorb_ring_element(ctx, b"whir.final.constant", &opening.proof.final_constant);
        let final_positions = transcript.derive_unique_positions(
            b"whir.final.query",
            current_domain.size(),
            self.public_params.final_repetitions,
        )?;
        if positions_to_sorted_usize(final_positions, current_domain.size())?
            != opening.proof.final_openings.queried_indices
        {
            return Ok(false);
        }
        if !ByteMerkleTree::verify(
            self.public_params.hash_id,
            current_root,
            &opening.proof.final_openings,
        )? {
            return Ok(false);
        }
        for payload in &opening.proof.final_openings.leaf_payloads {
            if ctx.deserialize(payload)? != opening.proof.final_constant {
                return Ok(false);
            }
        }

        Ok(true)
    }

    pub fn receive_commitment<H>(
        &self,
        verifier_state: &mut VerifierState<H>,
    ) -> VerificationResult<WhirGrCommitment>
    where
        H: DuplexSpongeInterface,
        Hash: ProverMessage<[H::U]>,
    {
        let oracle_root = verifier_state.prover_message()?;
        Ok(WhirGrCommitment { oracle_root })
    }

    pub fn verify_transcript<H>(
        &self,
        verifier_state: &mut VerifierState<H>,
        commitment: &WhirGrCommitment,
        point: &[GrElem],
    ) -> VerificationResult<GrElem>
    where
        H: DuplexSpongeInterface,
        WhirGrOpeningPayload: ProverMessage<[H::U]>,
    {
        let payload: WhirGrOpeningPayload = verifier_state.prover_message()?;
        let opening = payload
            .into_opening(&self.public_params.ctx)
            .map_err(|_| VerificationError)?;
        let verified = self
            .verify(commitment, point, &opening)
            .map_err(|_| VerificationError)?;
        verify!(verified);
        Ok(opening.value)
    }

    fn shape_valid(
        &self,
        commitment: &WhirGrCommitment,
        point: &[GrElem],
        opening: &WhirGrOpening,
    ) -> Result<bool> {
        validate_public_parameters(self.public_params)?;
        Ok(point.len() == self.public_params.variable_count as usize
            && commitment.oracle_root != Hash::default()
            && opening.proof.rounds.len() == self.public_params.layer_widths.len()
            && merkle_shape_valid(&opening.proof.final_openings)
            && opening.proof.rounds.iter().all(round_shape_valid))
    }
}

struct NextVerifierState {
    domain: Domain,
    root: Hash,
}

struct VerifierRoundInput<'a> {
    layer: u64,
    width: u64,
    live_variables: u64,
    current_domain: &'a Domain,
    current_root: Hash,
    round: &'a WhirGrRoundProof,
}

struct VerifierRoundState<'a> {
    constraint: &'a mut WhirConstraint,
    sigma: &'a mut GrElem,
}

fn verify_round(
    params: &WhirGrPublicParameters,
    transcript: &mut Transcript,
    input: &VerifierRoundInput<'_>,
    state: &mut VerifierRoundState<'_>,
) -> Result<Option<NextVerifierState>> {
    if input.round.sumcheck_polynomials.len() != input.width as usize
        || input.width > input.live_variables
        || input.round.g_root == Hash::default()
    {
        return Ok(None);
    }

    let ctx = &params.ctx;
    let mut alphas = Vec::with_capacity(input.width as usize);
    for j in 0..input.width {
        let h = &input.round.sumcheck_polynomials[j as usize];
        transcript.absorb_labeled_bytes(
            &indexed_label(b"whir.sumcheck.poly", input.layer, Some(j)),
            &serialize_sumcheck_polynomial(ctx, h),
        );
        let alpha = transcript
            .challenge_teichmuller(ctx, &indexed_label(b"whir.alpha", input.layer, Some(j)))?;
        if !check_sumcheck_identity(
            ctx,
            &params.ternary_grid,
            h,
            state.sigma,
            params.degree_bounds[input.layer as usize],
        )? {
            return Ok(None);
        }
        *state.sigma = sumcheck_next_sigma(ctx, h, &alpha);
        alphas.push(alpha);
    }

    transcript.absorb_labeled_bytes(
        &indexed_label(b"whir.g_root", input.layer, None),
        &input.round.g_root.0,
    );
    let fold_width = pow3_checked(input.width)?;
    let shift_domain_size = input.current_domain.size() / fold_width;
    let shift_positions = transcript.derive_unique_positions(
        &indexed_label(b"whir.shift", input.layer, None),
        shift_domain_size,
        params.shift_repetitions[input.layer as usize],
    )?;
    let parent_indices_by_shift =
        parent_indices_by_shift(input.current_domain.size(), input.width, &shift_positions)?;
    let expected_parent_indices = parent_indices_by_shift
        .iter()
        .flatten()
        .copied()
        .collect::<Vec<_>>();
    if positions_to_sorted_usize(expected_parent_indices, input.current_domain.size())?
        != input.round.virtual_fold_openings.queried_indices
    {
        return Ok(None);
    }
    if !ByteMerkleTree::verify(
        params.hash_id,
        input.current_root,
        &input.round.virtual_fold_openings,
    )? {
        return Ok(None);
    }

    let gamma =
        transcript.challenge_teichmuller(ctx, &indexed_label(b"whir.gamma", input.layer, None))?;
    let mut next_constraint = state.constraint.restrict_prefix(ctx, &alphas)?;
    let shift_domain = input.current_domain.pow_map(fold_width)?;
    let next_variable_count = input.live_variables - input.width;
    let mut gamma_power = gamma.clone();
    for (&shift_index, indices) in shift_positions.iter().zip(&parent_indices_by_shift) {
        let payloads = payloads_for_indices(&input.round.virtual_fold_openings, indices)?;
        let folded_value = evaluate_virtual_fold_query_from_payloads(
            input.current_domain,
            indices,
            &payloads,
            &alphas,
        )?;
        let shift_point = pow_m(
            ctx,
            &shift_domain.element(shift_index)?,
            next_variable_count,
        )?;
        next_constraint.add_shift_term(gamma_power.clone(), shift_point)?;
        *state.sigma = ctx.add(state.sigma, &ctx.mul(&gamma_power, &folded_value));
        gamma_power = ctx.mul(&gamma_power, &gamma);
    }

    *state.constraint = next_constraint;
    Ok(Some(NextVerifierState {
        domain: input.current_domain.pow_map(3)?,
        root: input.round.g_root,
    }))
}

fn parent_indices_by_shift(
    domain_size: u64,
    width: u64,
    shift_positions: &[u64],
) -> Result<Vec<Vec<u64>>> {
    shift_positions
        .iter()
        .map(|&shift_index| virtual_fold_query_indices(domain_size, width, shift_index))
        .collect()
}

fn payloads_for_indices(proof: &MerkleProof, indices: &[u64]) -> Result<Vec<Vec<u8>>> {
    let mut payloads = Vec::with_capacity(indices.len());
    for &index in indices {
        let index =
            usize::try_from(index).map_err(|_| GrError::ArithmeticOverflow("query index"))?;
        let position = proof
            .queried_indices
            .binary_search(&index)
            .map_err(|_| GrError::InvalidDomain("missing WHIR_GR Merkle payload for query"))?;
        payloads.push(proof.leaf_payloads[position].clone());
    }
    Ok(payloads)
}

fn evaluate_virtual_fold_query_from_payloads(
    domain: &Domain,
    parent_indices: &[u64],
    payloads: &[Vec<u8>],
    alphas: &[GrElem],
) -> Result<GrElem> {
    if payloads.len() != parent_indices.len() {
        return Err(GrError::InvalidDomain(
            "WHIR_GR virtual fold query requires one payload per parent index",
        ));
    }

    let ctx = domain.context();
    let mut points = Vec::with_capacity(parent_indices.len());
    let mut values = Vec::with_capacity(parent_indices.len());
    for (&index, payload) in parent_indices.iter().zip(payloads) {
        points.push(domain.element(index)?);
        values.push(ctx.deserialize(payload)?);
    }
    evaluate_repeated_ternary_fold_from_values(ctx, &points, &values, alphas)
}

fn round_shape_valid(round: &WhirGrRoundProof) -> bool {
    !round.sumcheck_polynomials.is_empty()
        && round.g_root != Hash::default()
        && merkle_shape_valid(&round.virtual_fold_openings)
        && round
            .sumcheck_polynomials
            .iter()
            .all(|polynomial| !polynomial.coefficients.is_empty())
}

const fn merkle_shape_valid(proof: &MerkleProof) -> bool {
    proof.leaf_count != 0
        && !proof.queried_indices.is_empty()
        && proof.queried_indices.len() == proof.leaf_payloads.len()
}
