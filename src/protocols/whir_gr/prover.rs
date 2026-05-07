use std::time::Instant;

use ark_std::rand::{CryptoRng, RngCore};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::{
    algebra::galois_ring::{Domain, GrContext, GrElem, GrError, Result},
    hash::Hash,
    protocols::whir_gr::{
        common::{
            WhirGrCommitment, WhirGrOpening, WhirGrProof, WhirGrProofHints, WhirGrPublicParameters,
            WhirGrRoundHints, WhirGrRoundProof,
        },
        constraint::{
            honest_sumcheck_polynomial_for_restricted, sumcheck_next_sigma, WhirConstraint,
        },
        folding::{
            evaluate_ordered_repeated_ternary_fold_batch_from_values, repeated_ternary_fold_table,
            virtual_fold_query_indices, virtual_fold_query_points,
        },
        merkle::{build_oracle_tree, ByteMerkleTree},
        multiquadratic::{
            pow2_checked, pow3_checked, pow_m, MultiQuadraticPolynomial, MultilinearPolynomial,
        },
        oracle_encoding::rs_encode_teichmuller_coset,
        serialization::{
            serialize_public_parameters, serialize_ring_vector, serialize_sumcheck_polynomial,
            WhirGrOpeningHintPayload, WhirGrOpeningProofPayload,
        },
        transcript::Transcript,
    },
    transcript::{DuplexSpongeInterface, NargSerialize, ProverMessage, ProverState},
};

#[derive(Clone, Debug)]
pub struct WhirGrCommitmentState {
    public_params: WhirGrPublicParameters,
    polynomial: MultiQuadraticPolynomial,
    initial_tree: ByteMerkleTree,
    initial_oracle: Vec<GrElem>,
    oracle_root: Hash,
}

#[doc(hidden)]
#[derive(Clone, Copy, Debug, Default)]
pub struct WhirGrCommitTimings {
    pub encode_oracle_ms: f64,
    pub merkle_ms: f64,
    pub to_multiquadratic_ms: f64,
}

#[doc(hidden)]
#[derive(Clone, Copy, Debug, Default)]
pub struct WhirGrOpenTimings {
    pub clone_ms: f64,
    pub init_ms: f64,
    pub sumcheck_ms: f64,
    pub sumcheck_constraint_plan_ms: f64,
    pub sumcheck_poly_restrict_ms: f64,
    pub sumcheck_poly_eval_ms: f64,
    pub sumcheck_accumulate_ms: f64,
    pub sumcheck_interpolate_ms: f64,
    pub restrict_ms: f64,
    pub encode_oracle_ms: f64,
    pub merkle_ms: f64,
    pub fold_ms: f64,
    pub fold_indices_ms: f64,
    pub fold_eval_ms: f64,
    pub fold_shift_points_ms: f64,
    pub merkle_open_ms: f64,
    pub constraint_ms: f64,
    pub final_ms: f64,
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

struct ShiftQueryInput<'a> {
    params: &'a WhirGrPublicParameters,
    current_domain: &'a Domain,
    current_oracle: &'a [GrElem],
    next_polynomial: &'a MultiQuadraticPolynomial,
    alphas: &'a [GrElem],
    shift_positions: &'a [u64],
    width: u64,
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
        let initial_tree = build_oracle_tree(
            self.public_params.hash_id,
            &self.public_params.ctx,
            &initial_oracle,
        )?;
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

    pub fn commit_transcript<H, R>(
        &self,
        prover_state: &mut ProverState<H, R>,
        polynomial: &MultiQuadraticPolynomial,
    ) -> Result<(WhirGrCommitment, WhirGrCommitmentState)>
    where
        H: DuplexSpongeInterface,
        R: RngCore + CryptoRng,
        Hash: ProverMessage<[H::U]>,
    {
        let (commitment, state) = self.commit(polynomial)?;
        prover_state.prover_message(&commitment.oracle_root);
        Ok((commitment, state))
    }

    pub fn commit_multilinear(
        &self,
        polynomial: &MultilinearPolynomial,
    ) -> Result<(WhirGrCommitment, WhirGrCommitmentState)> {
        let (commitment, state, _timings) = self.commit_multilinear_impl(polynomial, false)?;
        Ok((commitment, state))
    }

    #[doc(hidden)]
    pub fn commit_multilinear_profiled(
        &self,
        polynomial: &MultilinearPolynomial,
    ) -> Result<(WhirGrCommitment, WhirGrCommitmentState, WhirGrCommitTimings)> {
        self.commit_multilinear_impl(polynomial, true)
    }

    fn commit_multilinear_impl(
        &self,
        polynomial: &MultilinearPolynomial,
        capture_timings: bool,
    ) -> Result<(WhirGrCommitment, WhirGrCommitmentState, WhirGrCommitTimings)> {
        validate_public_parameters(self.public_params)?;
        if polynomial.variable_count() != self.public_params.variable_count {
            return Err(GrError::InvalidPolynomial(
                "WHIR_GR commit polynomial variable count mismatch",
            ));
        }

        let mut timings = WhirGrCommitTimings::default();
        let dense_start = capture_timings.then(Instant::now);
        let dense_polynomial = polynomial.to_multi_quadratic(&self.public_params.ctx)?;
        record_elapsed(&mut timings.to_multiquadratic_ms, dense_start);

        let encode_start = capture_timings.then(Instant::now);
        let initial_oracle = encode_multilinear_oracle_with_structured_encoding(
            self.public_params,
            &self.public_params.initial_domain,
            polynomial,
            &dense_polynomial,
        )?;
        record_elapsed(&mut timings.encode_oracle_ms, encode_start);

        let merkle_start = capture_timings.then(Instant::now);
        let initial_tree = build_oracle_tree(
            self.public_params.hash_id,
            &self.public_params.ctx,
            &initial_oracle,
        )?;
        record_elapsed(&mut timings.merkle_ms, merkle_start);
        let oracle_root = initial_tree.root();
        let commitment = WhirGrCommitment { oracle_root };

        let state = WhirGrCommitmentState {
            public_params: self.public_params.clone(),
            polynomial: dense_polynomial,
            initial_tree,
            initial_oracle,
            oracle_root,
        };
        Ok((commitment, state, timings))
    }

    pub fn commit_multilinear_transcript<H, R>(
        &self,
        prover_state: &mut ProverState<H, R>,
        polynomial: &MultilinearPolynomial,
    ) -> Result<(WhirGrCommitment, WhirGrCommitmentState)>
    where
        H: DuplexSpongeInterface,
        R: RngCore + CryptoRng,
        Hash: ProverMessage<[H::U]>,
    {
        self.commit_transcript(
            prover_state,
            &polynomial.to_multi_quadratic(&self.public_params.ctx)?,
        )
    }

    pub fn open(
        &self,
        commitment: &WhirGrCommitment,
        state: &WhirGrCommitmentState,
        point: &[GrElem],
    ) -> Result<WhirGrOpening> {
        let (opening, _timings) = self.open_impl(commitment, state, point, false)?;
        Ok(opening)
    }

    #[doc(hidden)]
    pub fn open_profiled(
        &self,
        commitment: &WhirGrCommitment,
        state: &WhirGrCommitmentState,
        point: &[GrElem],
    ) -> Result<(WhirGrOpening, WhirGrOpenTimings)> {
        self.open_impl(commitment, state, point, true)
    }

    fn open_impl(
        &self,
        commitment: &WhirGrCommitment,
        state: &WhirGrCommitmentState,
        point: &[GrElem],
        capture_timings: bool,
    ) -> Result<(WhirGrOpening, WhirGrOpenTimings)> {
        validate_open_inputs(self.public_params, commitment, state, point)?;

        let mut timings = WhirGrOpenTimings::default();
        let ctx = &self.public_params.ctx;
        let clone_start = capture_timings.then(Instant::now);
        let mut current_polynomial = state.polynomial.clone();
        let mut current_domain = self.public_params.initial_domain.clone();
        let mut current_oracle = state.initial_oracle.clone();
        let mut current_tree = state.initial_tree.clone();
        record_elapsed(&mut timings.clone_ms, clone_start);

        let init_start = capture_timings.then(Instant::now);
        let mut constraint = WhirConstraint::new(self.public_params.ternary_grid.clone());
        constraint.add_shift_term(ctx.one(), point.to_vec())?;
        let mut sigma = current_polynomial.evaluate(ctx, point)?;
        let mut transcript = opening_transcript(self.public_params, commitment, point, &sigma);
        record_elapsed(&mut timings.init_ms, init_start);

        let merkle_open_start = capture_timings.then(Instant::now);
        let (final_openings, final_opening_hint) = current_tree.open_compact(&[0])?;
        record_elapsed(&mut timings.merkle_open_ms, merkle_open_start);

        let mut opening = WhirGrOpening {
            value: sigma.clone(),
            proof: WhirGrProof {
                rounds: Vec::with_capacity(self.public_params.layer_widths.len()),
                final_constant: ctx.zero(),
                final_openings,
            },
            hints: WhirGrProofHints {
                rounds: Vec::with_capacity(self.public_params.layer_widths.len()),
                final_leaf_payloads: final_opening_hint.leaf_payloads,
            },
        };

        for (layer, &width) in self.public_params.layer_widths.iter().enumerate() {
            let mut round_state = ProverRoundState {
                polynomial: &mut current_polynomial,
                domain: &mut current_domain,
                oracle: &mut current_oracle,
                tree: &mut current_tree,
                constraint: &mut constraint,
                sigma: &mut sigma,
            };
            let (round, round_hints) = prove_round(
                self.public_params,
                &mut transcript,
                layer as u64,
                width,
                &mut round_state,
                &mut timings,
                capture_timings,
            )?;
            opening.proof.rounds.push(round);
            opening.hints.rounds.push(round_hints);
        }

        let final_start = capture_timings.then(Instant::now);
        opening.proof.final_constant = current_polynomial.evaluate(ctx, &[])?;
        transcript.absorb_ring_element(ctx, b"whir.final.constant", &opening.proof.final_constant);
        let final_positions = transcript.derive_unique_positions(
            b"whir.final.query",
            current_domain.size(),
            self.public_params.final_repetitions,
        )?;
        record_elapsed(&mut timings.final_ms, final_start);

        let merkle_open_start = capture_timings.then(Instant::now);
        let (final_openings, final_opening_hint) = current_tree.open_compact(
            &positions_to_sorted_usize(final_positions, current_domain.size())?,
        )?;
        opening.proof.final_openings = final_openings;
        opening.hints.final_leaf_payloads = final_opening_hint.leaf_payloads;
        record_elapsed(&mut timings.merkle_open_ms, merkle_open_start);
        Ok((opening, timings))
    }

    pub fn open_transcript<H, R>(
        &self,
        prover_state: &mut ProverState<H, R>,
        commitment: &WhirGrCommitment,
        state: &WhirGrCommitmentState,
        point: &[GrElem],
    ) -> Result<GrElem>
    where
        H: DuplexSpongeInterface,
        R: RngCore + CryptoRng,
        WhirGrOpeningProofPayload: ProverMessage<[H::U]>,
        WhirGrOpeningHintPayload: NargSerialize,
    {
        let opening = self.open(commitment, state, point)?;
        let value = opening.value.clone();
        prover_state.prover_message(&WhirGrOpeningProofPayload::from_opening(
            &self.public_params.ctx,
            &opening,
        ));
        prover_state.prover_hint(&WhirGrOpeningHintPayload::from_opening(&opening));
        Ok(value)
    }
}

fn prove_round(
    params: &WhirGrPublicParameters,
    transcript: &mut Transcript,
    layer: u64,
    width: u64,
    state: &mut ProverRoundState<'_>,
    timings: &mut WhirGrOpenTimings,
    capture_timings: bool,
) -> Result<(WhirGrRoundProof, WhirGrRoundHints)> {
    let ctx = &params.ctx;
    let merkle_open_start = capture_timings.then(Instant::now);
    let (virtual_fold_openings, virtual_fold_hint) = state.tree.open_compact(&[0])?;
    record_elapsed(&mut timings.merkle_open_ms, merkle_open_start);
    let mut round = WhirGrRoundProof {
        sumcheck_polynomials: Vec::with_capacity(width as usize),
        g_root: Hash::default(),
        virtual_fold_openings,
    };
    let mut round_hints = WhirGrRoundHints {
        virtual_fold_leaf_payloads: virtual_fold_hint.leaf_payloads,
    };

    let mut alphas = Vec::with_capacity(width as usize);
    let mut sumcheck_polynomial = state.polynomial.clone();
    let sumcheck_start = capture_timings.then(Instant::now);
    for j in 0..width {
        let (h, sumcheck_timings) = honest_sumcheck_polynomial_for_restricted(
            ctx,
            &sumcheck_polynomial,
            state.constraint,
            &alphas,
            capture_timings,
        )?;
        timings.sumcheck_constraint_plan_ms += sumcheck_timings.constraint_plan;
        timings.sumcheck_poly_eval_ms += sumcheck_timings.poly_eval;
        timings.sumcheck_accumulate_ms += sumcheck_timings.accumulate;
        timings.sumcheck_interpolate_ms += sumcheck_timings.interpolate;
        transcript.absorb_labeled_bytes(
            &indexed_label(b"whir.sumcheck.poly", layer, Some(j)),
            &serialize_sumcheck_polynomial(ctx, &h),
        );
        let alpha =
            transcript.challenge_teichmuller(ctx, &indexed_label(b"whir.alpha", layer, Some(j)))?;
        *state.sigma = sumcheck_next_sigma(ctx, &h, &alpha);
        let poly_restrict_start = capture_timings.then(Instant::now);
        let next_sumcheck_polynomial =
            sumcheck_polynomial.restrict_prefix(ctx, std::slice::from_ref(&alpha))?;
        record_elapsed(&mut timings.sumcheck_poly_restrict_ms, poly_restrict_start);
        alphas.push(alpha);
        sumcheck_polynomial = next_sumcheck_polynomial;
        round.sumcheck_polynomials.push(h);
    }
    record_elapsed(&mut timings.sumcheck_ms, sumcheck_start);

    let restrict_start = capture_timings.then(Instant::now);
    let next_polynomial = sumcheck_polynomial;
    let next_domain = state.domain.pow_map(3)?;
    record_elapsed(&mut timings.restrict_ms, restrict_start);

    let encode_start = capture_timings.then(Instant::now);
    let next_oracle = encode_oracle(params, &next_domain, &next_polynomial)?;
    record_elapsed(&mut timings.encode_oracle_ms, encode_start);

    let merkle_start = capture_timings.then(Instant::now);
    let next_tree = build_oracle_tree(params.hash_id, ctx, &next_oracle)?;
    record_elapsed(&mut timings.merkle_ms, merkle_start);
    round.g_root = next_tree.root();
    transcript.absorb_labeled_bytes(&indexed_label(b"whir.g_root", layer, None), &round.g_root.0);

    let shift_domain_size = state.domain.size() / pow3_checked(width)?;
    let shift_positions = transcript.derive_unique_positions(
        &indexed_label(b"whir.shift", layer, None),
        shift_domain_size,
        params.shift_repetitions[layer as usize],
    )?;
    let fold_start = capture_timings.then(Instant::now);
    let shift_input = ShiftQueryInput {
        params,
        current_domain: state.domain,
        current_oracle: state.oracle,
        next_polynomial: &next_polynomial,
        alphas: &alphas,
        shift_positions: &shift_positions,
        width,
    };
    let shift_data = shift_query_data(&shift_input, timings, capture_timings)?;
    record_elapsed(&mut timings.fold_ms, fold_start);

    let merkle_open_start = capture_timings.then(Instant::now);
    let (virtual_fold_openings, virtual_fold_hint) = state.tree.open_compact(
        &positions_to_sorted_usize(shift_data.parent_indices, state.domain.size())?,
    )?;
    round.virtual_fold_openings = virtual_fold_openings;
    round_hints.virtual_fold_leaf_payloads = virtual_fold_hint.leaf_payloads;
    record_elapsed(&mut timings.merkle_open_ms, merkle_open_start);

    let constraint_start = capture_timings.then(Instant::now);
    let gamma =
        transcript.challenge_teichmuller(ctx, &indexed_label(b"whir.gamma", layer, None))?;
    let mut next_constraint = state.constraint.restrict_prefix(ctx, &alphas)?;
    let mut gamma_power = gamma.clone();
    for (point, value) in shift_data.shift_points.iter().zip(&shift_data.shift_values) {
        next_constraint.add_shift_term(gamma_power.clone(), point.clone())?;
        *state.sigma = ctx.add(state.sigma, &ctx.mul(&gamma_power, value));
        gamma_power = ctx.mul(&gamma_power, &gamma);
    }
    record_elapsed(&mut timings.constraint_ms, constraint_start);

    *state.polynomial = next_polynomial;
    *state.domain = next_domain;
    *state.oracle = next_oracle;
    *state.tree = next_tree;
    *state.constraint = next_constraint;
    Ok((round, round_hints))
}

fn shift_query_data(
    input: &ShiftQueryInput<'_>,
    timings: &mut WhirGrOpenTimings,
    capture_timings: bool,
) -> Result<ShiftQueryData> {
    let ctx = &input.params.ctx;
    if input.current_oracle.len()
        != checked_usize(input.current_domain.size(), "current oracle size")?
    {
        return Err(GrError::InvalidDomain(
            "WHIR virtual fold query requires oracle size == domain size",
        ));
    }

    let fold_width = pow3_checked(input.width)?;
    let shift_domain = input.current_domain.pow_map(fold_width)?;
    let dense_shift_queries =
        input.shift_positions.len() as u64 >= input.current_domain.size() / fold_width / 2;
    let fold_eval_start = capture_timings.then(Instant::now);
    let folded_for_queries = if dense_shift_queries {
        repeated_ternary_fold_table(input.current_domain, input.current_oracle, input.alphas)?
    } else {
        Vec::new()
    };
    record_elapsed(&mut timings.fold_eval_ms, fold_eval_start);

    let mut parent_indices = Vec::new();
    let mut shift_points = Vec::with_capacity(input.shift_positions.len());
    let mut shift_values = if dense_shift_queries {
        Vec::with_capacity(input.shift_positions.len())
    } else {
        Vec::new()
    };
    let mut sparse_points = Vec::new();
    let mut sparse_values = Vec::new();
    for &shift_index in input.shift_positions {
        let indices_start = capture_timings.then(Instant::now);
        let indices =
            virtual_fold_query_indices(input.current_domain.size(), input.width, shift_index)?;
        parent_indices.extend_from_slice(&indices);
        record_elapsed(&mut timings.fold_indices_ms, indices_start);

        let fold_eval_start = capture_timings.then(Instant::now);
        if dense_shift_queries {
            shift_values.push(folded_for_queries[shift_index as usize].clone());
        } else {
            sparse_points.push(virtual_fold_query_points(
                input.current_domain,
                input.width,
                shift_index,
            )?);
            sparse_values.push(
                indices
                    .iter()
                    .map(|&index| input.current_oracle[index as usize].clone())
                    .collect::<Vec<_>>(),
            );
        }
        record_elapsed(&mut timings.fold_eval_ms, fold_eval_start);

        let shift_point_start = capture_timings.then(Instant::now);
        shift_points.push(pow_m(
            ctx,
            &shift_domain.element(shift_index)?,
            input.next_polynomial.variable_count(),
        )?);
        record_elapsed(&mut timings.fold_shift_points_ms, shift_point_start);
    }

    if !dense_shift_queries {
        let fold_eval_start = capture_timings.then(Instant::now);
        shift_values = evaluate_ordered_repeated_ternary_fold_batch_from_values(
            ctx,
            sparse_points,
            sparse_values,
            input.alphas,
        )?;
        record_elapsed(&mut timings.fold_eval_ms, fold_eval_start);
    }

    Ok(ShiftQueryData {
        parent_indices,
        shift_points,
        shift_values,
    })
}

pub(crate) fn encode_oracle(
    params: &WhirGrPublicParameters,
    domain: &Domain,
    polynomial: &MultiQuadraticPolynomial,
) -> Result<Vec<GrElem>> {
    if domain.context().config() != params.ctx.config() {
        return Err(GrError::DifferentRings);
    }
    if should_try_structured_oracle_encoding(domain, polynomial) {
        if let Some(oracle) =
            rs_encode_teichmuller_coset(&params.ctx, domain, polynomial.coefficients())?
        {
            return Ok(oracle);
        }
    }
    #[cfg(feature = "parallel")]
    {
        if domain.size() >= PARALLEL_ORACLE_THRESHOLD && rayon::current_num_threads() > 1 {
            return encode_oracle_parallel(params, domain, polynomial);
        }
    }
    encode_oracle_sequential(params, domain, polynomial)
}

#[cfg(test)]
pub(crate) fn encode_oracle_horner_reference(
    params: &WhirGrPublicParameters,
    domain: &Domain,
    polynomial: &MultiQuadraticPolynomial,
) -> Result<Vec<GrElem>> {
    if domain.context().config() != params.ctx.config() {
        return Err(GrError::DifferentRings);
    }
    encode_oracle_sequential(params, domain, polynomial)
}

fn should_try_structured_oracle_encoding(
    domain: &Domain,
    polynomial: &MultiQuadraticPolynomial,
) -> bool {
    let coefficient_count = polynomial.coefficients().len() as u64;
    domain.size() >= STRUCTURED_ORACLE_THRESHOLD
        && coefficient_count > 1
        && coefficient_count.saturating_mul(domain.size()) >= STRUCTURED_ORACLE_WORK_THRESHOLD
}

const STRUCTURED_ORACLE_THRESHOLD: u64 = 10_000;
const STRUCTURED_ORACLE_WORK_THRESHOLD: u64 = 16_384;

fn encode_oracle_sequential(
    params: &WhirGrPublicParameters,
    domain: &Domain,
    polynomial: &MultiQuadraticPolynomial,
) -> Result<Vec<GrElem>> {
    let size = checked_usize(domain.size(), "domain size")?;
    encode_oracle_chunk(params, domain, polynomial, 0, size)
}

#[cfg(feature = "parallel")]
fn encode_oracle_parallel(
    params: &WhirGrPublicParameters,
    domain: &Domain,
    polynomial: &MultiQuadraticPolynomial,
) -> Result<Vec<GrElem>> {
    let size = checked_usize(domain.size(), "domain size")?;
    if size == 0 {
        return Ok(Vec::new());
    }

    let target_chunks = rayon::current_num_threads().saturating_mul(4).max(1);
    let chunk_size = size.div_ceil(target_chunks).max(1);
    let starts = (0..size).step_by(chunk_size).collect::<Vec<_>>();
    let chunks = starts
        .par_iter()
        .map(|&begin| {
            let len = chunk_size.min(size - begin);
            encode_oracle_chunk(params, domain, polynomial, begin, len)
        })
        .collect::<Vec<_>>();

    let mut oracle = Vec::with_capacity(checked_usize(domain.size(), "domain size")?);
    for chunk in chunks {
        oracle.extend(chunk?);
    }
    Ok(oracle)
}

fn encode_oracle_chunk(
    params: &WhirGrPublicParameters,
    domain: &Domain,
    polynomial: &MultiQuadraticPolynomial,
    begin: usize,
    len: usize,
) -> Result<Vec<GrElem>> {
    let mut oracle = Vec::with_capacity(len);
    let begin = u64::try_from(begin).map_err(|_| GrError::ArithmeticOverflow("chunk begin"))?;
    for point in domain.iter_elements_from(begin)?.take(len) {
        oracle.push(polynomial.evaluate_pow(&params.ctx, &point));
    }
    Ok(oracle)
}

pub(crate) fn encode_multilinear_oracle(
    params: &WhirGrPublicParameters,
    domain: &Domain,
    polynomial: &MultilinearPolynomial,
) -> Result<Vec<GrElem>> {
    if domain.context().config() != params.ctx.config() {
        return Err(GrError::DifferentRings);
    }
    #[cfg(feature = "parallel")]
    {
        if domain.size() >= PARALLEL_ORACLE_THRESHOLD && rayon::current_num_threads() > 1 {
            return encode_multilinear_oracle_parallel(params, domain, polynomial);
        }
    }
    encode_multilinear_oracle_sequential(params, domain, polynomial)
}

fn encode_multilinear_oracle_with_structured_encoding(
    params: &WhirGrPublicParameters,
    domain: &Domain,
    polynomial: &MultilinearPolynomial,
    dense_polynomial: &MultiQuadraticPolynomial,
) -> Result<Vec<GrElem>> {
    if should_try_structured_oracle_encoding(domain, dense_polynomial) {
        if let Some(oracle) =
            rs_encode_teichmuller_coset(&params.ctx, domain, dense_polynomial.coefficients())?
        {
            return Ok(oracle);
        }
    }
    encode_multilinear_oracle(params, domain, polynomial)
}

const PARALLEL_ORACLE_THRESHOLD: u64 = 1024;

pub(crate) fn encode_multilinear_oracle_sequential(
    params: &WhirGrPublicParameters,
    domain: &Domain,
    polynomial: &MultilinearPolynomial,
) -> Result<Vec<GrElem>> {
    let size = checked_usize(domain.size(), "domain size")?;
    encode_multilinear_oracle_chunk(params, domain, polynomial, 0, size)
}

#[cfg(feature = "parallel")]
pub(crate) fn encode_multilinear_oracle_parallel(
    params: &WhirGrPublicParameters,
    domain: &Domain,
    polynomial: &MultilinearPolynomial,
) -> Result<Vec<GrElem>> {
    let size = checked_usize(domain.size(), "domain size")?;
    if size == 0 {
        return Ok(Vec::new());
    }

    let target_chunks = rayon::current_num_threads().saturating_mul(4).max(1);
    let chunk_size = size.div_ceil(target_chunks).max(1);
    let starts = (0..size).step_by(chunk_size).collect::<Vec<_>>();
    let chunks = starts
        .par_iter()
        .map(|&begin| {
            let len = chunk_size.min(size - begin);
            encode_multilinear_oracle_chunk(params, domain, polynomial, begin, len)
        })
        .collect::<Vec<_>>();

    let mut oracle = Vec::with_capacity(size);
    for chunk in chunks {
        oracle.extend(chunk?);
    }
    Ok(oracle)
}

fn encode_multilinear_oracle_chunk(
    params: &WhirGrPublicParameters,
    domain: &Domain,
    polynomial: &MultilinearPolynomial,
    begin: usize,
    len: usize,
) -> Result<Vec<GrElem>> {
    let mut oracle = Vec::with_capacity(len);
    let mut point_powers = Vec::with_capacity(checked_usize(
        polynomial.variable_count(),
        "variable count",
    )?);
    let mut evaluation_scratch = Vec::with_capacity(checked_usize(
        pow2_checked(polynomial.variable_count())?,
        "multilinear coefficient count",
    )?);
    let mut pow_square = params.ctx.zero();
    let mut pow_product = params.ctx.zero();
    let mut pow_scratch = vec![0; params.ctx.mul_scratch_len()];
    let mut scaled_high = params.ctx.zero();
    let mut next_value = params.ctx.zero();
    let mut eval_mul_scratch = vec![0; params.ctx.mul_scratch_len()];
    let begin = u64::try_from(begin).map_err(|_| GrError::ArithmeticOverflow("chunk begin"))?;
    for point in domain.iter_elements_from(begin)?.take(len) {
        pow_m_into(
            &params.ctx,
            &point,
            polynomial.variable_count(),
            &mut point_powers,
            &mut pow_square,
            &mut pow_product,
            &mut pow_scratch,
        )?;
        oracle.push(evaluate_multilinear_folded(
            &params.ctx,
            polynomial,
            &point_powers,
            &mut evaluation_scratch,
            &mut scaled_high,
            &mut next_value,
            &mut eval_mul_scratch,
        )?);
    }
    Ok(oracle)
}

fn pow_m_into(
    ctx: &GrContext,
    x: &GrElem,
    variable_count: u64,
    out: &mut Vec<GrElem>,
    square: &mut GrElem,
    product: &mut GrElem,
    scratch: &mut [u64],
) -> Result<()> {
    out.clear();
    out.reserve(checked_usize(variable_count, "variable count")?);
    let mut current = x.clone();
    for _ in 0..variable_count {
        out.push(current.clone());
        ctx.square_into(square, &current, scratch);
        ctx.mul_into(product, square, &current, scratch);
        std::mem::swap(&mut current, product);
    }
    Ok(())
}

fn evaluate_multilinear_folded(
    ctx: &GrContext,
    polynomial: &MultilinearPolynomial,
    point: &[GrElem],
    scratch: &mut Vec<GrElem>,
    scaled_high: &mut GrElem,
    next_value: &mut GrElem,
    mul_scratch: &mut [u64],
) -> Result<GrElem> {
    let variable_count = checked_usize(polynomial.variable_count(), "variable count")?;
    if point.len() != variable_count {
        return Err(GrError::InvalidPolynomial(
            "multilinear point length mismatch",
        ));
    }

    let coefficient_count = checked_usize(
        pow2_checked(polynomial.variable_count())?,
        "multilinear coefficient count",
    )?;
    if scratch.len() < coefficient_count {
        scratch.resize_with(coefficient_count, || ctx.zero());
    }
    for (target, coefficient) in scratch.iter_mut().zip(polynomial.coefficients()) {
        target.clone_from(coefficient);
    }
    let zero = ctx.zero();
    for target in scratch
        .iter_mut()
        .take(coefficient_count)
        .skip(polynomial.coefficients().len())
    {
        target.clone_from(&zero);
    }

    let mut active_len = coefficient_count;
    for coordinate in point {
        let next_len = active_len / 2;
        for index in 0..next_len {
            let low_index = 2 * index;
            let high_index = low_index + 1;
            if let Some(scalar) = base_scalar(&scratch[high_index]) {
                ctx.mul_base_scalar_into(scaled_high, coordinate, scalar);
            } else {
                ctx.mul_into(scaled_high, &scratch[high_index], coordinate, mul_scratch);
            }
            ctx.add_into(next_value, &scratch[low_index], scaled_high);
            std::mem::swap(&mut scratch[index], next_value);
        }
        active_len = next_len;
    }

    Ok(scratch.first().cloned().unwrap_or_else(|| ctx.zero()))
}

fn base_scalar(value: &GrElem) -> Option<u64> {
    let coefficients = value.coefficients();
    if coefficients
        .iter()
        .skip(1)
        .any(|&coefficient| coefficient != 0)
    {
        None
    } else {
        coefficients.first().copied()
    }
}

fn checked_usize(value: u64, label: &'static str) -> Result<usize> {
    usize::try_from(value).map_err(|_| GrError::ArithmeticOverflow(label))
}

fn record_elapsed(slot: &mut f64, start: Option<Instant>) {
    if let Some(start) = start {
        *slot += start.elapsed().as_secs_f64() * 1000.0;
    }
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

    #[cfg(feature = "parallel")]
    use crate::protocols::whir_gr::prover::encode_multilinear_oracle_parallel;
    use crate::{
        algebra::galois_ring::{teichmuller_generator, Domain, GrConfig, GrContext, GrElem},
        protocols::whir_gr::{
            bench_support,
            common::WhirGrPublicParameters,
            constraint::ternary_grid,
            merkle::build_oracle_tree,
            multiquadratic::{pow2_checked, MultiQuadraticPolynomial, MultilinearPolynomial},
            oracle_encoding::rs_encode_teichmuller_coset,
            prover::{
                encode_multilinear_oracle, encode_multilinear_oracle_sequential, encode_oracle,
                encode_oracle_horner_reference, WhirGrProver,
            },
            serialization::serialize_public_parameters,
            verifier::WhirGrVerifier,
        },
        transcript::{codecs::Empty, DomainSeparator, ProverState, VerifierState},
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

    fn multilinear_polynomial(
        ctx: &GrContext,
        variable_count: u64,
        seed: u64,
    ) -> MultilinearPolynomial {
        let coefficient_count = pow2_checked(variable_count).unwrap();
        let coefficients = (0..coefficient_count)
            .map(|index| ctx.from_u64((seed.wrapping_add(13 * index).wrapping_add(7)) % 29))
            .collect();
        MultilinearPolynomial::new(variable_count, coefficients).unwrap()
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
    fn profiled_opening_should_match_plain_opening() {
        let params = public_parameters(3, 1);
        let poly = polynomial(&params.ctx, 3);
        let point = open_point(&params.ctx, 3);
        let prover = WhirGrProver::new(&params);
        let (commitment, state) = prover.commit(&poly).unwrap();

        let plain_opening = prover.open(&commitment, &state, &point).unwrap();
        let (profiled_opening, _timings) =
            prover.open_profiled(&commitment, &state, &point).unwrap();

        assert_eq!(profiled_opening, plain_opening);
    }

    #[test]
    fn multilinear_oracle_encoding_should_match_dense_embedding() {
        for variable_count in 1..=3 {
            let params = public_parameters(variable_count, 1);
            for seed in [0, 19] {
                let multilinear = multilinear_polynomial(&params.ctx, variable_count, seed);
                let embedded = multilinear.to_multi_quadratic(&params.ctx).unwrap();
                assert_eq!(
                    encode_multilinear_oracle(&params, &params.initial_domain, &multilinear)
                        .unwrap(),
                    encode_oracle(&params, &params.initial_domain, &embedded).unwrap()
                );
                let prover = WhirGrProver::new(&params);
                assert_eq!(
                    prover.commit_multilinear(&multilinear).unwrap().0,
                    prover.commit(&embedded).unwrap().0
                );
            }
        }

        let case = bench_support::find_case("m4").unwrap();
        let params = bench_support::build_params(case).unwrap();
        for seed in [0, 19] {
            let multilinear =
                bench_support::multilinear_polynomial(&params.ctx, case.variable_count, seed)
                    .unwrap();
            let embedded = multilinear.to_multi_quadratic(&params.ctx).unwrap();
            assert_eq!(
                encode_multilinear_oracle(&params, &params.initial_domain, &multilinear).unwrap(),
                encode_oracle(&params, &params.initial_domain, &embedded).unwrap()
            );
            let prover = WhirGrProver::new(&params);
            assert_eq!(
                prover.commit_multilinear(&multilinear).unwrap().0,
                prover.commit(&embedded).unwrap().0
            );
        }
    }

    #[test]
    fn structured_oracle_encoding_should_match_horner_reference() {
        let case = bench_support::find_case("m6").unwrap();
        let params = bench_support::build_params(case).unwrap();
        let poly = polynomial(&params.ctx, 4);

        assert_eq!(
            rs_encode_teichmuller_coset(&params.ctx, &params.initial_domain, poly.coefficients())
                .unwrap()
                .unwrap(),
            encode_oracle_horner_reference(&params, &params.initial_domain, &poly).unwrap()
        );
    }

    #[test]
    fn structured_oracle_encoding_should_match_horner_reference_on_coset() {
        let params = public_parameters(3, 1);
        let offset = teichmuller_generator(&params.ctx).unwrap();
        let domain = Domain::teichmuller_coset(Arc::clone(&params.ctx), offset, 81).unwrap();
        let poly = polynomial(&params.ctx, 3);

        assert_eq!(
            rs_encode_teichmuller_coset(&params.ctx, &domain, poly.coefficients())
                .unwrap()
                .unwrap(),
            encode_oracle_horner_reference(&params, &domain, &poly).unwrap()
        );
    }

    #[test]
    fn structured_oracle_encoding_should_preserve_merkle_root() {
        let case = bench_support::find_case("m6").unwrap();
        let params = bench_support::build_params(case).unwrap();
        let poly = polynomial(&params.ctx, 4);
        let encoded =
            rs_encode_teichmuller_coset(&params.ctx, &params.initial_domain, poly.coefficients())
                .unwrap()
                .unwrap();
        let reference =
            encode_oracle_horner_reference(&params, &params.initial_domain, &poly).unwrap();

        assert_eq!(
            build_oracle_tree(params.hash_id, &params.ctx, &encoded)
                .unwrap()
                .root(),
            build_oracle_tree(params.hash_id, &params.ctx, &reference)
                .unwrap()
                .root()
        );
    }

    #[test]
    fn structured_oracle_encoding_should_cover_small_bench_domain_shapes() {
        let coefficient_count = 9;
        for case in bench_support::WHIR_GR_SMALL_CASES {
            let params = bench_support::build_params(case).unwrap();
            let coefficients = (0..coefficient_count)
                .map(|index| params.ctx.from_u64((case.variable_count + 17 * index) % 37))
                .collect();
            let poly = MultiQuadraticPolynomial::new(2, coefficients).unwrap();

            assert_eq!(
                rs_encode_teichmuller_coset(
                    &params.ctx,
                    &params.initial_domain,
                    poly.coefficients()
                )
                .unwrap()
                .unwrap(),
                encode_oracle_horner_reference(&params, &params.initial_domain, &poly).unwrap(),
                "{}",
                case.short_name()
            );
        }
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn parallel_multilinear_oracle_encoding_should_match_sequential() {
        let case = bench_support::find_case("m6").unwrap();
        let params = bench_support::build_params(case).unwrap();
        let multilinear =
            bench_support::multilinear_polynomial(&params.ctx, case.variable_count, 23).unwrap();
        let sequential =
            encode_multilinear_oracle_sequential(&params, &params.initial_domain, &multilinear)
                .unwrap();
        let parallel =
            encode_multilinear_oracle_parallel(&params, &params.initial_domain, &multilinear)
                .unwrap();

        assert_eq!(parallel, sequential);
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

    #[test]
    fn transcript_wrapped_opening_should_verify() {
        let params = public_parameters(2, 1);
        let poly = polynomial(&params.ctx, 2);
        let point = open_point(&params.ctx, 2);
        let prover = WhirGrProver::new(&params);
        let verifier = WhirGrVerifier::new(&params);
        let instance = Empty;
        let ds = DomainSeparator::protocol(&serialize_public_parameters(&params))
            .session(&"whir-gr transcript wrapper")
            .instance(&instance);

        let mut prover_state = ProverState::new_std(&ds);
        let (commitment, state) = prover.commit_transcript(&mut prover_state, &poly).unwrap();
        let value = prover
            .open_transcript(&mut prover_state, &commitment, &state, &point)
            .unwrap();
        let proof = prover_state.proof();

        let mut verifier_state = VerifierState::new_std(&ds, &proof);
        let received_commitment = verifier.receive_commitment(&mut verifier_state).unwrap();
        let verified_value = verifier
            .verify_transcript(&mut verifier_state, &received_commitment, &point)
            .unwrap();
        verifier_state.check_eof().unwrap();

        assert_eq!(received_commitment, commitment);
        assert_eq!(verified_value, value);
        assert_eq!(verified_value, poly.evaluate(&params.ctx, &point).unwrap());
    }
}
