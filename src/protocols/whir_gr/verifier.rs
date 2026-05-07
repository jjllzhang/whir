use std::time::Instant;

use crate::{
    algebra::galois_ring::{Domain, GrContext, GrElem, GrError, Result},
    hash::Hash,
    protocols::whir_gr::{
        common::{
            WhirGrCommitment, WhirGrOpening, WhirGrPublicParameters, WhirGrRoundHints,
            WhirGrRoundProof,
        },
        constraint::{check_sumcheck_identity, sumcheck_next_sigma, WhirConstraint},
        folding::{
            evaluate_ordered_repeated_ternary_fold_batch_from_values, virtual_fold_query_indices,
            virtual_fold_query_points,
        },
        merkle::{ByteMerkleTree, CompactMerkleProof, MerkleOpeningHint},
        multiquadratic::{pow3_checked, pow_m},
        prover::{
            indexed_label, opening_transcript, positions_to_sorted_usize,
            validate_public_parameters,
        },
        serialization::{
            serialize_sumcheck_polynomial, WhirGrOpeningHintPayload, WhirGrOpeningProofPayload,
        },
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

#[doc(hidden)]
#[derive(Clone, Copy, Debug, Default)]
pub struct WhirGrVerifyTimings {
    pub sumcheck_ms: f64,
    pub merkle_ms: f64,
    pub fold_ms: f64,
    pub constraint_ms: f64,
    pub final_ms: f64,
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
        let (accepted, _timings) = self.verify_impl(commitment, point, opening, false)?;
        Ok(accepted)
    }

    #[doc(hidden)]
    pub fn verify_profiled(
        &self,
        commitment: &WhirGrCommitment,
        point: &[GrElem],
        opening: &WhirGrOpening,
    ) -> Result<(bool, WhirGrVerifyTimings)> {
        self.verify_impl(commitment, point, opening, true)
    }

    fn verify_impl(
        &self,
        commitment: &WhirGrCommitment,
        point: &[GrElem],
        opening: &WhirGrOpening,
        capture_timings: bool,
    ) -> Result<(bool, WhirGrVerifyTimings)> {
        let mut timings = WhirGrVerifyTimings::default();
        if !self.shape_valid(commitment, point, opening)? {
            return Ok((false, timings));
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
                round_hints: &opening.hints.rounds[layer],
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
                &mut timings,
                capture_timings,
            )?
            else {
                return Ok((false, timings));
            };
            current_domain = next_state.domain;
            current_root = next_state.root;
            live_variables -= width;
        }

        let final_start = capture_timings.then(Instant::now);
        if constraint.evaluate_w(ctx, &opening.proof.final_constant, &[])? != sigma {
            record_elapsed(&mut timings.final_ms, final_start);
            return Ok((false, timings));
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
            record_elapsed(&mut timings.final_ms, final_start);
            return Ok((false, timings));
        }
        record_elapsed(&mut timings.final_ms, final_start);

        let merkle_start = capture_timings.then(Instant::now);
        let final_hint = MerkleOpeningHint {
            leaf_payloads: opening.hints.final_leaf_payloads.clone(),
        };
        if !ByteMerkleTree::verify_compact(
            self.public_params.hash_id,
            current_root,
            &opening.proof.final_openings,
            &final_hint,
        )? {
            record_elapsed(&mut timings.merkle_ms, merkle_start);
            return Ok((false, timings));
        }
        record_elapsed(&mut timings.merkle_ms, merkle_start);

        let final_start = capture_timings.then(Instant::now);
        for payload in &opening.hints.final_leaf_payloads {
            if ctx.deserialize(payload)? != opening.proof.final_constant {
                record_elapsed(&mut timings.final_ms, final_start);
                return Ok((false, timings));
            }
        }
        record_elapsed(&mut timings.final_ms, final_start);

        Ok((true, timings))
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
        WhirGrOpeningProofPayload: ProverMessage<[H::U]>,
    {
        let proof_payload: WhirGrOpeningProofPayload = verifier_state.prover_message()?;
        let hint_payload: WhirGrOpeningHintPayload = verifier_state.prover_hint()?;
        let (value, proof) = proof_payload
            .into_parts(&self.public_params.ctx)
            .map_err(|_| VerificationError)?;
        let hints = hint_payload.into_hints().map_err(|_| VerificationError)?;
        let opening = WhirGrOpening {
            value,
            proof,
            hints,
        };
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
            && opening.hints.rounds.len() == opening.proof.rounds.len()
            && merkle_shape_valid(&opening.proof.final_openings)
            && merkle_hint_shape_valid(
                &opening.proof.final_openings,
                &opening.hints.final_leaf_payloads,
            )
            && opening.proof.rounds.iter().all(round_shape_valid)
            && opening
                .proof
                .rounds
                .iter()
                .zip(&opening.hints.rounds)
                .all(|(round, hints)| {
                    merkle_hint_shape_valid(
                        &round.virtual_fold_openings,
                        &hints.virtual_fold_leaf_payloads,
                    )
                }))
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
    round_hints: &'a WhirGrRoundHints,
}

struct VerifierRoundState<'a> {
    constraint: &'a mut WhirConstraint,
    sigma: &'a mut GrElem,
}

struct FoldQueryCheck {
    shift_points: Vec<Vec<GrElem>>,
    folded_values: Vec<GrElem>,
}

fn verify_round(
    params: &WhirGrPublicParameters,
    transcript: &mut Transcript,
    input: &VerifierRoundInput<'_>,
    state: &mut VerifierRoundState<'_>,
    timings: &mut WhirGrVerifyTimings,
    capture_timings: bool,
) -> Result<Option<NextVerifierState>> {
    if input.round.sumcheck_polynomials.len() != input.width as usize
        || input.width > input.live_variables
        || input.round.g_root == Hash::default()
    {
        return Ok(None);
    }

    let ctx = &params.ctx;
    let mut alphas = Vec::with_capacity(input.width as usize);
    let sumcheck_start = capture_timings.then(Instant::now);
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
    record_elapsed(&mut timings.sumcheck_ms, sumcheck_start);

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
    let merkle_start = capture_timings.then(Instant::now);
    let virtual_fold_hint = MerkleOpeningHint {
        leaf_payloads: input.round_hints.virtual_fold_leaf_payloads.clone(),
    };
    if !ByteMerkleTree::verify_compact(
        params.hash_id,
        input.current_root,
        &input.round.virtual_fold_openings,
        &virtual_fold_hint,
    )? {
        record_elapsed(&mut timings.merkle_ms, merkle_start);
        return Ok(None);
    }
    record_elapsed(&mut timings.merkle_ms, merkle_start);

    let gamma =
        transcript.challenge_teichmuller(ctx, &indexed_label(b"whir.gamma", input.layer, None))?;
    let constraint_start = capture_timings.then(Instant::now);
    let mut next_constraint = state.constraint.restrict_prefix(ctx, &alphas)?;
    record_elapsed(&mut timings.constraint_ms, constraint_start);

    let fold_start = capture_timings.then(Instant::now);
    let fold_queries = check_virtual_fold_queries(
        input,
        &shift_positions,
        &parent_indices_by_shift,
        &alphas,
        fold_width,
    )?;
    record_elapsed(&mut timings.fold_ms, fold_start);

    let mut gamma_power = gamma.clone();
    for (shift_point, folded_value) in fold_queries
        .shift_points
        .into_iter()
        .zip(&fold_queries.folded_values)
    {
        let constraint_start = capture_timings.then(Instant::now);
        next_constraint.add_shift_term(gamma_power.clone(), shift_point)?;
        *state.sigma = ctx.add(state.sigma, &ctx.mul(&gamma_power, folded_value));
        gamma_power = ctx.mul(&gamma_power, &gamma);
        record_elapsed(&mut timings.constraint_ms, constraint_start);
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

fn check_virtual_fold_queries(
    input: &VerifierRoundInput<'_>,
    shift_positions: &[u64],
    parent_indices_by_shift: &[Vec<u64>],
    alphas: &[GrElem],
    fold_width: u64,
) -> Result<FoldQueryCheck> {
    let ctx = input.current_domain.context();
    let shift_domain = input.current_domain.pow_map(fold_width)?;
    let next_variable_count = input.live_variables - input.width;
    let mut query_points = Vec::with_capacity(shift_positions.len());
    let mut query_values = Vec::with_capacity(shift_positions.len());
    let mut shift_points = Vec::with_capacity(shift_positions.len());
    for (&shift_index, indices) in shift_positions.iter().zip(parent_indices_by_shift) {
        query_points.push(virtual_fold_query_points(
            input.current_domain,
            input.width,
            shift_index,
        )?);
        query_values.push(values_for_indices(
            ctx,
            &input.round.virtual_fold_openings,
            &input.round_hints.virtual_fold_leaf_payloads,
            indices,
        )?);
        shift_points.push(pow_m(
            ctx,
            &shift_domain.element(shift_index)?,
            next_variable_count,
        )?);
    }

    let folded_values = evaluate_ordered_repeated_ternary_fold_batch_from_values(
        ctx,
        query_points,
        query_values,
        alphas,
    )?;
    Ok(FoldQueryCheck {
        shift_points,
        folded_values,
    })
}

fn values_for_indices(
    ctx: &GrContext,
    proof: &CompactMerkleProof,
    leaf_payloads: &[Vec<u8>],
    indices: &[u64],
) -> Result<Vec<GrElem>> {
    let mut values = Vec::with_capacity(indices.len());
    for &index in indices {
        let index =
            usize::try_from(index).map_err(|_| GrError::ArithmeticOverflow("query index"))?;
        let position = proof
            .queried_indices
            .binary_search(&index)
            .map_err(|_| GrError::InvalidDomain("missing WHIR_GR Merkle payload for query"))?;
        values.push(ctx.deserialize(&leaf_payloads[position])?);
    }
    Ok(values)
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

fn merkle_shape_valid(proof: &CompactMerkleProof) -> bool {
    proof.leaf_count != 0
        && !proof.queried_indices.is_empty()
        && proof
            .queried_indices
            .windows(2)
            .all(|pair| pair[0] < pair[1])
}

const fn merkle_hint_shape_valid(proof: &CompactMerkleProof, leaf_payloads: &[Vec<u8>]) -> bool {
    proof.queried_indices.len() == leaf_payloads.len()
}

fn record_elapsed(slot: &mut f64, start: Option<Instant>) {
    if let Some(start) = start {
        *slot += start.elapsed().as_secs_f64() * 1000.0;
    }
}
