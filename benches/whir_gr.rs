use std::{
    collections::HashSet,
    sync::{Mutex, OnceLock},
};

use divan::{black_box, AllocProfiler, Bencher};
use whir::protocols::whir_gr::{
    bench_support::{
        commit_input, open_input, verify_input, WhirGrBenchCase, WHIR_GR_CASES, WHIR_GR_SMALL_CASES,
    },
    prover::WhirGrProver,
    serialization::{serialize_opening, serialize_proof},
    verifier::WhirGrVerifier,
};

#[global_allocator]
static ALLOC: AllocProfiler = AllocProfiler::system();

static PRINTED_PROOF_SIZES: OnceLock<Mutex<HashSet<&'static str>>> = OnceLock::new();

#[divan::bench(args = WHIR_GR_SMALL_CASES, sample_count = 1, sample_size = 1, ignore)]
fn whir_gr_commit_small(bencher: Bencher, case: &WhirGrBenchCase) {
    bench_commit(bencher, case);
}

#[divan::bench(args = WHIR_GR_CASES, sample_count = 1, sample_size = 1, ignore)]
fn whir_gr_commit(bencher: Bencher, case: &WhirGrBenchCase) {
    bench_commit(bencher, case);
}

fn bench_commit(bencher: Bencher, case: &WhirGrBenchCase) {
    bencher
        .with_inputs(|| commit_input(case).unwrap_or_else(|error| panic!("{error}")))
        .bench_values(|input| {
            let prover = WhirGrProver::new(&input.params);
            let (commitment, state) = prover
                .commit_multilinear(&input.polynomial)
                .unwrap_or_else(|error| panic!("{error}"));
            black_box((commitment, state));
        });
}

#[divan::bench(args = WHIR_GR_CASES, sample_count = 1, sample_size = 1, ignore)]
fn whir_gr_open(bencher: Bencher, case: &WhirGrBenchCase) {
    bencher
        .with_inputs(|| open_input(case).unwrap_or_else(|error| panic!("{error}")))
        .bench_values(|input| {
            let prover = WhirGrProver::new(&input.params);
            let opening = prover
                .open(&input.commitment, &input.state, &input.point)
                .unwrap_or_else(|error| panic!("{error}"));
            let serialized_proof = serialize_proof(&input.params.ctx, &opening.proof);
            let serialized_opening = serialize_opening(&input.params.ctx, &opening);
            print_proof_size_once(case, serialized_proof.len(), serialized_opening.len());
            black_box((serialized_proof, serialized_opening));
        });
}

fn print_proof_size_once(
    case: &WhirGrBenchCase,
    proof_size_bytes: usize,
    opening_size_bytes: usize,
) {
    let printed = PRINTED_PROOF_SIZES.get_or_init(|| Mutex::new(HashSet::new()));
    let mut printed = printed
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    if printed.insert(case.name) {
        eprintln!(
            "WHIR_GR_PROOF_SIZE case={} k_exp={} r={} n={} variable_count={} max_layer_width={} lambda_target={} rho0={}/{} proof_size_bytes={} opening_size_bytes={}",
            case.name,
            case.k_exp,
            case.r,
            case.n,
            case.variable_count,
            case.max_layer_width,
            case.lambda_target,
            case.rho0.numerator,
            case.rho0.denominator,
            proof_size_bytes,
            opening_size_bytes
        );
    }
}

#[divan::bench(args = WHIR_GR_CASES, sample_count = 1, sample_size = 1, ignore)]
fn whir_gr_verify(bencher: Bencher, case: &WhirGrBenchCase) {
    bencher
        .with_inputs(|| verify_input(case).unwrap_or_else(|error| panic!("{error}")))
        .bench_values(|input| {
            let verifier = WhirGrVerifier::new(&input.params);
            let accepted = verifier
                .verify(&input.commitment, &input.point, &input.opening)
                .unwrap_or_else(|error| panic!("{error}"));
            assert!(
                accepted,
                "WHIR_GR verifier rejected an honest benchmark proof"
            );
            black_box(accepted);
        });
}

fn main() {
    divan::main();
}
