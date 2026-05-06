use divan::{black_box, AllocProfiler, Bencher};
use whir::protocols::whir_gr::{
    bench_support::{
        commit_input, open_input, verify_input, WhirGrBenchCase, WHIR_GR_CASES, WHIR_GR_SMALL_CASES,
    },
    prover::WhirGrProver,
    serialization::serialize_opening,
    verifier::WhirGrVerifier,
};

#[global_allocator]
static ALLOC: AllocProfiler = AllocProfiler::system();

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
            black_box(serialize_opening(&input.params.ctx, &opening));
        });
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
