use std::{error::Error, sync::Arc};

use divan::{black_box, AllocProfiler, Bencher};
use whir::{
    algebra::galois_ring::{Domain, GrConfig, GrContext, GrElem},
    hash::BLAKE3,
    protocols::whir_gr::{
        common::{WhirGrCommitment, WhirGrOpening, WhirGrPublicParameters},
        constraint::ternary_grid,
        multiquadratic::{pow2_checked, MultilinearPolynomial},
        prover::{WhirGrCommitmentState, WhirGrProver},
        serialization::serialize_opening,
        soundness::{
            select_whir_unique_decoding_parameters, WhirRational, WhirUniqueDecodingInputs,
        },
        verifier::WhirGrVerifier,
    },
};

#[global_allocator]
static ALLOC: AllocProfiler = AllocProfiler::system();

#[derive(Clone, Copy, Debug)]
struct WhirGrBenchCase {
    name: &'static str,
    k_exp: u32,
    r: u64,
    n: u64,
    variable_count: u64,
    max_layer_width: u64,
    lambda_target: u64,
    rho0: WhirRational,
}

struct CommitInput {
    params: WhirGrPublicParameters,
    polynomial: MultilinearPolynomial,
}

struct OpenInput {
    params: WhirGrPublicParameters,
    commitment: WhirGrCommitment,
    state: WhirGrCommitmentState,
    point: Vec<GrElem>,
}

struct VerifyInput {
    params: WhirGrPublicParameters,
    commitment: WhirGrCommitment,
    point: Vec<GrElem>,
    opening: WhirGrOpening,
}

const WHIR_GR_CASES: &[WhirGrBenchCase] = &[
    WhirGrBenchCase {
        name: "gr216_r162_m4_multilinear",
        k_exp: 16,
        r: 162,
        n: 189,
        variable_count: 4,
        max_layer_width: 3,
        lambda_target: 128,
        rho0: WhirRational {
            numerator: 1,
            denominator: 2,
        },
    },
    WhirGrBenchCase {
        name: "gr216_r162_m5_multilinear",
        k_exp: 16,
        r: 162,
        n: 513,
        variable_count: 5,
        max_layer_width: 3,
        lambda_target: 128,
        rho0: WhirRational {
            numerator: 1,
            denominator: 2,
        },
    },
    WhirGrBenchCase {
        name: "gr216_r162_m6_multilinear",
        k_exp: 16,
        r: 162,
        n: 1539,
        variable_count: 6,
        max_layer_width: 3,
        lambda_target: 128,
        rho0: WhirRational {
            numerator: 1,
            denominator: 2,
        },
    },
    WhirGrBenchCase {
        name: "gr216_r162_m7_multilinear",
        k_exp: 16,
        r: 162,
        n: 4617,
        variable_count: 7,
        max_layer_width: 3,
        lambda_target: 128,
        rho0: WhirRational {
            numerator: 1,
            denominator: 2,
        },
    },
    WhirGrBenchCase {
        name: "gr216_r162_m8_multilinear",
        k_exp: 16,
        r: 162,
        n: 13203,
        variable_count: 8,
        max_layer_width: 3,
        lambda_target: 128,
        rho0: WhirRational {
            numerator: 1,
            denominator: 2,
        },
    },
    WhirGrBenchCase {
        name: "gr216_r162_m9_multilinear",
        k_exp: 16,
        r: 162,
        n: 39609,
        variable_count: 9,
        max_layer_width: 3,
        lambda_target: 128,
        rho0: WhirRational {
            numerator: 1,
            denominator: 2,
        },
    },
    WhirGrBenchCase {
        name: "gr216_r162_m10_multilinear",
        k_exp: 16,
        r: 162,
        n: 124_173,
        variable_count: 10,
        max_layer_width: 3,
        lambda_target: 128,
        rho0: WhirRational {
            numerator: 1,
            denominator: 2,
        },
    },
];

#[divan::bench(args = WHIR_GR_CASES, sample_count = 1, sample_size = 1, ignore)]
fn whir_gr_commit(bencher: Bencher, case: &WhirGrBenchCase) {
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

fn commit_input(case: &WhirGrBenchCase) -> Result<CommitInput, Box<dyn Error>> {
    let params = build_params(case)?;
    let polynomial = multilinear_polynomial(&params.ctx, case.variable_count, 0)?;
    Ok(CommitInput { params, polynomial })
}

fn open_input(case: &WhirGrBenchCase) -> Result<OpenInput, Box<dyn Error>> {
    let input = commit_input(case)?;
    let prover = WhirGrProver::new(&input.params);
    let (commitment, state) = prover.commit_multilinear(&input.polynomial)?;
    let point = open_point(&input.params.ctx, case.variable_count, 0);
    Ok(OpenInput {
        params: input.params,
        commitment,
        state,
        point,
    })
}

fn verify_input(case: &WhirGrBenchCase) -> Result<VerifyInput, Box<dyn Error>> {
    let input = open_input(case)?;
    let prover = WhirGrProver::new(&input.params);
    let opening = prover.open(&input.commitment, &input.state, &input.point)?;
    Ok(VerifyInput {
        params: input.params,
        commitment: input.commitment,
        point: input.point,
        opening,
    })
}

fn build_params(case: &WhirGrBenchCase) -> Result<WhirGrPublicParameters, Box<dyn Error>> {
    let selection = select_whir_unique_decoding_parameters(&WhirUniqueDecodingInputs {
        lambda_target: case.lambda_target,
        ring_exponent: u64::from(case.k_exp),
        variable_count: case.variable_count,
        max_layer_width: case.max_layer_width,
        rho0: case.rho0,
        fixed_extension_degree: case.r,
        ..WhirUniqueDecodingInputs::default()
    })?;
    if !selection.feasible {
        return Err(format!(
            "{}: WHIR_GR selector found no feasible parameters: {}",
            case.name,
            selection.notes.join("; ")
        )
        .into());
    }
    if selection.selected_r != case.r {
        return Err(format!(
            "{}: selector chose r={}, expected {}",
            case.name, selection.selected_r, case.r
        )
        .into());
    }
    if selection.public_params.initial_domain_size != case.n {
        return Err(format!(
            "{}: selector chose n={}, expected {}",
            case.name, selection.public_params.initial_domain_size, case.n
        )
        .into());
    }

    let ctx = Arc::new(GrContext::new(GrConfig {
        p: 2,
        k_exp: case.k_exp,
        r: selection.selected_r as usize,
    })?);
    let domain = Domain::teichmuller_subgroup(
        Arc::clone(&ctx),
        selection.public_params.initial_domain_size,
    )?;
    let omega = ctx.pow(domain.root(), u128::from(domain.size() / 3));
    let grid = ternary_grid(&ctx, &omega)?;
    let mut params =
        WhirGrPublicParameters::new(Arc::clone(&ctx), domain, case.variable_count, omega, grid);
    params.layer_widths = selection.public_params.layer_widths;
    params.shift_repetitions = selection.public_params.shift_repetitions;
    params.final_repetitions = selection.public_params.final_repetitions;
    params.degree_bounds = selection.public_params.degree_bounds;
    params.lambda_target = case.lambda_target;
    params.hash_id = BLAKE3;
    Ok(params)
}

fn multilinear_polynomial(
    ctx: &GrContext,
    variable_count: u64,
    seed: u64,
) -> Result<MultilinearPolynomial, Box<dyn Error>> {
    let coefficient_count = pow2_checked(variable_count)?;
    let coefficients = (0..coefficient_count)
        .map(|index| ctx.from_u64((seed.wrapping_add(13 * index).wrapping_add(7)) % 29))
        .collect();
    Ok(MultilinearPolynomial::new(variable_count, coefficients)?)
}

fn open_point(ctx: &GrContext, variable_count: u64, seed: u64) -> Vec<GrElem> {
    (0..variable_count)
        .map(|index| ctx.from_u64((seed.wrapping_add(7).wrapping_add(3 * index)) % 31))
        .collect()
}

fn main() {
    divan::main();
}
