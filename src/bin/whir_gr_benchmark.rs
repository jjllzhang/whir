use std::{error::Error, sync::Arc, time::Instant};

use clap::{Parser, ValueEnum};
use whir::{
    algebra::galois_ring::{Domain, GrConfig, GrContext},
    protocols::whir_gr::{
        common::WhirGrPublicParameters,
        constraint::ternary_grid,
        multiquadratic::{
            pow2_checked, pow3_checked, MultiQuadraticPolynomial, MultilinearPolynomial,
        },
        prover::WhirGrProver,
        serialization::serialize_opening,
        soundness::{
            select_whir_unique_decoding_parameters, WhirRational, WhirUniqueDecodingInputs,
        },
        verifier::WhirGrVerifier,
    },
};

#[derive(Clone, Copy, Debug, ValueEnum)]
enum PolynomialKind {
    Multiquadratic,
    Multilinear,
}

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "Benchmark the Rust WHIR_GR unique-decoding prototype"
)]
struct Args {
    #[arg(long, default_value_t = 2)]
    p: u64,

    #[arg(long = "k-exp", default_value_t = 16)]
    k_exp: u32,

    #[arg(long)]
    r: Option<u64>,

    #[arg(long)]
    n: Option<u64>,

    #[arg(long = "m", default_value_t = 3)]
    variable_count: u64,

    #[arg(long = "bmax", default_value_t = 1)]
    max_layer_width: u64,

    #[arg(long = "lambda", default_value_t = 32)]
    lambda_target: u64,

    #[arg(long = "rho0", default_value = "1/3")]
    rho0: String,

    #[arg(long = "repetitions", default_value_t = 1)]
    repetitions: u64,

    #[arg(long = "seed", default_value_t = 0)]
    seed: u64,

    #[arg(long = "polynomial", value_enum, default_value_t = PolynomialKind::Multilinear)]
    polynomial_kind: PolynomialKind,

    #[arg(long = "max-extension-degree", default_value_t = 0)]
    max_extension_degree: u64,

    #[arg(long = "max-domain-size", default_value_t = 0)]
    max_domain_size: u64,

    #[arg(long = "csv-header", default_value_t = false)]
    csv_header: bool,
}

struct BenchmarkCase {
    params: WhirGrPublicParameters,
    selection_effective_security_bits: u64,
    rate: String,
    polynomial_kind: PolynomialKind,
    seed: u64,
}

struct TimingTotals {
    commit_ms: f64,
    open_ms: f64,
    verify_ms: f64,
    serialized_bytes_actual: usize,
}

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    if args.p != 2 {
        return Err("WHIR_GR prototype currently requires --p 2".into());
    }
    if args.repetitions == 0 {
        return Err("--repetitions must be greater than zero".into());
    }

    let case = build_case(&args)?;
    let totals = run_benchmark(&case, args.repetitions)?;
    if args.csv_header {
        println!(
            "protocol,p,k_exp,r,n,rate,lambda,effective_security_bits,commit_ms,open_ms,verify_ms,serialized_bytes_actual,polynomial,seed,repetitions"
        );
    }
    println!(
        "whir_gr_ud,{},{},{},{},{},{},{},{:.6},{:.6},{:.6},{},{:?},{},{}",
        args.p,
        args.k_exp,
        case.params.ctx.config().r,
        case.params.initial_domain.size(),
        csv_escape(&case.rate),
        args.lambda_target,
        case.selection_effective_security_bits,
        totals.commit_ms / args.repetitions as f64,
        totals.open_ms / args.repetitions as f64,
        totals.verify_ms / args.repetitions as f64,
        totals.serialized_bytes_actual,
        case.polynomial_kind,
        case.seed,
        args.repetitions
    );
    Ok(())
}

fn build_case(args: &Args) -> Result<BenchmarkCase, Box<dyn Error>> {
    let rho0 = parse_rational(&args.rho0)?;
    let selection = select_whir_unique_decoding_parameters(&WhirUniqueDecodingInputs {
        lambda_target: args.lambda_target,
        ring_exponent: u64::from(args.k_exp),
        variable_count: args.variable_count,
        max_layer_width: args.max_layer_width,
        rho0,
        fixed_extension_degree: args.r.unwrap_or(0),
        max_extension_degree: args.max_extension_degree,
        max_domain_size: args.max_domain_size,
        ..WhirUniqueDecodingInputs::default()
    })?;
    if !selection.feasible {
        return Err(format!(
            "WHIR_GR selector found no feasible parameters: {}",
            selection.notes.join("; ")
        )
        .into());
    }
    if let Some(expected_n) = args.n {
        if selection.public_params.initial_domain_size != expected_n {
            return Err(format!(
                "--n mismatch: selector chose {}, requested {expected_n}",
                selection.public_params.initial_domain_size
            )
            .into());
        }
    }

    let ctx = Arc::new(GrContext::new(GrConfig {
        p: 2,
        k_exp: args.k_exp,
        r: selection.selected_r as usize,
    })?);
    let domain = Domain::teichmuller_subgroup(
        Arc::clone(&ctx),
        selection.public_params.initial_domain_size,
    )?;
    let omega = ctx.pow(domain.root(), u128::from(domain.size() / 3));
    let grid = ternary_grid(&ctx, &omega)?;
    let mut params =
        WhirGrPublicParameters::new(Arc::clone(&ctx), domain, args.variable_count, omega, grid);
    params.layer_widths = selection.public_params.layer_widths;
    params.shift_repetitions = selection.public_params.shift_repetitions;
    params.final_repetitions = selection.public_params.final_repetitions;
    params.degree_bounds = selection.public_params.degree_bounds;
    params.lambda_target = args.lambda_target;

    let rate = reduced_ratio(
        pow3_checked(args.variable_count)?,
        params.initial_domain.size(),
    );
    Ok(BenchmarkCase {
        params,
        selection_effective_security_bits: selection.effective_security_bits,
        rate,
        polynomial_kind: args.polynomial_kind,
        seed: args.seed,
    })
}

fn run_benchmark(case: &BenchmarkCase, repetitions: u64) -> Result<TimingTotals, Box<dyn Error>> {
    let ctx = &case.params.ctx;
    let prover = WhirGrProver::new(&case.params);
    let verifier = WhirGrVerifier::new(&case.params);
    let point = open_point(ctx, case.params.variable_count, case.seed);
    let mut totals = TimingTotals {
        commit_ms: 0.0,
        open_ms: 0.0,
        verify_ms: 0.0,
        serialized_bytes_actual: 0,
    };

    for _ in 0..repetitions {
        let commit_start = Instant::now();
        let (commitment, state) = match case.polynomial_kind {
            PolynomialKind::Multiquadratic => {
                let polynomial =
                    multiquadratic_polynomial(ctx, case.params.variable_count, case.seed)?;
                prover.commit(&polynomial)?
            }
            PolynomialKind::Multilinear => {
                let polynomial =
                    multilinear_polynomial(ctx, case.params.variable_count, case.seed)?;
                prover.commit_multilinear(&polynomial)?
            }
        };
        totals.commit_ms += commit_start.elapsed().as_secs_f64() * 1_000.0;

        let open_start = Instant::now();
        let opening = prover.open(&commitment, &state, &point)?;
        totals.open_ms += open_start.elapsed().as_secs_f64() * 1_000.0;

        let verify_start = Instant::now();
        if !verifier.verify(&commitment, &point, &opening)? {
            return Err("WHIR_GR verifier rejected an honest benchmark proof".into());
        }
        totals.verify_ms += verify_start.elapsed().as_secs_f64() * 1_000.0;
        totals.serialized_bytes_actual = serialize_opening(ctx, &opening).len();
    }
    Ok(totals)
}

fn multiquadratic_polynomial(
    ctx: &GrContext,
    variable_count: u64,
    seed: u64,
) -> Result<MultiQuadraticPolynomial, Box<dyn Error>> {
    let coefficient_count = pow3_checked(variable_count)?;
    let coefficients = (0..coefficient_count)
        .map(|index| ctx.from_u64((seed.wrapping_add(11 * index).wrapping_add(5)) % 23))
        .collect();
    Ok(MultiQuadraticPolynomial::new(variable_count, coefficients)?)
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

fn open_point(
    ctx: &GrContext,
    variable_count: u64,
    seed: u64,
) -> Vec<whir::algebra::galois_ring::GrElem> {
    (0..variable_count)
        .map(|index| ctx.from_u64((seed.wrapping_add(7).wrapping_add(3 * index)) % 31))
        .collect()
}

fn parse_rational(value: &str) -> Result<WhirRational, Box<dyn Error>> {
    let (numerator, denominator) = value
        .split_once('/')
        .ok_or("--rho0 must have the form numerator/denominator")?;
    Ok(WhirRational {
        numerator: numerator.parse()?,
        denominator: denominator.parse()?,
    })
}

fn reduced_ratio(numerator: u64, denominator: u64) -> String {
    let divisor = gcd(numerator, denominator);
    format!("{}/{}", numerator / divisor, denominator / divisor)
}

const fn gcd(mut lhs: u64, mut rhs: u64) -> u64 {
    while rhs != 0 {
        let remainder = lhs % rhs;
        lhs = rhs;
        rhs = remainder;
    }
    lhs
}

fn csv_escape(value: &str) -> String {
    if value.contains([',', '"', '\n']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_owned()
    }
}
