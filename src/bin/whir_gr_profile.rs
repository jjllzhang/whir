use std::{error::Error, time::Instant};

use clap::{Parser, ValueEnum};
use serde_json::json;
use whir::protocols::whir_gr::{
    bench_support::{commit_input, find_case, open_input, verify_input, WhirGrBenchCase},
    prover::WhirGrProver,
    serialization::serialize_opening,
    verifier::WhirGrVerifier,
};

#[derive(Clone, Copy, Debug, ValueEnum)]
enum ProfilePhase {
    Commit,
    Open,
    Verify,
    Roundtrip,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum OutputFormat {
    Csv,
    Text,
    Json,
}

#[derive(Debug, Parser)]
struct Args {
    #[arg(long = "case")]
    case_name: String,
    #[arg(long, value_enum)]
    phase: ProfilePhase,
    #[arg(long, default_value_t = 1)]
    reps: u64,
    #[arg(long, default_value_t = false)]
    allocator_stats: bool,
    #[arg(long, value_enum, default_value_t = OutputFormat::Text)]
    format: OutputFormat,
}

#[derive(Default)]
struct ProfileRow {
    commit_ms: Option<f64>,
    open_ms: Option<f64>,
    verify_ms: Option<f64>,
    encode_oracle_ms: Option<f64>,
    merkle_ms: Option<f64>,
    to_multiquadratic_ms: Option<f64>,
    open_fold_ms: Option<f64>,
    open_clone_ms: Option<f64>,
    open_init_ms: Option<f64>,
    open_sumcheck_ms: Option<f64>,
    open_restrict_ms: Option<f64>,
    open_merkle_open_ms: Option<f64>,
    open_constraint_ms: Option<f64>,
    open_final_ms: Option<f64>,
    verify_algebra_ms: Option<f64>,
    serialized_opening_bytes: Option<usize>,
    accepted: Option<bool>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    if args.reps == 0 {
        return Err("--reps must be greater than zero".into());
    }
    if args.allocator_stats {
        eprintln!(
            "--allocator-stats is accepted for CLI stability; use the Divan bench for allocator tallies"
        );
    }

    let case = find_case(&args.case_name)
        .ok_or_else(|| format!("unknown WHIR_GR bench case '{}'", args.case_name))?;
    let row = match args.phase {
        ProfilePhase::Commit => profile_commit(case, args.reps)?,
        ProfilePhase::Open => profile_open(case, args.reps)?,
        ProfilePhase::Verify => profile_verify(case, args.reps)?,
        ProfilePhase::Roundtrip => profile_roundtrip(case, args.reps)?,
    };
    write_row(&args, case, &row);
    Ok(())
}

fn profile_commit(case: &WhirGrBenchCase, reps: u64) -> Result<ProfileRow, Box<dyn Error>> {
    let input = commit_input(case)?;
    let prover = WhirGrProver::new(&input.params);
    let mut encode_oracle_ms = 0.0;
    let mut merkle_ms = 0.0;
    let mut to_multiquadratic_ms = 0.0;
    let start = Instant::now();
    for _ in 0..reps {
        let (_commitment, _state, timings) =
            prover.commit_multilinear_profiled(&input.polynomial)?;
        encode_oracle_ms += timings.encode_oracle_ms;
        merkle_ms += timings.merkle_ms;
        to_multiquadratic_ms += timings.to_multiquadratic_ms;
    }
    Ok(ProfileRow {
        commit_ms: Some(mean_ms(start, reps)),
        encode_oracle_ms: Some(encode_oracle_ms / reps as f64),
        merkle_ms: Some(merkle_ms / reps as f64),
        to_multiquadratic_ms: Some(to_multiquadratic_ms / reps as f64),
        ..ProfileRow::default()
    })
}

fn profile_open(case: &WhirGrBenchCase, reps: u64) -> Result<ProfileRow, Box<dyn Error>> {
    let input = open_input(case)?;
    let prover = WhirGrProver::new(&input.params);
    let mut serialized_opening_bytes = None;
    let mut clone_ms = 0.0;
    let mut init_ms = 0.0;
    let mut sumcheck_ms = 0.0;
    let mut restrict_ms = 0.0;
    let mut encode_oracle_ms = 0.0;
    let mut merkle_ms = 0.0;
    let mut fold_ms = 0.0;
    let mut merkle_open_ms = 0.0;
    let mut constraint_ms = 0.0;
    let mut final_ms = 0.0;
    let start = Instant::now();
    for _ in 0..reps {
        let (opening, timings) =
            prover.open_profiled(&input.commitment, &input.state, &input.point)?;
        serialized_opening_bytes = Some(serialize_opening(&input.params.ctx, &opening).len());
        clone_ms += timings.clone_ms;
        init_ms += timings.init_ms;
        sumcheck_ms += timings.sumcheck_ms;
        restrict_ms += timings.restrict_ms;
        encode_oracle_ms += timings.encode_oracle_ms;
        merkle_ms += timings.merkle_ms;
        fold_ms += timings.fold_ms;
        merkle_open_ms += timings.merkle_open_ms;
        constraint_ms += timings.constraint_ms;
        final_ms += timings.final_ms;
    }
    Ok(ProfileRow {
        open_ms: Some(mean_ms(start, reps)),
        encode_oracle_ms: Some(encode_oracle_ms / reps as f64),
        merkle_ms: Some(merkle_ms / reps as f64),
        open_fold_ms: Some(fold_ms / reps as f64),
        open_clone_ms: Some(clone_ms / reps as f64),
        open_init_ms: Some(init_ms / reps as f64),
        open_sumcheck_ms: Some(sumcheck_ms / reps as f64),
        open_restrict_ms: Some(restrict_ms / reps as f64),
        open_merkle_open_ms: Some(merkle_open_ms / reps as f64),
        open_constraint_ms: Some(constraint_ms / reps as f64),
        open_final_ms: Some(final_ms / reps as f64),
        serialized_opening_bytes,
        ..ProfileRow::default()
    })
}

fn profile_verify(case: &WhirGrBenchCase, reps: u64) -> Result<ProfileRow, Box<dyn Error>> {
    let input = verify_input(case)?;
    let verifier = WhirGrVerifier::new(&input.params);
    let serialized_opening_bytes = serialize_opening(&input.params.ctx, &input.opening).len();
    let mut accepted = true;
    let start = Instant::now();
    for _ in 0..reps {
        accepted &= verifier.verify(&input.commitment, &input.point, &input.opening)?;
    }
    Ok(ProfileRow {
        verify_ms: Some(mean_ms(start, reps)),
        serialized_opening_bytes: Some(serialized_opening_bytes),
        accepted: Some(accepted),
        ..ProfileRow::default()
    })
}

fn profile_roundtrip(case: &WhirGrBenchCase, reps: u64) -> Result<ProfileRow, Box<dyn Error>> {
    let input = commit_input(case)?;
    let point = whir::protocols::whir_gr::bench_support::open_point(
        &input.params.ctx,
        case.variable_count,
        0,
    );
    let prover = WhirGrProver::new(&input.params);
    let verifier = WhirGrVerifier::new(&input.params);
    let mut commit_ms = 0.0;
    let mut open_ms = 0.0;
    let mut verify_ms = 0.0;
    let mut row_encode_oracle_ms = 0.0;
    let mut row_merkle_ms = 0.0;
    let mut row_to_multiquadratic_ms = 0.0;
    let mut row_open_clone_ms = 0.0;
    let mut row_open_init_ms = 0.0;
    let mut row_open_sumcheck_ms = 0.0;
    let mut row_open_restrict_ms = 0.0;
    let mut row_open_fold_ms = 0.0;
    let mut row_open_merkle_open_ms = 0.0;
    let mut row_open_constraint_ms = 0.0;
    let mut row_open_final_ms = 0.0;
    let mut serialized_opening_bytes = None;
    let mut accepted = true;

    for _ in 0..reps {
        let start = Instant::now();
        let (commitment, state, timings) = prover.commit_multilinear_profiled(&input.polynomial)?;
        commit_ms += elapsed_ms(start);
        let encode_oracle_ms = timings.encode_oracle_ms;
        let merkle_ms = timings.merkle_ms;
        let to_multiquadratic_ms = timings.to_multiquadratic_ms;

        let start = Instant::now();
        let (opening, open_timings) = prover.open_profiled(&commitment, &state, &point)?;
        let serialized = serialize_opening(&input.params.ctx, &opening);
        serialized_opening_bytes = Some(serialized.len());
        open_ms += elapsed_ms(start);

        let start = Instant::now();
        accepted &= verifier.verify(&commitment, &point, &opening)?;
        verify_ms += elapsed_ms(start);

        row_encode_oracle_ms += encode_oracle_ms;
        row_merkle_ms += merkle_ms;
        row_to_multiquadratic_ms += to_multiquadratic_ms;

        row_open_clone_ms += open_timings.clone_ms;
        row_open_init_ms += open_timings.init_ms;
        row_open_sumcheck_ms += open_timings.sumcheck_ms;
        row_open_restrict_ms += open_timings.restrict_ms;
        row_open_fold_ms += open_timings.fold_ms;
        row_open_merkle_open_ms += open_timings.merkle_open_ms;
        row_open_constraint_ms += open_timings.constraint_ms;
        row_open_final_ms += open_timings.final_ms;
    }

    let reps = reps as f64;
    Ok(ProfileRow {
        commit_ms: Some(commit_ms / reps),
        open_ms: Some(open_ms / reps),
        verify_ms: Some(verify_ms / reps),
        encode_oracle_ms: Some(row_encode_oracle_ms / reps),
        merkle_ms: Some(row_merkle_ms / reps),
        to_multiquadratic_ms: Some(row_to_multiquadratic_ms / reps),
        open_fold_ms: Some(row_open_fold_ms / reps),
        open_clone_ms: Some(row_open_clone_ms / reps),
        open_init_ms: Some(row_open_init_ms / reps),
        open_sumcheck_ms: Some(row_open_sumcheck_ms / reps),
        open_restrict_ms: Some(row_open_restrict_ms / reps),
        open_merkle_open_ms: Some(row_open_merkle_open_ms / reps),
        open_constraint_ms: Some(row_open_constraint_ms / reps),
        open_final_ms: Some(row_open_final_ms / reps),
        serialized_opening_bytes,
        accepted: Some(accepted),
        ..ProfileRow::default()
    })
}

fn mean_ms(start: Instant, reps: u64) -> f64 {
    elapsed_ms(start) / reps as f64
}

fn elapsed_ms(start: Instant) -> f64 {
    start.elapsed().as_secs_f64() * 1000.0
}

fn write_row(args: &Args, case: &WhirGrBenchCase, row: &ProfileRow) {
    match args.format {
        OutputFormat::Csv => write_csv(args, case, row),
        OutputFormat::Text => write_text(args, case, row),
        OutputFormat::Json => write_json(args, case, row),
    }
}

fn write_csv(args: &Args, case: &WhirGrBenchCase, row: &ProfileRow) {
    println!(
        "case,k_exp,r,n,variable_count,max_layer_width,lambda_target,rho0,phase,reps,commit_ms,open_ms,verify_ms,encode_oracle_ms,merkle_ms,to_multiquadratic_ms,open_fold_ms,open_clone_ms,open_init_ms,open_sumcheck_ms,open_restrict_ms,open_merkle_open_ms,open_constraint_ms,open_final_ms,verify_algebra_ms,serialized_opening_bytes,accepted"
    );
    let fields = [
        case.name.to_string(),
        case.k_exp.to_string(),
        case.r.to_string(),
        case.n.to_string(),
        case.variable_count.to_string(),
        case.max_layer_width.to_string(),
        case.lambda_target.to_string(),
        format!("{}/{}", case.rho0.numerator, case.rho0.denominator),
        format!("{:?}", args.phase),
        args.reps.to_string(),
        fmt_optional_f64(row.commit_ms),
        fmt_optional_f64(row.open_ms),
        fmt_optional_f64(row.verify_ms),
        fmt_optional_f64(row.encode_oracle_ms),
        fmt_optional_f64(row.merkle_ms),
        fmt_optional_f64(row.to_multiquadratic_ms),
        fmt_optional_f64(row.open_fold_ms),
        fmt_optional_f64(row.open_clone_ms),
        fmt_optional_f64(row.open_init_ms),
        fmt_optional_f64(row.open_sumcheck_ms),
        fmt_optional_f64(row.open_restrict_ms),
        fmt_optional_f64(row.open_merkle_open_ms),
        fmt_optional_f64(row.open_constraint_ms),
        fmt_optional_f64(row.open_final_ms),
        fmt_optional_f64(row.verify_algebra_ms),
        fmt_optional_usize(row.serialized_opening_bytes),
        fmt_optional_bool(row.accepted),
    ];
    println!("{}", fields.join(","));
}

fn write_text(args: &Args, case: &WhirGrBenchCase, row: &ProfileRow) {
    println!("case={}", case.name);
    println!("k_exp={}", case.k_exp);
    println!("r={}", case.r);
    println!("n={}", case.n);
    println!("variable_count={}", case.variable_count);
    println!("max_layer_width={}", case.max_layer_width);
    println!("lambda_target={}", case.lambda_target);
    println!("rho0={}/{}", case.rho0.numerator, case.rho0.denominator);
    println!("phase={:?}", args.phase);
    println!("reps={}", args.reps);
    print_optional_f64("commit_ms", row.commit_ms);
    print_optional_f64("open_ms", row.open_ms);
    print_optional_f64("verify_ms", row.verify_ms);
    print_optional_f64("encode_oracle_ms", row.encode_oracle_ms);
    print_optional_f64("merkle_ms", row.merkle_ms);
    print_optional_f64("to_multiquadratic_ms", row.to_multiquadratic_ms);
    print_optional_f64("open_fold_ms", row.open_fold_ms);
    print_optional_f64("open_clone_ms", row.open_clone_ms);
    print_optional_f64("open_init_ms", row.open_init_ms);
    print_optional_f64("open_sumcheck_ms", row.open_sumcheck_ms);
    print_optional_f64("open_restrict_ms", row.open_restrict_ms);
    print_optional_f64("open_merkle_open_ms", row.open_merkle_open_ms);
    print_optional_f64("open_constraint_ms", row.open_constraint_ms);
    print_optional_f64("open_final_ms", row.open_final_ms);
    print_optional_f64("verify_algebra_ms", row.verify_algebra_ms);
    println!(
        "serialized_opening_bytes={}",
        fmt_optional_usize(row.serialized_opening_bytes)
    );
    println!("accepted={}", fmt_optional_bool(row.accepted));
}

fn write_json(args: &Args, case: &WhirGrBenchCase, row: &ProfileRow) {
    println!(
        "{}",
        json!({
            "case": case.name,
            "k_exp": case.k_exp,
            "r": case.r,
            "n": case.n,
            "variable_count": case.variable_count,
            "max_layer_width": case.max_layer_width,
            "lambda_target": case.lambda_target,
            "rho0": {
                "numerator": case.rho0.numerator,
                "denominator": case.rho0.denominator,
            },
            "phase": format!("{:?}", args.phase),
            "reps": args.reps,
            "commit_ms": row.commit_ms,
            "open_ms": row.open_ms,
            "verify_ms": row.verify_ms,
            "encode_oracle_ms": row.encode_oracle_ms,
            "merkle_ms": row.merkle_ms,
            "to_multiquadratic_ms": row.to_multiquadratic_ms,
            "open_fold_ms": row.open_fold_ms,
            "open_clone_ms": row.open_clone_ms,
            "open_init_ms": row.open_init_ms,
            "open_sumcheck_ms": row.open_sumcheck_ms,
            "open_restrict_ms": row.open_restrict_ms,
            "open_merkle_open_ms": row.open_merkle_open_ms,
            "open_constraint_ms": row.open_constraint_ms,
            "open_final_ms": row.open_final_ms,
            "verify_algebra_ms": row.verify_algebra_ms,
            "serialized_opening_bytes": row.serialized_opening_bytes,
            "accepted": row.accepted,
        })
    );
}

fn print_optional_f64(label: &str, value: Option<f64>) {
    println!("{label}={}", fmt_optional_f64(value));
}

fn fmt_optional_f64(value: Option<f64>) -> String {
    value.map(|value| format!("{value:.6}")).unwrap_or_default()
}

fn fmt_optional_usize(value: Option<usize>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn fmt_optional_bool(value: Option<bool>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}
