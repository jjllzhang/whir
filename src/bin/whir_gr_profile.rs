use std::{error::Error, time::Instant};

use clap::{Parser, ValueEnum};
use serde_json::json;
use whir::protocols::whir_gr::{
    bench_support::{
        commit_bench_polynomial_profiled, commit_input, find_case, find_case_with_polynomial,
        open_input, verify_input, WhirGrBenchCase, WhirGrPolynomialKind,
    },
    prover::{WhirGrCommitTimings, WhirGrOpenTimings, WhirGrProver},
    serialization::serialize_opening,
    verifier::{WhirGrVerifier, WhirGrVerifyTimings},
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
    SummaryCsv,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum PolynomialKindArg {
    Multilinear,
    Multiquadratic,
}

impl From<PolynomialKindArg> for WhirGrPolynomialKind {
    fn from(value: PolynomialKindArg) -> Self {
        match value {
            PolynomialKindArg::Multilinear => Self::Multilinear,
            PolynomialKindArg::Multiquadratic => Self::MultiQuadratic,
        }
    }
}

#[derive(Debug, Parser)]
struct Args {
    #[arg(long = "case")]
    case_name: String,
    #[arg(long, value_enum)]
    polynomial: Option<PolynomialKindArg>,
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
    rayon_threads: usize,
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
    open_sumcheck_constraint_plan_ms: Option<f64>,
    open_sumcheck_poly_restrict_ms: Option<f64>,
    open_sumcheck_poly_eval_ms: Option<f64>,
    open_sumcheck_accumulate_ms: Option<f64>,
    open_sumcheck_interpolate_ms: Option<f64>,
    open_restrict_ms: Option<f64>,
    open_merkle_open_ms: Option<f64>,
    open_fold_indices_ms: Option<f64>,
    open_fold_eval_ms: Option<f64>,
    open_fold_shift_points_ms: Option<f64>,
    open_constraint_ms: Option<f64>,
    open_final_ms: Option<f64>,
    verify_algebra_ms: Option<f64>,
    verify_sumcheck_ms: Option<f64>,
    verify_merkle_ms: Option<f64>,
    verify_fold_ms: Option<f64>,
    verify_constraint_ms: Option<f64>,
    verify_final_ms: Option<f64>,
    serialized_opening_bytes: Option<usize>,
    accepted: Option<bool>,
}

#[derive(Default)]
struct RoundtripTotals {
    commit_ms: f64,
    open_ms: f64,
    verify_ms: f64,
    encode_oracle_ms: f64,
    merkle_ms: f64,
    to_multiquadratic_ms: f64,
    open_clone_ms: f64,
    open_init_ms: f64,
    open_sumcheck_ms: f64,
    open_sumcheck_constraint_plan_ms: f64,
    open_sumcheck_poly_restrict_ms: f64,
    open_sumcheck_poly_eval_ms: f64,
    open_sumcheck_accumulate_ms: f64,
    open_sumcheck_interpolate_ms: f64,
    open_restrict_ms: f64,
    open_fold_ms: f64,
    open_fold_indices_ms: f64,
    open_fold_eval_ms: f64,
    open_fold_shift_points_ms: f64,
    open_merkle_open_ms: f64,
    open_constraint_ms: f64,
    open_final_ms: f64,
    verify_sumcheck_ms: f64,
    verify_merkle_ms: f64,
    verify_fold_ms: f64,
    verify_constraint_ms: f64,
    verify_final_ms: f64,
    serialized_opening_bytes: Option<usize>,
    accepted: bool,
}

impl RoundtripTotals {
    fn add_commit(&mut self, elapsed_ms: f64, timings: &WhirGrCommitTimings) {
        self.commit_ms += elapsed_ms;
        self.encode_oracle_ms += timings.encode_oracle_ms;
        self.merkle_ms += timings.merkle_ms;
        self.to_multiquadratic_ms += timings.to_multiquadratic_ms;
    }

    fn add_open(&mut self, elapsed_ms: f64, serialized_bytes: usize, timings: &WhirGrOpenTimings) {
        self.open_ms += elapsed_ms;
        self.serialized_opening_bytes = Some(serialized_bytes);
        self.open_clone_ms += timings.clone_ms;
        self.open_init_ms += timings.init_ms;
        self.open_sumcheck_ms += timings.sumcheck_ms;
        self.open_sumcheck_constraint_plan_ms += timings.sumcheck_constraint_plan_ms;
        self.open_sumcheck_poly_restrict_ms += timings.sumcheck_poly_restrict_ms;
        self.open_sumcheck_poly_eval_ms += timings.sumcheck_poly_eval_ms;
        self.open_sumcheck_accumulate_ms += timings.sumcheck_accumulate_ms;
        self.open_sumcheck_interpolate_ms += timings.sumcheck_interpolate_ms;
        self.open_restrict_ms += timings.restrict_ms;
        self.open_fold_ms += timings.fold_ms;
        self.open_fold_indices_ms += timings.fold_indices_ms;
        self.open_fold_eval_ms += timings.fold_eval_ms;
        self.open_fold_shift_points_ms += timings.fold_shift_points_ms;
        self.open_merkle_open_ms += timings.merkle_open_ms;
        self.open_constraint_ms += timings.constraint_ms;
        self.open_final_ms += timings.final_ms;
    }

    fn add_verify(&mut self, elapsed_ms: f64, accepted: bool, timings: &WhirGrVerifyTimings) {
        self.verify_ms += elapsed_ms;
        self.accepted &= accepted;
        self.verify_sumcheck_ms += timings.sumcheck_ms;
        self.verify_merkle_ms += timings.merkle_ms;
        self.verify_fold_ms += timings.fold_ms;
        self.verify_constraint_ms += timings.constraint_ms;
        self.verify_final_ms += timings.final_ms;
    }

    fn into_row(self, reps: f64) -> ProfileRow {
        ProfileRow {
            rayon_threads: rayon_threads(),
            commit_ms: Some(self.commit_ms / reps),
            open_ms: Some(self.open_ms / reps),
            verify_ms: Some(self.verify_ms / reps),
            encode_oracle_ms: Some(self.encode_oracle_ms / reps),
            merkle_ms: Some(self.merkle_ms / reps),
            to_multiquadratic_ms: Some(self.to_multiquadratic_ms / reps),
            open_fold_ms: Some(self.open_fold_ms / reps),
            open_clone_ms: Some(self.open_clone_ms / reps),
            open_init_ms: Some(self.open_init_ms / reps),
            open_sumcheck_ms: Some(self.open_sumcheck_ms / reps),
            open_sumcheck_constraint_plan_ms: Some(self.open_sumcheck_constraint_plan_ms / reps),
            open_sumcheck_poly_restrict_ms: Some(self.open_sumcheck_poly_restrict_ms / reps),
            open_sumcheck_poly_eval_ms: Some(self.open_sumcheck_poly_eval_ms / reps),
            open_sumcheck_accumulate_ms: Some(self.open_sumcheck_accumulate_ms / reps),
            open_sumcheck_interpolate_ms: Some(self.open_sumcheck_interpolate_ms / reps),
            open_restrict_ms: Some(self.open_restrict_ms / reps),
            open_merkle_open_ms: Some(self.open_merkle_open_ms / reps),
            open_fold_indices_ms: Some(self.open_fold_indices_ms / reps),
            open_fold_eval_ms: Some(self.open_fold_eval_ms / reps),
            open_fold_shift_points_ms: Some(self.open_fold_shift_points_ms / reps),
            open_constraint_ms: Some(self.open_constraint_ms / reps),
            open_final_ms: Some(self.open_final_ms / reps),
            verify_algebra_ms: Some(self.verify_algebra_ms() / reps),
            verify_sumcheck_ms: Some(self.verify_sumcheck_ms / reps),
            verify_merkle_ms: Some(self.verify_merkle_ms / reps),
            verify_fold_ms: Some(self.verify_fold_ms / reps),
            verify_constraint_ms: Some(self.verify_constraint_ms / reps),
            verify_final_ms: Some(self.verify_final_ms / reps),
            serialized_opening_bytes: self.serialized_opening_bytes,
            accepted: Some(self.accepted),
        }
    }

    const fn verify_algebra_ms(&self) -> f64 {
        self.verify_sumcheck_ms
            + self.verify_fold_ms
            + self.verify_constraint_ms
            + self.verify_final_ms
    }
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

    let case = find_profile_case(&args)?;
    let row = match args.phase {
        ProfilePhase::Commit => profile_commit(case, args.reps)?,
        ProfilePhase::Open => profile_open(case, args.reps)?,
        ProfilePhase::Verify => profile_verify(case, args.reps)?,
        ProfilePhase::Roundtrip => profile_roundtrip(case, args.reps)?,
    };
    write_row(&args, case, &row);
    Ok(())
}

fn find_profile_case(args: &Args) -> Result<&'static WhirGrBenchCase, Box<dyn Error>> {
    if let Some(polynomial_kind) = args.polynomial {
        return find_case_with_polynomial(&args.case_name, polynomial_kind.into()).ok_or_else(
            || {
                format!(
                    "unknown WHIR_GR {:?} bench case '{}'",
                    polynomial_kind, args.case_name
                )
                .into()
            },
        );
    }

    find_case(&args.case_name)
        .ok_or_else(|| format!("unknown WHIR_GR bench case '{}'", args.case_name).into())
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
            commit_bench_polynomial_profiled(&prover, &input.polynomial)?;
        encode_oracle_ms += timings.encode_oracle_ms;
        merkle_ms += timings.merkle_ms;
        to_multiquadratic_ms += timings.to_multiquadratic_ms;
    }
    Ok(ProfileRow {
        rayon_threads: rayon_threads(),
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
    let mut sumcheck_constraint_plan_ms = 0.0;
    let mut sumcheck_poly_restrict_ms = 0.0;
    let mut sumcheck_poly_eval_ms = 0.0;
    let mut sumcheck_accumulate_ms = 0.0;
    let mut sumcheck_interpolate_ms = 0.0;
    let mut restrict_ms = 0.0;
    let mut encode_oracle_ms = 0.0;
    let mut merkle_ms = 0.0;
    let mut fold_ms = 0.0;
    let mut fold_indices_ms = 0.0;
    let mut fold_eval_ms = 0.0;
    let mut fold_shift_points_ms = 0.0;
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
        sumcheck_constraint_plan_ms += timings.sumcheck_constraint_plan_ms;
        sumcheck_poly_restrict_ms += timings.sumcheck_poly_restrict_ms;
        sumcheck_poly_eval_ms += timings.sumcheck_poly_eval_ms;
        sumcheck_accumulate_ms += timings.sumcheck_accumulate_ms;
        sumcheck_interpolate_ms += timings.sumcheck_interpolate_ms;
        restrict_ms += timings.restrict_ms;
        encode_oracle_ms += timings.encode_oracle_ms;
        merkle_ms += timings.merkle_ms;
        fold_ms += timings.fold_ms;
        fold_indices_ms += timings.fold_indices_ms;
        fold_eval_ms += timings.fold_eval_ms;
        fold_shift_points_ms += timings.fold_shift_points_ms;
        merkle_open_ms += timings.merkle_open_ms;
        constraint_ms += timings.constraint_ms;
        final_ms += timings.final_ms;
    }
    Ok(ProfileRow {
        rayon_threads: rayon_threads(),
        open_ms: Some(mean_ms(start, reps)),
        encode_oracle_ms: Some(encode_oracle_ms / reps as f64),
        merkle_ms: Some(merkle_ms / reps as f64),
        open_fold_ms: Some(fold_ms / reps as f64),
        open_clone_ms: Some(clone_ms / reps as f64),
        open_init_ms: Some(init_ms / reps as f64),
        open_sumcheck_ms: Some(sumcheck_ms / reps as f64),
        open_sumcheck_constraint_plan_ms: Some(sumcheck_constraint_plan_ms / reps as f64),
        open_sumcheck_poly_restrict_ms: Some(sumcheck_poly_restrict_ms / reps as f64),
        open_sumcheck_poly_eval_ms: Some(sumcheck_poly_eval_ms / reps as f64),
        open_sumcheck_accumulate_ms: Some(sumcheck_accumulate_ms / reps as f64),
        open_sumcheck_interpolate_ms: Some(sumcheck_interpolate_ms / reps as f64),
        open_restrict_ms: Some(restrict_ms / reps as f64),
        open_merkle_open_ms: Some(merkle_open_ms / reps as f64),
        open_fold_indices_ms: Some(fold_indices_ms / reps as f64),
        open_fold_eval_ms: Some(fold_eval_ms / reps as f64),
        open_fold_shift_points_ms: Some(fold_shift_points_ms / reps as f64),
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
    let mut sumcheck_ms = 0.0;
    let mut merkle_ms = 0.0;
    let mut fold_ms = 0.0;
    let mut constraint_ms = 0.0;
    let mut final_ms = 0.0;
    let start = Instant::now();
    for _ in 0..reps {
        let (verified, timings) =
            verifier.verify_profiled(&input.commitment, &input.point, &input.opening)?;
        accepted &= verified;
        sumcheck_ms += timings.sumcheck_ms;
        merkle_ms += timings.merkle_ms;
        fold_ms += timings.fold_ms;
        constraint_ms += timings.constraint_ms;
        final_ms += timings.final_ms;
    }
    let reps = reps as f64;
    Ok(ProfileRow {
        rayon_threads: rayon_threads(),
        verify_ms: Some(elapsed_ms(start) / reps),
        verify_algebra_ms: Some((sumcheck_ms + fold_ms + constraint_ms + final_ms) / reps),
        verify_sumcheck_ms: Some(sumcheck_ms / reps),
        verify_merkle_ms: Some(merkle_ms / reps),
        verify_fold_ms: Some(fold_ms / reps),
        verify_constraint_ms: Some(constraint_ms / reps),
        verify_final_ms: Some(final_ms / reps),
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
    let mut totals = RoundtripTotals {
        accepted: true,
        ..RoundtripTotals::default()
    };

    for _ in 0..reps {
        let start = Instant::now();
        let (commitment, state, timings) =
            commit_bench_polynomial_profiled(&prover, &input.polynomial)?;
        totals.add_commit(elapsed_ms(start), &timings);

        let start = Instant::now();
        let (opening, open_timings) = prover.open_profiled(&commitment, &state, &point)?;
        let serialized = serialize_opening(&input.params.ctx, &opening);
        totals.add_open(elapsed_ms(start), serialized.len(), &open_timings);

        let start = Instant::now();
        let (verified, verify_timings) = verifier.verify_profiled(&commitment, &point, &opening)?;
        totals.add_verify(elapsed_ms(start), verified, &verify_timings);
    }

    Ok(totals.into_row(reps as f64))
}

fn mean_ms(start: Instant, reps: u64) -> f64 {
    elapsed_ms(start) / reps as f64
}

fn elapsed_ms(start: Instant) -> f64 {
    start.elapsed().as_secs_f64() * 1000.0
}

fn rayon_threads() -> usize {
    #[cfg(feature = "parallel")]
    {
        rayon::current_num_threads()
    }
    #[cfg(not(feature = "parallel"))]
    {
        1
    }
}

fn write_row(args: &Args, case: &WhirGrBenchCase, row: &ProfileRow) {
    match args.format {
        OutputFormat::Csv => write_csv(args, case, row),
        OutputFormat::Text => write_text(args, case, row),
        OutputFormat::Json => write_json(args, case, row),
        OutputFormat::SummaryCsv => write_summary_csv(case, row),
    }
}

fn write_csv(args: &Args, case: &WhirGrBenchCase, row: &ProfileRow) {
    println!(
        "case,k_exp,r,n,variable_count,max_layer_width,lambda_target,rho0,phase,reps,rayon_threads,commit_ms,open_ms,verify_ms,encode_oracle_ms,merkle_ms,to_multiquadratic_ms,open_fold_ms,open_fold_indices_ms,open_fold_eval_ms,open_fold_shift_points_ms,open_clone_ms,open_init_ms,open_sumcheck_ms,open_sumcheck_constraint_plan_ms,open_sumcheck_poly_restrict_ms,open_sumcheck_poly_eval_ms,open_sumcheck_accumulate_ms,open_sumcheck_interpolate_ms,open_restrict_ms,open_merkle_open_ms,open_constraint_ms,open_final_ms,verify_algebra_ms,verify_sumcheck_ms,verify_merkle_ms,verify_fold_ms,verify_constraint_ms,verify_final_ms,serialized_opening_bytes,accepted"
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
        row.rayon_threads.to_string(),
        fmt_optional_f64(row.commit_ms),
        fmt_optional_f64(row.open_ms),
        fmt_optional_f64(row.verify_ms),
        fmt_optional_f64(row.encode_oracle_ms),
        fmt_optional_f64(row.merkle_ms),
        fmt_optional_f64(row.to_multiquadratic_ms),
        fmt_optional_f64(row.open_fold_ms),
        fmt_optional_f64(row.open_fold_indices_ms),
        fmt_optional_f64(row.open_fold_eval_ms),
        fmt_optional_f64(row.open_fold_shift_points_ms),
        fmt_optional_f64(row.open_clone_ms),
        fmt_optional_f64(row.open_init_ms),
        fmt_optional_f64(row.open_sumcheck_ms),
        fmt_optional_f64(row.open_sumcheck_constraint_plan_ms),
        fmt_optional_f64(row.open_sumcheck_poly_restrict_ms),
        fmt_optional_f64(row.open_sumcheck_poly_eval_ms),
        fmt_optional_f64(row.open_sumcheck_accumulate_ms),
        fmt_optional_f64(row.open_sumcheck_interpolate_ms),
        fmt_optional_f64(row.open_restrict_ms),
        fmt_optional_f64(row.open_merkle_open_ms),
        fmt_optional_f64(row.open_constraint_ms),
        fmt_optional_f64(row.open_final_ms),
        fmt_optional_f64(row.verify_algebra_ms),
        fmt_optional_f64(row.verify_sumcheck_ms),
        fmt_optional_f64(row.verify_merkle_ms),
        fmt_optional_f64(row.verify_fold_ms),
        fmt_optional_f64(row.verify_constraint_ms),
        fmt_optional_f64(row.verify_final_ms),
        fmt_optional_usize(row.serialized_opening_bytes),
        fmt_optional_bool(row.accepted),
    ];
    println!("{}", fields.join(","));
}

fn write_summary_csv(case: &WhirGrBenchCase, row: &ProfileRow) {
    println!(
        "case,ring,k_exp,r,lambda_target,variable_count,max_layer_width,poly_dim,message_length,n_0,rho,commit_ms,open_ms,prove_ms,verify_ms,proof_size_bytes,proof_size_kb"
    );
    let message_length = pow_u64(3, case.variable_count);
    let proof_size_kb = row
        .serialized_opening_bytes
        .map(|bytes| format!("{:.3}", bytes as f64 / 1024.0))
        .unwrap_or_default();
    let fields = [
        case.name.to_string(),
        format!("\"GR(2^{},{})\"", case.k_exp, case.r),
        case.k_exp.to_string(),
        case.r.to_string(),
        case.lambda_target.to_string(),
        case.variable_count.to_string(),
        case.max_layer_width.to_string(),
        source_coefficient_count(case).to_string(),
        message_length.to_string(),
        case.n.to_string(),
        reduced_fraction(message_length, case.n),
        fmt_optional_f64(row.commit_ms),
        fmt_optional_f64(row.open_ms),
        fmt_optional_f64(sum_options(row.commit_ms, row.open_ms)),
        fmt_optional_f64(row.verify_ms),
        fmt_optional_usize(row.serialized_opening_bytes),
        proof_size_kb,
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
    println!("rayon_threads={}", row.rayon_threads);
    print_optional_f64("commit_ms", row.commit_ms);
    print_optional_f64("open_ms", row.open_ms);
    print_optional_f64("verify_ms", row.verify_ms);
    print_optional_f64("encode_oracle_ms", row.encode_oracle_ms);
    print_optional_f64("merkle_ms", row.merkle_ms);
    print_optional_f64("to_multiquadratic_ms", row.to_multiquadratic_ms);
    print_optional_f64("open_fold_ms", row.open_fold_ms);
    print_optional_f64("open_fold_indices_ms", row.open_fold_indices_ms);
    print_optional_f64("open_fold_eval_ms", row.open_fold_eval_ms);
    print_optional_f64("open_fold_shift_points_ms", row.open_fold_shift_points_ms);
    print_optional_f64("open_clone_ms", row.open_clone_ms);
    print_optional_f64("open_init_ms", row.open_init_ms);
    print_optional_f64("open_sumcheck_ms", row.open_sumcheck_ms);
    print_optional_f64(
        "open_sumcheck_constraint_plan_ms",
        row.open_sumcheck_constraint_plan_ms,
    );
    print_optional_f64(
        "open_sumcheck_poly_restrict_ms",
        row.open_sumcheck_poly_restrict_ms,
    );
    print_optional_f64("open_sumcheck_poly_eval_ms", row.open_sumcheck_poly_eval_ms);
    print_optional_f64(
        "open_sumcheck_accumulate_ms",
        row.open_sumcheck_accumulate_ms,
    );
    print_optional_f64(
        "open_sumcheck_interpolate_ms",
        row.open_sumcheck_interpolate_ms,
    );
    print_optional_f64("open_restrict_ms", row.open_restrict_ms);
    print_optional_f64("open_merkle_open_ms", row.open_merkle_open_ms);
    print_optional_f64("open_constraint_ms", row.open_constraint_ms);
    print_optional_f64("open_final_ms", row.open_final_ms);
    print_optional_f64("verify_algebra_ms", row.verify_algebra_ms);
    print_optional_f64("verify_sumcheck_ms", row.verify_sumcheck_ms);
    print_optional_f64("verify_merkle_ms", row.verify_merkle_ms);
    print_optional_f64("verify_fold_ms", row.verify_fold_ms);
    print_optional_f64("verify_constraint_ms", row.verify_constraint_ms);
    print_optional_f64("verify_final_ms", row.verify_final_ms);
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
            "rayon_threads": row.rayon_threads,
            "commit_ms": row.commit_ms,
            "open_ms": row.open_ms,
            "verify_ms": row.verify_ms,
            "encode_oracle_ms": row.encode_oracle_ms,
            "merkle_ms": row.merkle_ms,
            "to_multiquadratic_ms": row.to_multiquadratic_ms,
            "open_fold_ms": row.open_fold_ms,
            "open_fold_indices_ms": row.open_fold_indices_ms,
            "open_fold_eval_ms": row.open_fold_eval_ms,
            "open_fold_shift_points_ms": row.open_fold_shift_points_ms,
            "open_clone_ms": row.open_clone_ms,
            "open_init_ms": row.open_init_ms,
            "open_sumcheck_ms": row.open_sumcheck_ms,
            "open_sumcheck_constraint_plan_ms": row.open_sumcheck_constraint_plan_ms,
            "open_sumcheck_poly_restrict_ms": row.open_sumcheck_poly_restrict_ms,
            "open_sumcheck_poly_eval_ms": row.open_sumcheck_poly_eval_ms,
            "open_sumcheck_accumulate_ms": row.open_sumcheck_accumulate_ms,
            "open_sumcheck_interpolate_ms": row.open_sumcheck_interpolate_ms,
            "open_restrict_ms": row.open_restrict_ms,
            "open_merkle_open_ms": row.open_merkle_open_ms,
            "open_constraint_ms": row.open_constraint_ms,
            "open_final_ms": row.open_final_ms,
            "verify_algebra_ms": row.verify_algebra_ms,
            "verify_sumcheck_ms": row.verify_sumcheck_ms,
            "verify_merkle_ms": row.verify_merkle_ms,
            "verify_fold_ms": row.verify_fold_ms,
            "verify_constraint_ms": row.verify_constraint_ms,
            "verify_final_ms": row.verify_final_ms,
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

fn source_coefficient_count(case: &WhirGrBenchCase) -> u64 {
    match case.polynomial_kind {
        WhirGrPolynomialKind::Multilinear => pow_u64(2, case.variable_count),
        WhirGrPolynomialKind::MultiQuadratic => pow_u64(3, case.variable_count),
    }
}

fn pow_u64(base: u64, exponent: u64) -> u64 {
    let mut out = 1;
    for _ in 0..exponent {
        out *= base;
    }
    out
}

fn reduced_fraction(numerator: u64, denominator: u64) -> String {
    let divisor = gcd(numerator, denominator);
    format!("{}/{}", numerator / divisor, denominator / divisor)
}

const fn gcd(mut left: u64, mut right: u64) -> u64 {
    while right != 0 {
        let next = left % right;
        left = right;
        right = next;
    }
    left
}

fn sum_options(left: Option<f64>, right: Option<f64>) -> Option<f64> {
    Some(left? + right?)
}
