use std::{error::Error, sync::Arc};

use crate::{
    algebra::galois_ring::{Domain, GrConfig, GrContext, GrElem},
    hash::BLAKE3,
    protocols::whir_gr::{
        common::{WhirGrCommitment, WhirGrOpening, WhirGrPublicParameters},
        constraint::ternary_grid,
        multiquadratic::{pow2_checked, MultilinearPolynomial},
        prover::{WhirGrCommitmentState, WhirGrProver},
        soundness::{
            select_whir_unique_decoding_parameters, WhirRational, WhirUniqueDecodingInputs,
        },
    },
};

#[doc(hidden)]
#[derive(Clone, Copy, Debug)]
pub struct WhirGrBenchCase {
    pub name: &'static str,
    pub k_exp: u32,
    pub r: u64,
    pub n: u64,
    pub variable_count: u64,
    pub max_layer_width: u64,
    pub lambda_target: u64,
    pub rho0: WhirRational,
}

#[doc(hidden)]
pub struct CommitInput {
    pub params: WhirGrPublicParameters,
    pub polynomial: MultilinearPolynomial,
}

#[doc(hidden)]
pub struct OpenInput {
    pub params: WhirGrPublicParameters,
    pub commitment: WhirGrCommitment,
    pub state: WhirGrCommitmentState,
    pub point: Vec<GrElem>,
}

#[doc(hidden)]
pub struct VerifyInput {
    pub params: WhirGrPublicParameters,
    pub commitment: WhirGrCommitment,
    pub point: Vec<GrElem>,
    pub opening: WhirGrOpening,
}

#[doc(hidden)]
pub const WHIR_GR_M4: WhirGrBenchCase = WhirGrBenchCase {
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
};

#[doc(hidden)]
pub const WHIR_GR_M5: WhirGrBenchCase = WhirGrBenchCase {
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
};

#[doc(hidden)]
pub const WHIR_GR_M6: WhirGrBenchCase = WhirGrBenchCase {
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
};

#[doc(hidden)]
pub const WHIR_GR_M7: WhirGrBenchCase = WhirGrBenchCase {
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
};

#[doc(hidden)]
pub const WHIR_GR_M8: WhirGrBenchCase = WhirGrBenchCase {
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
};

#[doc(hidden)]
pub const WHIR_GR_M9: WhirGrBenchCase = WhirGrBenchCase {
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
};

#[doc(hidden)]
pub const WHIR_GR_M10: WhirGrBenchCase = WhirGrBenchCase {
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
};

#[doc(hidden)]
pub const WHIR_GR_CASES: &[WhirGrBenchCase] = &[
    WHIR_GR_M4,
    WHIR_GR_M5,
    WHIR_GR_M6,
    WHIR_GR_M7,
    WHIR_GR_M8,
    WHIR_GR_M9,
    WHIR_GR_M10,
];

#[doc(hidden)]
pub const WHIR_GR_SMALL_CASES: &[WhirGrBenchCase] = &[WHIR_GR_M4, WHIR_GR_M5, WHIR_GR_M6];

impl WhirGrBenchCase {
    pub const fn short_name(self) -> &'static str {
        match self.variable_count {
            4 => "m4",
            5 => "m5",
            6 => "m6",
            7 => "m7",
            8 => "m8",
            9 => "m9",
            10 => "m10",
            _ => self.name,
        }
    }
}

#[doc(hidden)]
pub fn find_case(name: &str) -> Option<&'static WhirGrBenchCase> {
    WHIR_GR_CASES
        .iter()
        .find(|case| case.short_name() == name || case.name == name)
}

#[doc(hidden)]
pub fn commit_input(case: &WhirGrBenchCase) -> Result<CommitInput, Box<dyn Error>> {
    let params = build_params(case)?;
    let polynomial = multilinear_polynomial(&params.ctx, case.variable_count, 0)?;
    Ok(CommitInput { params, polynomial })
}

#[doc(hidden)]
pub fn open_input(case: &WhirGrBenchCase) -> Result<OpenInput, Box<dyn Error>> {
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

#[doc(hidden)]
pub fn verify_input(case: &WhirGrBenchCase) -> Result<VerifyInput, Box<dyn Error>> {
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

#[doc(hidden)]
pub fn build_params(case: &WhirGrBenchCase) -> Result<WhirGrPublicParameters, Box<dyn Error>> {
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

#[doc(hidden)]
pub fn multilinear_polynomial(
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

#[doc(hidden)]
pub fn open_point(ctx: &GrContext, variable_count: u64, seed: u64) -> Vec<GrElem> {
    (0..variable_count)
        .map(|index| ctx.from_u64((seed.wrapping_add(7).wrapping_add(3 * index)) % 31))
        .collect()
}
