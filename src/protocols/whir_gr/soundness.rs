use crate::algebra::galois_ring::{GrError, Result};

const WHIR_SUMCHECK_DEGREE_BOUND: u64 = 4;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WhirRational {
    pub numerator: u64,
    pub denominator: u64,
}

impl Default for WhirRational {
    fn default() -> Self {
        Self {
            numerator: 1,
            denominator: 3,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WhirUniqueDecodingInputs {
    pub lambda_target: u64,
    pub ring_exponent: u64,
    pub variable_count: u64,
    pub max_layer_width: u64,
    pub rho0: WhirRational,
    pub fixed_extension_degree: u64,
    pub max_extension_degree: u64,
    pub max_domain_size: u64,
    pub max_n0_search_steps: u64,
}

impl Default for WhirUniqueDecodingInputs {
    fn default() -> Self {
        Self {
            lambda_target: 128,
            ring_exponent: 16,
            variable_count: 0,
            max_layer_width: 1,
            rho0: WhirRational::default(),
            fixed_extension_degree: 0,
            max_extension_degree: 0,
            max_domain_size: 0,
            max_n0_search_steps: 100_000,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct WhirUniqueDecodingLayer {
    pub layer_index: u64,
    pub variable_count: u64,
    pub width: u64,
    pub domain_size: u64,
    pub rate_numerator: u64,
    pub rate_denominator: u64,
    pub rate: f64,
    pub delta: f64,
    pub repetition_count: u64,
    pub sumcheck_degree_bound: u64,
    pub folding_algebra_bound: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct WhirUniqueDecodingPublicParameters {
    pub base_prime: u64,
    pub ring_exponent: u64,
    pub extension_degree: u64,
    pub initial_domain_size: u64,
    pub variable_count: u64,
    pub layer_widths: Vec<u64>,
    pub shift_repetitions: Vec<u64>,
    pub final_repetitions: u64,
    pub degree_bounds: Vec<u64>,
    pub rates: Vec<f64>,
    pub deltas: Vec<f64>,
    pub lambda_target: u64,
}

impl Default for WhirUniqueDecodingPublicParameters {
    fn default() -> Self {
        Self {
            base_prime: 2,
            ring_exponent: 0,
            extension_degree: 0,
            initial_domain_size: 0,
            variable_count: 0,
            layer_widths: Vec::new(),
            shift_repetitions: Vec::new(),
            final_repetitions: 0,
            degree_bounds: Vec::new(),
            rates: Vec::new(),
            deltas: Vec::new(),
            lambda_target: 0,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct WhirUniqueDecodingSelection {
    pub feasible: bool,
    pub public_params: WhirUniqueDecodingPublicParameters,
    pub required_3_adic_power: u64,
    pub rdom: u64,
    pub rsec: u64,
    pub selected_r: u64,
    pub repetition_security_bits: u64,
    pub effective_security_bits: u64,
    pub algebraic_bound: String,
    pub algebraic_error_log2: f64,
    pub total_error_log2: f64,
    pub layers: Vec<WhirUniqueDecodingLayer>,
    pub notes: Vec<String>,
}

impl Default for WhirUniqueDecodingSelection {
    fn default() -> Self {
        Self {
            feasible: false,
            public_params: WhirUniqueDecodingPublicParameters::default(),
            required_3_adic_power: 0,
            rdom: 0,
            rsec: 0,
            selected_r: 0,
            repetition_security_bits: 0,
            effective_security_bits: 0,
            algebraic_bound: String::new(),
            algebraic_error_log2: 0.0,
            total_error_log2: 0.0,
            layers: Vec::new(),
            notes: Vec::new(),
        }
    }
}

struct CandidateAnalysis {
    selection: WhirUniqueDecodingSelection,
    valid: bool,
}

struct LayerAnalysis {
    summary: WhirUniqueDecodingLayer,
    algebraic_bound_term: u128,
    repetition_log2: f64,
}

pub fn select_whir_unique_decoding_parameters(
    inputs: &WhirUniqueDecodingInputs,
) -> Result<WhirUniqueDecodingSelection> {
    validate_inputs(inputs)?;

    let mut result = WhirUniqueDecodingSelection::default();
    result.notes.push(
        "WHIR selector targets only the p=2 unique-decoding GR PCS mode; optional OOD rounds are redundant consistency constraints and are not the finite-field WHIR list-decoding/OOD-uniqueness argument."
            .to_owned(),
    );
    result.notes.push(
        "Layer domains use n_i = n0 / 3^i, matching the Phase 6 H_i.pow_map(3) oracle chain; width-b_i shift domains require n_i divisible by 3^b_i."
            .to_owned(),
    );
    result.notes.push(
        "The ring exponent s is carried into public parameters; the Section 9 algebraic field-size search depends on the Teichmuller size 2^r."
            .to_owned(),
    );
    result
        .notes
        .push("Unique-decoding thresholds use the half-gap value delta_i=(1-rho_i)/2.".to_owned());
    if inputs.fixed_extension_degree != 0 {
        result.notes.push(format!(
            "WHIR selector uses caller-fixed extension degree r={} and searches only compatible domains.",
            inputs.fixed_extension_degree
        ));
    }

    let layer_widths = layer_widths(inputs.variable_count, inputs.max_layer_width);
    let required_power = required_three_adic_power(&layer_widths)?;
    let required_divisor = pow3_checked(required_power)?;
    let pow3_m = pow3_checked(inputs.variable_count)?;
    let lower_bound = lower_bound_for_n0(inputs, pow3_m)?;
    if lower_bound <= 1 {
        return Err(GrError::InvalidDomain(
            "WHIR n0 lower bound must exceed one",
        ));
    }

    let mut candidate = first_odd_multiple_at_least(lower_bound, required_divisor)?;
    let step = checked_mul(2, required_divisor, "WHIR n0 step")?;
    for _ in 0..inputs.max_n0_search_steps {
        if inputs.max_domain_size != 0 && candidate > inputs.max_domain_size {
            result
                .notes
                .push("no WHIR n0 candidate fits the max_domain_size guard".to_owned());
            return Ok(result);
        }

        let analysis = analyze_candidate(inputs, &layer_widths, candidate)?;
        if analysis.valid {
            if !analysis.selection.feasible {
                for note in analysis.selection.notes {
                    result
                        .notes
                        .push(format!("skipping WHIR n0 candidate {candidate}: {note}"));
                }
                result.notes.push(
                    "skipping a WHIR n0 candidate because its soundness envelope is still infeasible"
                        .to_owned(),
                );
            } else if inputs.max_extension_degree == 0
                || analysis.selection.selected_r <= inputs.max_extension_degree
            {
                let mut selection = analysis.selection;
                selection.notes.splice(0..0, result.notes);
                return Ok(selection);
            } else {
                result.notes.push(
                    "skipping a WHIR n0 candidate because selected r exceeds the max_extension_degree guard"
                        .to_owned(),
                );
            }
        }

        candidate = checked_add(candidate, step, "WHIR n0 search")?;
    }

    result
        .notes
        .push("WHIR n0 search exhausted max_n0_search_steps".to_owned());
    Ok(result)
}

pub fn multiplicative_order_mod_odd(modulus: u64, base: u64) -> Result<u64> {
    if modulus <= 1 || modulus.is_multiple_of(2) {
        return Err(GrError::InvalidDomain(
            "multiplicative_order_mod_odd requires an odd modulus greater than one",
        ));
    }
    if gcd(base, modulus) != 1 {
        return Err(GrError::InvalidDomain(
            "multiplicative_order_mod_odd requires coprime base and modulus",
        ));
    }

    let mut order = euler_phi(modulus)?;
    for factor in unique_prime_factors(order) {
        while order.is_multiple_of(factor) && pow_mod(base, order / factor, modulus) == 1 {
            order /= factor;
        }
    }
    Ok(order)
}

pub fn domain_divides_teichmuller_group(domain_size: u64, extension_degree: u64) -> bool {
    if domain_size == 0 || domain_size.is_multiple_of(2) {
        return false;
    }
    if domain_size == 1 {
        return true;
    }
    pow_mod(2, extension_degree, domain_size) == 1
}

fn analyze_candidate(
    inputs: &WhirUniqueDecodingInputs,
    layer_widths: &[u64],
    n0: u64,
) -> Result<CandidateAnalysis> {
    let mut candidate = CandidateAnalysis {
        selection: WhirUniqueDecodingSelection::default(),
        valid: false,
    };
    let selection = &mut candidate.selection;
    selection.public_params.base_prime = 2;
    selection.public_params.ring_exponent = inputs.ring_exponent;
    selection.public_params.initial_domain_size = n0;
    selection.public_params.variable_count = inputs.variable_count;
    selection.public_params.lambda_target = inputs.lambda_target;
    selection.required_3_adic_power = required_three_adic_power(layer_widths)?;
    selection.repetition_security_bits = checked_add(
        checked_add(inputs.lambda_target, 1, "WHIR repetition target bits")?,
        ceil_log2_u64(checked_add(
            layer_widths.len() as u64,
            1,
            "WHIR repetition layer count",
        )?),
        "WHIR repetition target bits",
    )?;

    let mut algebraic_bound = 0u128;
    let mut remaining_variables = inputs.variable_count;
    let mut log2_error_terms = Vec::new();

    for (layer, &width) in layer_widths.iter().enumerate() {
        let Some(layer_analysis) =
            analyze_layer(selection, layer as u64, width, n0, remaining_variables)?
        else {
            return Ok(candidate);
        };
        algebraic_bound = checked_add_u128(
            algebraic_bound,
            layer_analysis.algebraic_bound_term,
            "WHIR algebraic bound",
        )?;
        log2_error_terms.push(layer_analysis.repetition_log2);

        selection.public_params.layer_widths.push(width);
        selection
            .public_params
            .shift_repetitions
            .push(layer_analysis.summary.repetition_count);
        selection
            .public_params
            .degree_bounds
            .push(WHIR_SUMCHECK_DEGREE_BOUND);
        selection
            .public_params
            .rates
            .push(layer_analysis.summary.rate);
        selection
            .public_params
            .deltas
            .push(layer_analysis.summary.delta);
        selection.layers.push(layer_analysis.summary);

        remaining_variables -= width;
    }

    if !add_final_repetition(
        selection,
        n0,
        layer_widths.len() as u64,
        &mut log2_error_terms,
    )? {
        return Ok(candidate);
    }
    if algebraic_bound == 0 {
        selection
            .notes
            .push("WHIR algebraic bound is unexpectedly zero".to_owned());
        return Ok(candidate);
    }

    finalize_candidate_soundness(
        selection,
        inputs,
        n0,
        algebraic_bound,
        &mut log2_error_terms,
    )?;
    candidate.valid = true;
    Ok(candidate)
}

fn add_final_repetition(
    selection: &mut WhirUniqueDecodingSelection,
    n0: u64,
    layer_count: u64,
    log2_error_terms: &mut Vec<f64>,
) -> Result<bool> {
    let final_domain_divisor = pow3_checked(layer_count)?;
    if !n0.is_multiple_of(final_domain_divisor) {
        selection.notes.push(
            "candidate n0 is not divisible by the final WHIR layer domain divisor".to_owned(),
        );
        return Ok(false);
    }

    let final_domain_size = n0 / final_domain_divisor;
    if final_domain_size <= 1 {
        selection
            .notes
            .push("candidate final WHIR domain is too small for constant checks".to_owned());
        return Ok(false);
    }

    let final_rho = 1.0 / final_domain_size as f64;
    let final_delta = 0.5 * (1.0 - final_rho);
    let final_repetitions =
        repetition_count_for_bits(final_delta, selection.repetition_security_bits)?;
    selection.public_params.final_repetitions = final_repetitions;
    log2_error_terms.push(final_repetitions as f64 * (1.0 - final_delta).log2());
    Ok(true)
}

fn finalize_candidate_soundness(
    selection: &mut WhirUniqueDecodingSelection,
    inputs: &WhirUniqueDecodingInputs,
    n0: u64,
    algebraic_bound: u128,
    log2_error_terms: &mut Vec<f64>,
) -> Result<()> {
    selection.rdom = multiplicative_order_mod_odd(n0, 2)?;
    let algebraic_bits = ceil_log2_u128(algebraic_bound);
    selection.rsec = checked_add(
        checked_add(inputs.lambda_target, 1, "WHIR rsec")?,
        algebraic_bits,
        "WHIR rsec",
    )?;
    selection.selected_r = if inputs.fixed_extension_degree != 0 {
        inputs.fixed_extension_degree
    } else {
        checked_mul(
            selection.rdom,
            ceil_div(selection.rsec, selection.rdom)?,
            "WHIR selected r",
        )?
    };
    selection.public_params.extension_degree = selection.selected_r;
    selection.algebraic_bound = algebraic_bound.to_string();
    selection.algebraic_error_log2 = log2_u128(algebraic_bound)? - selection.selected_r as f64;
    log2_error_terms.push(selection.algebraic_error_log2);
    selection.total_error_log2 = log2_sum(log2_error_terms);
    selection.effective_security_bits = security_bits_from_log2_error(selection.total_error_log2);
    selection.feasible = selection.effective_security_bits >= inputs.lambda_target
        && selection.algebraic_error_log2 < -(inputs.lambda_target as f64)
        && domain_divides_teichmuller_group(n0, selection.selected_r);

    if !selection.feasible {
        add_infeasible_candidate_notes(selection, inputs.lambda_target, n0);
    }
    Ok(())
}

fn add_infeasible_candidate_notes(
    selection: &mut WhirUniqueDecodingSelection,
    lambda_target: u64,
    n0: u64,
) {
    if !domain_divides_teichmuller_group(n0, selection.selected_r) {
        selection
            .notes
            .push("candidate n0 does not divide the fixed Teichmuller group size 2^r-1".to_owned());
    }
    if selection.algebraic_error_log2 >= -(lambda_target as f64) {
        selection
            .notes
            .push("candidate fixed r is too small for the WHIR algebraic error target".to_owned());
    }
    selection
        .notes
        .push("candidate does not meet the WHIR unique-decoding soundness target".to_owned());
}

fn analyze_layer(
    selection: &mut WhirUniqueDecodingSelection,
    layer: u64,
    width: u64,
    n0: u64,
    remaining_variables: u64,
) -> Result<Option<LayerAnalysis>> {
    let layer_divisor = pow3_checked(layer)?;
    if !n0.is_multiple_of(layer_divisor) {
        selection.notes.push(
            "candidate n0 is not divisible by the current WHIR layer domain divisor".to_owned(),
        );
        return Ok(None);
    }

    let domain_size = n0 / layer_divisor;
    let shift_divisor = pow3_checked(width)?;
    if !domain_size.is_multiple_of(shift_divisor) {
        selection.notes.push(
            "candidate n0 does not leave enough 3-adic divisibility for a WHIR shift domain"
                .to_owned(),
        );
        return Ok(None);
    }

    let rate_numerator = pow3_checked(remaining_variables)?;
    let rate_denominator = domain_size;
    if rate_numerator >= rate_denominator {
        selection
            .notes
            .push("candidate leaves the WHIR unique-decoding rate regime rho_i < 1".to_owned());
        return Ok(None);
    }

    let rho = rate_numerator as f64 / rate_denominator as f64;
    let delta = 0.5 * (1.0 - rho);
    if !(delta > 0.0 && delta <= 0.5 * (1.0 - rho)) {
        selection
            .notes
            .push("candidate leaves the WHIR half-gap unique-decoding regime".to_owned());
        return Ok(None);
    }

    let repetitions = repetition_count_for_bits(delta, selection.repetition_security_bits)?;
    let afold = folding_algebra_bound(width, domain_size)?;
    let algebraic_bound_term = checked_add_u128(
        checked_add_u128(
            checked_add_u128(
                afold,
                u128::from(width) * u128::from(WHIR_SUMCHECK_DEGREE_BOUND),
                "WHIR algebraic bound",
            )?,
            u128::from(repetitions),
            "WHIR algebraic bound",
        )?,
        1,
        "WHIR algebraic bound",
    )?;

    Ok(Some(LayerAnalysis {
        summary: WhirUniqueDecodingLayer {
            layer_index: layer,
            variable_count: remaining_variables,
            width,
            domain_size,
            rate_numerator,
            rate_denominator,
            rate: rho,
            delta,
            repetition_count: repetitions,
            sumcheck_degree_bound: WHIR_SUMCHECK_DEGREE_BOUND,
            folding_algebra_bound: afold.to_string(),
        },
        algebraic_bound_term,
        repetition_log2: repetitions as f64 * (1.0 - delta).log2(),
    }))
}

fn folding_algebra_bound(width: u64, domain_size: u64) -> Result<u128> {
    let mut afold = 0u128;
    for j in 0..width {
        let denominator = pow3_checked(j + 1)?;
        let folded_domain = domain_size / denominator;
        let folded = u128::from(folded_domain);
        afold = checked_add_u128(afold, 2 * folded * folded, "WHIR folding bound")?;
    }
    Ok(afold)
}

const fn validate_inputs(inputs: &WhirUniqueDecodingInputs) -> Result<()> {
    if inputs.lambda_target == 0 {
        return Err(GrError::InvalidDomain(
            "WHIR lambda_target must be non-zero",
        ));
    }
    if inputs.ring_exponent == 0 {
        return Err(GrError::InvalidDomain(
            "WHIR ring exponent s must be non-zero",
        ));
    }
    if inputs.variable_count == 0 {
        return Err(GrError::InvalidDomain(
            "WHIR variable count m must be non-zero",
        ));
    }
    if inputs.max_layer_width == 0 {
        return Err(GrError::InvalidDomain("WHIR bmax must be non-zero"));
    }
    if inputs.max_n0_search_steps == 0 {
        return Err(GrError::InvalidDomain(
            "WHIR n0 search step guard must be non-zero",
        ));
    }
    if inputs.fixed_extension_degree != 0
        && inputs.max_extension_degree != 0
        && inputs.fixed_extension_degree > inputs.max_extension_degree
    {
        return Err(GrError::InvalidDomain(
            "WHIR fixed_extension_degree exceeds max_extension_degree guard",
        ));
    }
    validate_open_unit_rational(inputs.rho0, "WHIR rho0")
}

const fn validate_open_unit_rational(value: WhirRational, label: &'static str) -> Result<()> {
    if value.denominator == 0 {
        return Err(GrError::InvalidDomain("WHIR rational denominator is zero"));
    }
    if value.numerator == 0 {
        return Err(GrError::InvalidDomain(label));
    }
    if value.numerator >= value.denominator {
        return Err(GrError::InvalidDomain(label));
    }
    Ok(())
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

fn required_three_adic_power(layer_widths: &[u64]) -> Result<u64> {
    let mut required = 0;
    for (layer, &width) in layer_widths.iter().enumerate() {
        required = required.max(checked_add(layer as u64, width, "required 3-adic power")?);
    }
    Ok(required)
}

fn lower_bound_for_n0(inputs: &WhirUniqueDecodingInputs, pow3_m: u64) -> Result<u64> {
    let numerator = u128::from(pow3_m) * u128::from(inputs.rho0.denominator);
    ceil_div_u128_to_u64(numerator, inputs.rho0.numerator)
}

fn first_odd_multiple_at_least(lower_bound: u64, divisor: u64) -> Result<u64> {
    let mut quotient = ceil_div(lower_bound, divisor)?;
    if quotient.is_multiple_of(2) {
        quotient = checked_add(quotient, 1, "odd multiple quotient")?;
    }
    checked_mul(quotient, divisor, "n0 candidate")
}

#[expect(
    clippy::cast_sign_loss,
    clippy::while_float,
    reason = "C++ WHIR selector computes repetitions from checked finite positive floating estimates"
)]
fn repetition_count_for_bits(delta: f64, target_bits: u64) -> Result<u64> {
    if !(delta > 0.0 && delta < 1.0) {
        return Err(GrError::InvalidDomain("WHIR delta must lie in (0, 1)"));
    }
    let denominator = -(1.0 - delta).log2();
    if !denominator.is_finite() || denominator <= 0.0 {
        return Err(GrError::InvalidDomain(
            "WHIR repetition denominator is invalid",
        ));
    }

    let estimate = ((target_bits as f64) / denominator).ceil();
    if !estimate.is_finite() || estimate < 1.0 || estimate > u64::MAX as f64 {
        return Err(GrError::ArithmeticOverflow("WHIR repetition count"));
    }

    let mut repetitions = estimate as u64;
    while repetitions as f64 * denominator < target_bits as f64 {
        repetitions = checked_add(repetitions, 1, "WHIR repetition count")?;
    }
    while repetitions > 1 && (repetitions - 1) as f64 * denominator >= target_bits as f64 {
        repetitions -= 1;
    }
    Ok(repetitions)
}

fn log2_sum(log2_terms: &[f64]) -> f64 {
    let max_log = log2_terms
        .iter()
        .copied()
        .filter(|term| term.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    if !max_log.is_finite() {
        return max_log;
    }

    let scaled_sum = log2_terms
        .iter()
        .copied()
        .filter(|term| term.is_finite())
        .map(|term| 2.0f64.powf(term - max_log))
        .sum::<f64>();
    max_log + scaled_sum.log2()
}

#[expect(
    clippy::cast_sign_loss,
    reason = "security bits are floored only after checked non-negative finite bounds"
)]
fn security_bits_from_log2_error(log2_error: f64) -> u64 {
    if !log2_error.is_finite() || log2_error >= 0.0 {
        return 0;
    }
    let bits = -log2_error;
    if bits >= u64::MAX as f64 {
        return u64::MAX;
    }
    bits.floor() as u64
}

fn log2_u128(value: u128) -> Result<f64> {
    if value == 0 {
        return Err(GrError::InvalidDomain(
            "log2 is undefined for non-positive integers",
        ));
    }
    let bit_count = u128::BITS - value.leading_zeros();
    let kept_bits = bit_count.min(64);
    let top = value >> (bit_count - kept_bits);
    let exponent =
        i32::try_from(kept_bits - 1).map_err(|_| GrError::ArithmeticOverflow("log2 kept bits"))?;
    let scaled = top as f64 / 2.0f64.powi(exponent);
    Ok(f64::from(bit_count - 1) + scaled.log2())
}

const fn ceil_log2_u64(mut value: u64) -> u64 {
    if value <= 1 {
        return 0;
    }
    let mut bits = 0;
    value -= 1;
    while value != 0 {
        bits += 1;
        value >>= 1;
    }
    bits
}

fn ceil_log2_u128(value: u128) -> u64 {
    if value <= 1 {
        0
    } else {
        u64::from(u128::BITS - (value - 1).leading_zeros())
    }
}

fn unique_prime_factors(mut value: u64) -> Vec<u64> {
    let mut factors = Vec::new();
    if value.is_multiple_of(2) {
        factors.push(2);
        while value.is_multiple_of(2) {
            value >>= 1;
        }
    }
    let mut divisor = 3;
    while divisor <= value / divisor {
        if value.is_multiple_of(divisor) {
            factors.push(divisor);
            while value.is_multiple_of(divisor) {
                value /= divisor;
            }
        }
        divisor += 2;
    }
    if value > 1 {
        factors.push(value);
    }
    factors
}

fn euler_phi(value: u64) -> Result<u64> {
    let mut phi = value;
    for factor in unique_prime_factors(value) {
        phi = checked_mul(phi / factor, factor - 1, "Euler phi")?;
    }
    Ok(phi)
}

fn pow_mod(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    if modulus == 1 {
        return 0;
    }
    let mut result = 1 % modulus;
    base %= modulus;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = mul_mod(result, base, modulus);
        }
        exponent >>= 1;
        if exponent != 0 {
            base = mul_mod(base, base, modulus);
        }
    }
    result
}

fn mul_mod(lhs: u64, rhs: u64, modulus: u64) -> u64 {
    ((u128::from(lhs) * u128::from(rhs)) % u128::from(modulus)) as u64
}

const fn gcd(mut lhs: u64, mut rhs: u64) -> u64 {
    while rhs != 0 {
        let remainder = lhs % rhs;
        lhs = rhs;
        rhs = remainder;
    }
    lhs
}

fn pow3_checked(exponent: u64) -> Result<u64> {
    let mut out = 1u64;
    for _ in 0..exponent {
        out = checked_mul(out, 3, "power of 3")?;
    }
    Ok(out)
}

fn ceil_div(numerator: u64, denominator: u64) -> Result<u64> {
    if denominator == 0 {
        return Err(GrError::InvalidDomain("division by zero"));
    }
    Ok(numerator / denominator + u64::from(!numerator.is_multiple_of(denominator)))
}

fn ceil_div_u128_to_u64(numerator: u128, denominator: u64) -> Result<u64> {
    if denominator == 0 {
        return Err(GrError::InvalidDomain("division by zero"));
    }
    let denominator = u128::from(denominator);
    let quotient = numerator / denominator + u128::from(!numerator.is_multiple_of(denominator));
    u64::try_from(quotient).map_err(|_| GrError::ArithmeticOverflow("ceil division result"))
}

fn checked_add(lhs: u64, rhs: u64, label: &'static str) -> Result<u64> {
    lhs.checked_add(rhs)
        .ok_or(GrError::ArithmeticOverflow(label))
}

fn checked_mul(lhs: u64, rhs: u64, label: &'static str) -> Result<u64> {
    lhs.checked_mul(rhs)
        .ok_or(GrError::ArithmeticOverflow(label))
}

fn checked_add_u128(lhs: u128, rhs: u128, label: &'static str) -> Result<u128> {
    lhs.checked_add(rhs)
        .ok_or(GrError::ArithmeticOverflow(label))
}

#[cfg(test)]
mod tests {
    use crate::{
        algebra::galois_ring::GrError,
        protocols::whir_gr::soundness::{
            domain_divides_teichmuller_group, multiplicative_order_mod_odd,
            select_whir_unique_decoding_parameters, WhirRational, WhirUniqueDecodingInputs,
        },
    };

    fn small_inputs() -> WhirUniqueDecodingInputs {
        WhirUniqueDecodingInputs {
            lambda_target: 32,
            ring_exponent: 16,
            variable_count: 3,
            max_layer_width: 1,
            rho0: WhirRational {
                numerator: 1,
                denominator: 3,
            },
            ..WhirUniqueDecodingInputs::default()
        }
    }

    #[test]
    fn invalid_rates_should_reject() {
        let mut inputs = small_inputs();
        inputs.rho0 = WhirRational {
            numerator: 0,
            denominator: 1,
        };
        assert!(matches!(
            select_whir_unique_decoding_parameters(&inputs),
            Err(GrError::InvalidDomain(_))
        ));

        inputs = small_inputs();
        inputs.rho0 = WhirRational {
            numerator: 1,
            denominator: 1,
        };
        assert!(matches!(
            select_whir_unique_decoding_parameters(&inputs),
            Err(GrError::InvalidDomain(_))
        ));
    }

    #[test]
    fn small_smoke_params_should_match_cpp_selector() {
        let selected = select_whir_unique_decoding_parameters(&small_inputs()).unwrap();

        assert!(selected.feasible);
        assert_eq!(selected.public_params.base_prime, 2);
        assert_eq!(selected.public_params.ring_exponent, 16);
        assert_eq!(selected.public_params.initial_domain_size, 81);
        assert_eq!(selected.public_params.layer_widths.len(), 3);
        assert_eq!(selected.public_params.shift_repetitions.len(), 3);
        assert!(selected.public_params.final_repetitions > 0);
        assert_eq!(selected.repetition_security_bits, 35);
        assert!(selected.effective_security_bits >= 32);
    }

    #[test]
    fn selected_n0_should_divide_selected_teichmuller_group() {
        let selected = select_whir_unique_decoding_parameters(&small_inputs()).unwrap();

        assert_ne!(selected.rdom, 0);
        assert_eq!(selected.selected_r % selected.rdom, 0);
        assert!(domain_divides_teichmuller_group(
            selected.public_params.initial_domain_size,
            selected.selected_r
        ));
        assert_eq!(
            multiplicative_order_mod_odd(selected.public_params.initial_domain_size, 2).unwrap(),
            selected.rdom
        );
    }

    #[test]
    fn fixed_extension_degree_should_match_cpp_selector() {
        let mut inputs = small_inputs();
        inputs.fixed_extension_degree = 54;
        let selected = select_whir_unique_decoding_parameters(&inputs).unwrap();

        assert!(selected.feasible);
        assert_eq!(selected.selected_r, 54);
        assert_eq!(selected.public_params.extension_degree, 54);
        assert_eq!(selected.public_params.initial_domain_size, 81);
        assert!(domain_divides_teichmuller_group(
            selected.public_params.initial_domain_size,
            selected.selected_r
        ));
    }

    #[test]
    fn incompatible_fixed_extension_degree_should_report_infeasible() {
        let mut inputs = small_inputs();
        inputs.fixed_extension_degree = 55;
        inputs.max_n0_search_steps = 3;
        let selected = select_whir_unique_decoding_parameters(&inputs).unwrap();

        assert!(!selected.feasible);
    }

    #[test]
    fn conflicting_fixed_extension_guard_should_reject() {
        let mut inputs = small_inputs();
        inputs.fixed_extension_degree = 54;
        inputs.max_extension_degree = 53;

        assert!(matches!(
            select_whir_unique_decoding_parameters(&inputs),
            Err(GrError::InvalidDomain(_))
        ));
    }

    #[test]
    fn required_domains_should_divide_through_round_chain() {
        let mut inputs = small_inputs();
        inputs.variable_count = 4;
        inputs.max_layer_width = 2;
        let selected = select_whir_unique_decoding_parameters(&inputs).unwrap();

        assert!(selected.feasible);
        assert_eq!(selected.required_3_adic_power, 3);
        assert_eq!(selected.public_params.initial_domain_size, 243);
        assert_eq!(selected.layers.len(), 2);
        assert_eq!(selected.layers[0].width, 2);
        assert_eq!(selected.layers[1].width, 2);

        for layer in selected.layers {
            assert_eq!(layer.domain_size % 3u64.pow(layer.width as u32), 0);
            assert!(layer.rate > 0.0);
            assert!(layer.rate < 1.0);
            assert!(layer.delta > 0.0);
            assert!((layer.delta - 0.5 * (1.0 - layer.rate)).abs() < 1.0e-18);
        }
    }

    #[test]
    fn algebraic_bound_should_be_below_target() {
        let selected = select_whir_unique_decoding_parameters(&small_inputs()).unwrap();

        assert!(selected.feasible);
        assert!(selected.algebraic_error_log2 < -(small_inputs().lambda_target as f64));
        assert!(selected.total_error_log2 <= -(small_inputs().lambda_target as f64));
        assert!(selected.effective_security_bits >= small_inputs().lambda_target);
    }

    #[test]
    fn benchmark_guards_should_report_infeasible() {
        let mut inputs = small_inputs();
        inputs.max_domain_size = 27;
        let selected = select_whir_unique_decoding_parameters(&inputs).unwrap();

        assert!(!selected.feasible);
        assert!(selected
            .notes
            .iter()
            .any(|note| note.contains("max_domain_size")));

        inputs = small_inputs();
        inputs.max_extension_degree = 1;
        inputs.max_n0_search_steps = 2;
        let selected = select_whir_unique_decoding_parameters(&inputs).unwrap();

        assert!(!selected.feasible);
        assert!(selected
            .notes
            .iter()
            .any(|note| note.contains("max_extension_degree")));
    }

    #[test]
    fn whir_preset_rows_should_match_cpp_selector() {
        let expected = [
            (4, 189),
            (5, 513),
            (6, 1_539),
            (7, 4_617),
            (8, 13_203),
            (9, 39_609),
            (10, 124_173),
        ];

        for (variable_count, expected_n0) in expected {
            let inputs = WhirUniqueDecodingInputs {
                lambda_target: 128,
                ring_exponent: 16,
                variable_count,
                max_layer_width: 3,
                rho0: WhirRational {
                    numerator: 1,
                    denominator: 2,
                },
                fixed_extension_degree: 162,
                ..WhirUniqueDecodingInputs::default()
            };
            let selected = select_whir_unique_decoding_parameters(&inputs).unwrap();

            assert!(selected.feasible);
            assert_eq!(selected.selected_r, 162);
            assert_eq!(selected.public_params.initial_domain_size, expected_n0);
        }
    }
}
