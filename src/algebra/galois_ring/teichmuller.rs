use super::{GrContext, GrElem, GrError, Result};

pub fn teichmuller_generator(ctx: &GrContext) -> Result<GrElem> {
    if ctx.config().r == 1 {
        return Ok(ctx.one());
    }

    for attempt in 0..2048 {
        let candidate = teichmuller_projection(ctx, &deterministic_candidate(ctx, attempt)?)?;
        if candidate != ctx.zero() && is_teichmuller_element(ctx, &candidate) {
            return Ok(candidate);
        }
    }

    Err(GrError::InvalidDomain(
        "failed to find deterministic Teichmuller generator candidate",
    ))
}

pub fn teichmuller_group_order_words(ctx: &GrContext) -> Vec<u64> {
    mersenne_words(ctx.config().r)
}

pub fn teichmuller_subgroup_size_supported(ctx: &GrContext, size: u64) -> bool {
    if size == 0 {
        return false;
    }
    if size == 1 {
        return true;
    }
    pow2_mod(ctx.config().r, size) == 1 % size
}

pub fn is_teichmuller_element(ctx: &GrContext, element: &GrElem) -> bool {
    if *element == ctx.zero() {
        return true;
    }
    if !ctx.is_unit(element) {
        return false;
    }
    ctx.pow_words_le(element, &teichmuller_group_order_words(ctx)) == ctx.one()
}

pub fn teichmuller_element_by_index(ctx: &GrContext, index: u128) -> Result<GrElem> {
    if let Some(size) = teichmuller_set_size_u128(ctx) {
        if index >= size {
            return Err(GrError::IndexOutOfRange {
                index: u64::MAX,
                size: u64::MAX,
            });
        }
    }

    if index == 0 {
        return Ok(ctx.zero());
    }

    Ok(ctx.pow(&teichmuller_generator(ctx)?, index - 1))
}

pub fn teichmuller_subgroup_generator(ctx: &GrContext, size: u64) -> Result<GrElem> {
    if !teichmuller_subgroup_size_supported(ctx, size) {
        return Err(GrError::InvalidSubgroupSize { size });
    }
    if size == 1 {
        return Ok(ctx.one());
    }

    let exponent_words = mersenne_div_u64_words(ctx.config().r, size)?;
    for attempt in 0..4096 {
        let base = deterministic_candidate(ctx, attempt)?;
        if base == ctx.zero() {
            continue;
        }
        let projected = teichmuller_projection(ctx, &base)?;
        if projected == ctx.zero() {
            continue;
        }
        let candidate = ctx.pow_words_le(&projected, &exponent_words);
        if has_exact_multiplicative_order(ctx, &candidate, size) {
            return Ok(candidate);
        }
    }

    Err(GrError::InvalidDomain(
        "failed to find Teichmuller subgroup generator",
    ))
}

pub fn has_exact_multiplicative_order(ctx: &GrContext, element: &GrElem, order: u64) -> bool {
    if order == 0 || !ctx.is_unit(element) {
        return false;
    }
    if ctx.pow(element, u128::from(order)) != ctx.one() {
        return false;
    }
    for divisor in prime_divisors(order) {
        if ctx.pow(element, u128::from(order / divisor)) == ctx.one() {
            return false;
        }
    }
    true
}

pub fn generate_teichmuller_subgroup(ctx: &GrContext, size: u64) -> Result<Vec<GrElem>> {
    if size == 0 {
        return Err(GrError::InvalidSubgroupSize { size });
    }

    let root = teichmuller_subgroup_generator(ctx, size)?;
    let mut values = Vec::with_capacity(size as usize);
    let mut current = ctx.one();
    for _ in 0..size {
        values.push(current.clone());
        current = ctx.mul(&current, &root);
    }
    Ok(values)
}

fn teichmuller_projection(ctx: &GrContext, base: &GrElem) -> Result<GrElem> {
    if ctx.config().k_exp == 0 {
        return Err(GrError::ZeroPrecision);
    }
    Ok(ctx.pow(base, 1u128 << (ctx.config().k_exp - 1)))
}

fn deterministic_candidate(ctx: &GrContext, attempt: u64) -> Result<GrElem> {
    match attempt {
        0 => Ok(ctx.x()),
        1 => Ok(ctx.add(&ctx.x(), &ctx.one())),
        2 => Ok(ctx.add(&ctx.x(), &ctx.from_u64(2))),
        _ => {
            let mut bytes = vec![0; ctx.elem_bytes()];
            let mut state = 0x5445_4943_484d_5551u64
                ^ u64::from(ctx.config().k_exp).rotate_left(17)
                ^ (ctx.config().r as u64).rotate_left(31)
                ^ attempt;
            for chunk in bytes.chunks_mut(8) {
                let word = splitmix64(&mut state).to_le_bytes();
                chunk.copy_from_slice(&word[..chunk.len()]);
            }
            ctx.deserialize(&bytes)
        }
    }
}

const fn teichmuller_set_size_u128(ctx: &GrContext) -> Option<u128> {
    if ctx.config().r >= u128::BITS as usize {
        None
    } else {
        Some(1u128 << ctx.config().r)
    }
}

fn mersenne_words(bit_count: usize) -> Vec<u64> {
    if bit_count == 0 {
        return Vec::new();
    }
    let mut words = vec![u64::MAX; bit_count.div_ceil(u64::BITS as usize)];
    let used_bits = bit_count % u64::BITS as usize;
    if used_bits != 0 {
        let last = words.len() - 1;
        words[last] = (1u64 << used_bits) - 1;
    }
    words
}

fn mersenne_div_u64_words(bit_count: usize, divisor: u64) -> Result<Vec<u64>> {
    if divisor == 0 {
        return Err(GrError::InvalidSubgroupSize { size: divisor });
    }

    let mut quotient = vec![0; bit_count.div_ceil(u64::BITS as usize)];
    let mut remainder = 0u128;
    let divisor = u128::from(divisor);
    for bit_index in (0..bit_count).rev() {
        remainder = (remainder << 1) | 1;
        if remainder >= divisor {
            let word_index = bit_index / u64::BITS as usize;
            let inner_bit = bit_index % u64::BITS as usize;
            quotient[word_index] |= 1u64 << inner_bit;
            remainder -= divisor;
        }
    }

    if remainder != 0 {
        return Err(GrError::InvalidSubgroupSize {
            size: divisor as u64,
        });
    }
    Ok(quotient)
}

fn pow2_mod(exponent: usize, modulus: u64) -> u64 {
    if modulus == 1 {
        return 0;
    }

    let mut result = 1 % modulus;
    let mut base = 2 % modulus;
    let mut remaining = exponent;
    while remaining != 0 {
        if remaining & 1 == 1 {
            result = ((u128::from(result) * u128::from(base)) % u128::from(modulus)) as u64;
        }
        remaining >>= 1;
        if remaining != 0 {
            base = ((u128::from(base) * u128::from(base)) % u128::from(modulus)) as u64;
        }
    }
    result
}

const fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn prime_divisors(mut value: u64) -> Vec<u64> {
    let mut divisors = Vec::new();
    let mut factor = 2;
    while factor * factor <= value {
        if value.is_multiple_of(factor) {
            divisors.push(factor);
            while value.is_multiple_of(factor) {
                value /= factor;
            }
        }
        factor += 1;
    }
    if value > 1 {
        divisors.push(value);
    }
    divisors
}

#[cfg(test)]
mod tests {
    use super::{
        generate_teichmuller_subgroup, has_exact_multiplicative_order, is_teichmuller_element,
        teichmuller_element_by_index, teichmuller_subgroup_generator,
        teichmuller_subgroup_size_supported,
    };
    use crate::algebra::galois_ring::{GrConfig, GrContext, GrError};

    fn ctx_r6() -> GrContext {
        GrContext::new(GrConfig {
            p: 2,
            k_exp: 16,
            r: 6,
        })
        .unwrap()
    }

    #[test]
    fn subgroup_size_support_should_match_mersenne_divisibility() {
        let ctx = ctx_r6();

        assert!(teichmuller_subgroup_size_supported(&ctx, 9));
        assert!(teichmuller_subgroup_size_supported(&ctx, 63));
        assert!(!teichmuller_subgroup_size_supported(&ctx, 10));
    }

    #[test]
    fn subgroup_generator_should_have_exact_order() {
        let ctx = ctx_r6();

        let root = teichmuller_subgroup_generator(&ctx, 9).unwrap();

        assert!(has_exact_multiplicative_order(&ctx, &root, 9));
        assert_eq!(ctx.pow(&root, 9), ctx.one());
    }

    #[test]
    fn generated_subgroup_should_enumerate_powers() {
        let ctx = ctx_r6();

        let subgroup = generate_teichmuller_subgroup(&ctx, 9).unwrap();

        assert_eq!(subgroup.len(), 9);
        assert_eq!(subgroup[0], ctx.one());
        assert_eq!(subgroup[8], ctx.pow(&subgroup[1], 8));
    }

    #[test]
    fn teichmuller_membership_should_reject_nilpotent_drift() {
        let ctx = ctx_r6();
        let root = teichmuller_subgroup_generator(&ctx, 9).unwrap();
        let two = ctx.from_u64(2);

        assert!(is_teichmuller_element(&ctx, &ctx.zero()));
        assert!(is_teichmuller_element(&ctx, &ctx.one()));
        assert!(is_teichmuller_element(&ctx, &root));
        assert!(!is_teichmuller_element(&ctx, &two));
    }

    #[test]
    fn teichmuller_element_by_index_should_start_with_zero_and_one() {
        let ctx = ctx_r6();

        assert_eq!(teichmuller_element_by_index(&ctx, 0).unwrap(), ctx.zero());
        assert_eq!(teichmuller_element_by_index(&ctx, 1).unwrap(), ctx.one());
    }

    #[test]
    fn invalid_subgroup_size_should_reject() {
        let ctx = ctx_r6();

        let result = teichmuller_subgroup_generator(&ctx, 10);

        assert!(matches!(result, Err(GrError::InvalidSubgroupSize { .. })));
    }
}
