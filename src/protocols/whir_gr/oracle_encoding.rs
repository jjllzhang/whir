#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::algebra::galois_ring::{clear_elem, Domain, GrContext, GrElem, GrError, Result};

const MAX_LOCAL_RADIX: u64 = 257;
const PARALLEL_DFT_THRESHOLD: usize = 1024;

pub(in crate::protocols::whir_gr) fn rs_encode_teichmuller_coset(
    ctx: &GrContext,
    domain: &Domain,
    coefficients: &[GrElem],
) -> Result<Option<Vec<GrElem>>> {
    if domain.context().config() != ctx.config() {
        return Err(GrError::DifferentRings);
    }
    if !domain.is_teichmuller_subset() {
        return Ok(None);
    }

    let size = checked_usize(domain.size(), "domain size")?;
    if coefficients.len() > size {
        return Ok(None);
    }

    let Some(factors) = supported_radix_factorization(domain.size()) else {
        return Ok(None);
    };

    let scaled = scale_coefficients_for_offset(ctx, domain, coefficients, size);
    Ok(Some(dft_mixed_radix(
        ctx,
        &scaled,
        domain.root(),
        &factors,
    )?))
}

fn scale_coefficients_for_offset(
    ctx: &GrContext,
    domain: &Domain,
    coefficients: &[GrElem],
    size: usize,
) -> Vec<GrElem> {
    let mut scaled = vec![ctx.zero(); size];
    if coefficients.is_empty() {
        return scaled;
    }

    if *domain.offset() == ctx.one() {
        scaled[..coefficients.len()].clone_from_slice(coefficients);
        return scaled;
    }

    let mut offset_power = ctx.one();
    let mut next_power = ctx.zero();
    let mut scratch = vec![0; ctx.mul_scratch_len()];
    for (index, coefficient) in coefficients.iter().enumerate() {
        ctx.mul_into(&mut scaled[index], coefficient, &offset_power, &mut scratch);
        if index + 1 < coefficients.len() {
            ctx.mul_into(
                &mut next_power,
                &offset_power,
                domain.offset(),
                &mut scratch,
            );
            std::mem::swap(&mut offset_power, &mut next_power);
        }
    }
    scaled
}

fn dft_mixed_radix(
    ctx: &GrContext,
    values: &[GrElem],
    root: &GrElem,
    factors: &[usize],
) -> Result<Vec<GrElem>> {
    if values.len() == 1 {
        return Ok(values.to_vec());
    }

    let Some((&radix, remaining_factors)) = factors.split_first() else {
        return Err(GrError::InvalidDomain("missing DFT factorization"));
    };
    if radix == 0 || !values.len().is_multiple_of(radix) {
        return Err(GrError::InvalidDomain("invalid DFT factorization"));
    }

    #[cfg(feature = "parallel")]
    {
        if values.len() >= PARALLEL_DFT_THRESHOLD && rayon::current_num_threads() > 1 {
            return dft_mixed_radix_parallel(ctx, values, root, radix, remaining_factors);
        }
    }

    dft_mixed_radix_sequential(ctx, values, root, radix, remaining_factors)
}

#[cfg(feature = "parallel")]
fn dft_mixed_radix_parallel(
    ctx: &GrContext,
    values: &[GrElem],
    root: &GrElem,
    radix: usize,
    remaining_factors: &[usize],
) -> Result<Vec<GrElem>> {
    let branches = (0..radix)
        .into_par_iter()
        .map(|output_residue| {
            dft_output_residue(ctx, values, root, radix, remaining_factors, output_residue)
                .map(|branch| (output_residue, branch))
        })
        .collect::<Vec<_>>();

    let mut out = vec![ctx.zero(); values.len()];
    for branch in branches {
        let (output_residue, values_for_residue) = branch?;
        scatter_output_residue(&mut out, radix, output_residue, values_for_residue);
    }
    Ok(out)
}

fn dft_mixed_radix_sequential(
    ctx: &GrContext,
    values: &[GrElem],
    root: &GrElem,
    radix: usize,
    remaining_factors: &[usize],
) -> Result<Vec<GrElem>> {
    let mut out = vec![ctx.zero(); values.len()];
    for output_residue in 0..radix {
        let values_for_residue =
            dft_output_residue(ctx, values, root, radix, remaining_factors, output_residue)?;
        scatter_output_residue(&mut out, radix, output_residue, values_for_residue);
    }
    Ok(out)
}

fn dft_output_residue(
    ctx: &GrContext,
    values: &[GrElem],
    root: &GrElem,
    radix: usize,
    remaining_factors: &[usize],
    output_residue: usize,
) -> Result<Vec<GrElem>> {
    let inner_size = values.len() / radix;
    let root_for_radix = ctx.pow(root, inner_size as u128);
    let root_for_inner = ctx.pow(root, radix as u128);
    let small_step = ctx.pow(&root_for_radix, output_residue as u128);
    let twiddle_step = ctx.pow(root, output_residue as u128);

    let mut sequence = vec![ctx.zero(); inner_size];
    let mut twiddle = ctx.one();
    let mut next_twiddle = ctx.zero();
    let mut sum = ctx.zero();
    let mut next_sum = ctx.zero();
    let mut product = ctx.zero();
    let mut small_power = ctx.one();
    let mut next_small_power = ctx.zero();
    let mut mul_scratch = vec![0; ctx.mul_scratch_len()];

    for (inner_index, output) in sequence.iter_mut().enumerate() {
        clear_elem(&mut sum);
        if output_residue != 0 {
            small_power.clone_from(&ctx.one());
        }
        for radix_index in 0..radix {
            let value = &values[inner_index + inner_size * radix_index];
            if output_residue == 0 || radix_index == 0 {
                ctx.add_into(&mut next_sum, &sum, value);
            } else {
                ctx.mul_into(&mut product, value, &small_power, &mut mul_scratch);
                ctx.add_into(&mut next_sum, &sum, &product);
            }
            std::mem::swap(&mut sum, &mut next_sum);

            if output_residue != 0 && radix_index + 1 < radix {
                ctx.mul_into(
                    &mut next_small_power,
                    &small_power,
                    &small_step,
                    &mut mul_scratch,
                );
                std::mem::swap(&mut small_power, &mut next_small_power);
            }
        }

        if output_residue == 0 {
            output.clone_from(&sum);
        } else {
            ctx.mul_into(output, &sum, &twiddle, &mut mul_scratch);
            if inner_index + 1 < inner_size {
                ctx.mul_into(&mut next_twiddle, &twiddle, &twiddle_step, &mut mul_scratch);
                std::mem::swap(&mut twiddle, &mut next_twiddle);
            }
        }
    }

    dft_mixed_radix(ctx, &sequence, &root_for_inner, remaining_factors)
}

fn scatter_output_residue(
    out: &mut [GrElem],
    radix: usize,
    output_residue: usize,
    values_for_residue: Vec<GrElem>,
) {
    for (inner_index, value) in values_for_residue.into_iter().enumerate() {
        out[output_residue + radix * inner_index] = value;
    }
}

fn supported_radix_factorization(size: u64) -> Option<Vec<usize>> {
    if size == 0 {
        return None;
    }

    let mut remaining = size;
    let mut factors = Vec::new();
    while remaining.is_multiple_of(3) {
        factors.push(3);
        remaining /= 3;
    }

    let mut divisor = 5;
    while remaining != 1 && divisor <= MAX_LOCAL_RADIX {
        while remaining.is_multiple_of(divisor) {
            factors.push(divisor as usize);
            remaining /= divisor;
        }
        divisor += 2;
    }

    if remaining == 1 {
        Some(factors)
    } else {
        None
    }
}

fn checked_usize(value: u64, label: &'static str) -> Result<usize> {
    usize::try_from(value).map_err(|_| GrError::ArithmeticOverflow(label))
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::{rs_encode_teichmuller_coset, supported_radix_factorization};
    use crate::{
        algebra::galois_ring::{teichmuller_generator, Domain, GrConfig, GrContext},
        protocols::whir_gr::{
            bench_support,
            multiquadratic::{pow3_checked, MultiQuadraticPolynomial},
        },
    };

    fn context() -> Arc<GrContext> {
        Arc::new(
            GrContext::new(GrConfig {
                p: 2,
                k_exp: 16,
                r: 18,
            })
            .unwrap(),
        )
    }

    fn polynomial(ctx: &GrContext, variable_count: u64) -> MultiQuadraticPolynomial {
        let coefficient_count = pow3_checked(variable_count).unwrap();
        let coefficients = (0..coefficient_count)
            .map(|index| ctx.from_u64((7 + 11 * index) % 31))
            .collect();
        MultiQuadraticPolynomial::new(variable_count, coefficients).unwrap()
    }

    fn direct_encode(
        ctx: &GrContext,
        domain: &Domain,
        polynomial: &MultiQuadraticPolynomial,
    ) -> Vec<crate::algebra::galois_ring::GrElem> {
        domain
            .iter_elements()
            .map(|point| polynomial.evaluate_pow(ctx, &point))
            .collect()
    }

    #[test]
    fn supported_radix_factorization_should_cover_current_bench_domains() {
        for case in bench_support::WHIR_GR_CASES {
            let factors = supported_radix_factorization(case.n).unwrap();
            assert_eq!(
                factors.iter().product::<usize>() as u64,
                case.n,
                "{}",
                case.short_name()
            );
        }
    }

    #[test]
    fn rs_encode_should_match_direct_evaluation_on_subgroup() {
        let ctx = context();
        let domain = Domain::teichmuller_subgroup(Arc::clone(&ctx), 27).unwrap();
        let polynomial = polynomial(&ctx, 3);

        let encoded = rs_encode_teichmuller_coset(&ctx, &domain, polynomial.coefficients())
            .unwrap()
            .unwrap();

        assert_eq!(encoded, direct_encode(&ctx, &domain, &polynomial));
    }

    #[test]
    fn rs_encode_should_match_direct_evaluation_on_nontrivial_coset() {
        let ctx = context();
        let offset = teichmuller_generator(&ctx).unwrap();
        let domain = Domain::teichmuller_coset(Arc::clone(&ctx), offset, 27).unwrap();
        let polynomial = polynomial(&ctx, 2);

        let encoded = rs_encode_teichmuller_coset(&ctx, &domain, polynomial.coefficients())
            .unwrap()
            .unwrap();

        assert_eq!(encoded, direct_encode(&ctx, &domain, &polynomial));
    }
}
