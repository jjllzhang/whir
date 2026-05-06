#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::algebra::galois_ring::{
    clear_elem, Domain, GrContext, GrElem, GrError, GrScratch, Result,
};

type FoldQueryBatches = Vec<Vec<GrElem>>;

const TERNARY_FOLD_SCRATCH_ELEMS: usize = 12;
const PARALLEL_FOLD_GROUP_THRESHOLD: usize = 64;

pub fn fold_eval(
    ctx: &GrContext,
    fiber_points: &[GrElem],
    fiber_values: &[GrElem],
    alpha: &GrElem,
) -> Result<GrElem> {
    if fiber_points.is_empty() {
        return Err(GrError::InvalidDomain("fold_eval requires non-empty fiber"));
    }
    if fiber_points.len() != fiber_values.len() {
        return Err(GrError::InvalidDomain(
            "fold_eval requires equal-sized fiber inputs",
        ));
    }
    if fiber_points.len() == 3 {
        return fold_eval_ternary(ctx, fiber_points, fiber_values, alpha);
    }

    let mut denominator_products = Vec::with_capacity(fiber_points.len());
    for (index, point) in fiber_points.iter().enumerate() {
        let mut denominator = ctx.one();
        for (other_index, other) in fiber_points.iter().enumerate() {
            if index == other_index {
                continue;
            }
            let difference = ctx.sub(point, other);
            if !ctx.is_unit(&difference) {
                return Err(GrError::InvalidDomain(
                    "fold_eval requires fiber points with unit differences",
                ));
            }
            denominator = ctx.mul(&denominator, &difference);
        }
        denominator_products.push(denominator);
    }
    let denominator_inverses = ctx.batch_inv(&denominator_products)?;

    let mut differences = Vec::with_capacity(fiber_points.len());
    for point in fiber_points {
        differences.push(ctx.sub(alpha, point));
    }

    let mut prefix_products = Vec::with_capacity(fiber_points.len() + 1);
    let mut prefix_product = ctx.one();
    prefix_products.push(prefix_product.clone());
    for difference in &differences {
        prefix_product = ctx.mul(&prefix_product, difference);
        prefix_products.push(prefix_product.clone());
    }

    let mut suffix_products = vec![ctx.one(); fiber_points.len() + 1];
    for index in (0..fiber_points.len()).rev() {
        suffix_products[index] = ctx.mul(&suffix_products[index + 1], &differences[index]);
    }

    let mut out = ctx.zero();
    for index in 0..fiber_points.len() {
        let numerator = ctx.mul(&prefix_products[index], &suffix_products[index + 1]);
        let basis_value = ctx.mul(&numerator, &denominator_inverses[index]);
        out = ctx.add(&out, &ctx.mul(&fiber_values[index], &basis_value));
    }
    Ok(out)
}

fn fold_eval_ternary(
    ctx: &GrContext,
    fiber_points: &[GrElem],
    fiber_values: &[GrElem],
    alpha: &GrElem,
) -> Result<GrElem> {
    let p0 = &fiber_points[0];
    let p1 = &fiber_points[1];
    let p2 = &fiber_points[2];
    let d01 = ctx.sub(p0, p1);
    let d02 = ctx.sub(p0, p2);
    let d10 = ctx.sub(p1, p0);
    let d12 = ctx.sub(p1, p2);
    let d20 = ctx.sub(p2, p0);
    let d21 = ctx.sub(p2, p1);
    if [&d01, &d02, &d10, &d12, &d20, &d21]
        .iter()
        .any(|difference| !ctx.is_unit(difference))
    {
        return Err(GrError::InvalidDomain(
            "fold_eval requires fiber points with unit differences",
        ));
    }

    let denominators = [
        ctx.mul(&d01, &d02),
        ctx.mul(&d10, &d12),
        ctx.mul(&d20, &d21),
    ];
    let denominator_inverses = ctx.batch_inv(&denominators)?;
    let a0 = ctx.sub(alpha, p0);
    let a1 = ctx.sub(alpha, p1);
    let a2 = ctx.sub(alpha, p2);
    let numerators = [ctx.mul(&a1, &a2), ctx.mul(&a0, &a2), ctx.mul(&a0, &a1)];

    let mut out = ctx.zero();
    for index in 0..3 {
        let basis_value = ctx.mul(&numerators[index], &denominator_inverses[index]);
        out = ctx.add(&out, &ctx.mul(&fiber_values[index], &basis_value));
    }
    Ok(out)
}

pub fn ternary_fold_table(
    domain: &Domain,
    evals: &[GrElem],
    alpha: &GrElem,
) -> Result<Vec<GrElem>> {
    if evals.len() != checked_size(domain.size(), "domain size")? {
        return Err(GrError::InvalidDomain(
            "ternary_fold_table requires eval count == domain size",
        ));
    }
    if !domain.size().is_multiple_of(3) {
        return Err(GrError::InvalidDomain(
            "ternary_fold_table requires domain size divisible by 3",
        ));
    }

    let ctx = domain.context();
    let folded_size = domain.size() / 3;
    let folded_size_usize = checked_size(folded_size, "folded size")?;
    let omega = ctx.pow(domain.root(), folded_size.into());
    let omega_squared = ctx.square(&omega);
    let one = ctx.one();
    let denominator_constants = [
        ctx.mul(&ctx.sub(&one, &omega), &ctx.sub(&one, &omega_squared)),
        ctx.mul(&ctx.sub(&omega, &one), &ctx.sub(&omega, &omega_squared)),
        ctx.mul(
            &ctx.sub(&omega_squared, &one),
            &ctx.sub(&omega_squared, &omega),
        ),
    ];
    let denominator_constant_inverses = ctx.batch_inv(&denominator_constants)?;
    let base_points = domain
        .iter_elements()
        .take(folded_size_usize)
        .collect::<Vec<_>>();
    let base_inverses = ctx.batch_inv(&base_points)?;

    #[cfg(feature = "parallel")]
    {
        if should_parallel_fold(folded_size_usize) {
            return (0..folded_size_usize)
                .into_par_iter()
                .map(|base| {
                    let mut scratch = GrScratch::with_elements(ctx, TERNARY_FOLD_SCRATCH_ELEMS);
                    Ok(fold_ordered_ternary_table_value_with_scratch(
                        ctx,
                        TernaryFoldTableInput {
                            base_point: &base_points[base],
                            base_inverse: &base_inverses[base],
                            evals,
                            base,
                            folded_size: folded_size_usize,
                            omega: &omega,
                            omega_squared: &omega_squared,
                            denominator_constant_inverses: &denominator_constant_inverses,
                            alpha,
                        },
                        &mut scratch,
                    ))
                })
                .collect();
        }
    }

    let mut scratch = GrScratch::with_elements(ctx, TERNARY_FOLD_SCRATCH_ELEMS);
    let mut out = Vec::with_capacity(folded_size_usize);
    for (base, base_point) in base_points.iter().enumerate() {
        out.push(fold_ordered_ternary_table_value_with_scratch(
            ctx,
            TernaryFoldTableInput {
                base_point,
                base_inverse: &base_inverses[base],
                evals,
                base,
                folded_size: folded_size_usize,
                omega: &omega,
                omega_squared: &omega_squared,
                denominator_constant_inverses: &denominator_constant_inverses,
                alpha,
            },
            &mut scratch,
        ));
    }
    Ok(out)
}

pub fn repeated_ternary_fold_table(
    domain: &Domain,
    evals: &[GrElem],
    alphas: &[GrElem],
) -> Result<Vec<GrElem>> {
    if evals.len() != checked_size(domain.size(), "domain size")? {
        return Err(GrError::InvalidDomain(
            "repeated_ternary_fold_table requires eval count == domain size",
        ));
    }

    let mut current_domain = domain.clone();
    let mut current_evals = evals.to_vec();
    for alpha in alphas {
        current_evals = ternary_fold_table(&current_domain, &current_evals, alpha)?;
        current_domain = current_domain.pow_map(3)?;
    }
    Ok(current_evals)
}

pub fn virtual_fold_query_indices(domain_size: u64, b: u64, child_index: u64) -> Result<Vec<u64>> {
    if domain_size == 0 {
        return Err(GrError::InvalidDomain(
            "virtual_fold_query_indices requires domain_size > 0",
        ));
    }

    let fold_width = pow3_checked(b, "virtual fold width")?;
    if !domain_size.is_multiple_of(fold_width) {
        return Err(GrError::InvalidDomain(
            "virtual_fold_query_indices requires 3^b dividing domain_size",
        ));
    }

    let child_count = domain_size / fold_width;
    if child_index >= child_count {
        return Err(GrError::IndexOutOfRange {
            index: child_index,
            size: child_count,
        });
    }

    build_virtual_fold_query_indices(domain_size, b, child_index)
}

pub(crate) fn virtual_fold_query_points(
    domain: &Domain,
    b: u64,
    child_index: u64,
) -> Result<Vec<GrElem>> {
    let fold_width = pow3_checked(b, "virtual fold width")?;
    if !domain.size().is_multiple_of(fold_width) {
        return Err(GrError::InvalidDomain(
            "virtual_fold_query_points requires 3^b dividing domain size",
        ));
    }
    let child_count = domain.size() / fold_width;
    if child_index >= child_count {
        return Err(GrError::IndexOutOfRange {
            index: child_index,
            size: child_count,
        });
    }

    build_virtual_fold_query_points(domain, b, child_index)
}

pub fn evaluate_repeated_ternary_fold_from_values(
    ctx: &GrContext,
    points: &[GrElem],
    values: &[GrElem],
    alphas: &[GrElem],
) -> Result<GrElem> {
    evaluate_repeated_ternary_fold_from_values_impl(ctx, points, values, alphas, true)
}

#[cfg(test)]
pub(crate) fn evaluate_ordered_repeated_ternary_fold_from_values(
    ctx: &GrContext,
    points: &[GrElem],
    values: &[GrElem],
    alphas: &[GrElem],
) -> Result<GrElem> {
    evaluate_repeated_ternary_fold_from_values_impl(ctx, points, values, alphas, false)
}

pub(crate) fn evaluate_ordered_repeated_ternary_fold_batch_from_values(
    ctx: &GrContext,
    mut point_batches: Vec<Vec<GrElem>>,
    mut value_batches: Vec<Vec<GrElem>>,
    alphas: &[GrElem],
) -> Result<Vec<GrElem>> {
    if point_batches.len() != value_batches.len() {
        return Err(GrError::InvalidDomain(
            "batched repeated ternary fold requires equal batch counts",
        ));
    }

    let expected_size = pow3_checked(alphas.len() as u64, "repeated ternary sparse input size")?;
    let expected_len = checked_size(expected_size, "sparse input size")?;
    for (points, values) in point_batches.iter().zip(&value_batches) {
        if points.len() != values.len() {
            return Err(GrError::InvalidDomain(
                "batched repeated ternary fold requires equal-sized query inputs",
            ));
        }
        if points.len() != expected_len {
            return Err(GrError::InvalidDomain(
                "batched repeated ternary fold requires 3^b values per query",
            ));
        }
        if values.is_empty() {
            return Err(GrError::InvalidDomain(
                "batched repeated ternary fold requires non-empty query inputs",
            ));
        }
    }

    for alpha in alphas {
        let (next_points, next_values) =
            fold_ordered_ternary_level_fast_batches(ctx, &point_batches, &value_batches, alpha)?;
        point_batches = next_points;
        value_batches = next_values;
    }

    value_batches
        .into_iter()
        .map(|mut values| {
            if values.len() != 1 {
                return Err(GrError::InvalidDomain(
                    "batched fold result did not collapse to one value",
                ));
            }
            Ok(values.remove(0))
        })
        .collect()
}

fn evaluate_repeated_ternary_fold_from_values_impl(
    ctx: &GrContext,
    points: &[GrElem],
    values: &[GrElem],
    alphas: &[GrElem],
    validate_cube_fibers: bool,
) -> Result<GrElem> {
    if points.len() != values.len() {
        return Err(GrError::InvalidDomain(
            "evaluate_repeated_ternary_fold_from_values requires equal-sized inputs",
        ));
    }

    let expected_size = pow3_checked(alphas.len() as u64, "repeated ternary sparse input size")?;
    if points.len() != checked_size(expected_size, "sparse input size")? {
        return Err(GrError::InvalidDomain(
            "evaluate_repeated_ternary_fold_from_values requires 3^b values",
        ));
    }
    if values.is_empty() {
        return Err(GrError::InvalidDomain(
            "evaluate_repeated_ternary_fold_from_values requires non-empty inputs",
        ));
    }

    let mut current_points = points.to_vec();
    let mut current_values = values.to_vec();

    for alpha in alphas {
        let (next_points, next_values) = fold_ternary_level_batched(
            ctx,
            &current_points,
            &current_values,
            alpha,
            validate_cube_fibers,
        )?;
        current_points = next_points;
        current_values = next_values;
    }

    current_values
        .into_iter()
        .next()
        .ok_or(GrError::InvalidDomain("fold result is unexpectedly empty"))
}

fn fold_ternary_level_batched(
    ctx: &GrContext,
    points: &[GrElem],
    values: &[GrElem],
    alpha: &GrElem,
    validate_cube_fibers: bool,
) -> Result<(Vec<GrElem>, Vec<GrElem>)> {
    if points.len() != values.len() {
        return Err(GrError::InvalidDomain(
            "fold_ternary_level_batched requires equal-sized inputs",
        ));
    }
    if !points.len().is_multiple_of(3) {
        return Err(GrError::InvalidDomain(
            "fold_ternary_level_batched saw non-ternary level",
        ));
    }
    if !validate_cube_fibers {
        return fold_ordered_ternary_level_fast(ctx, points, values, alpha);
    }

    let next_size = points.len() / 3;
    let mut next_points = Vec::with_capacity(next_size);
    let mut denominator_products = Vec::with_capacity(points.len());
    for group in 0..next_size {
        let base = group * 3;
        let fiber_points = &points[base..base + 3];
        if validate_cube_fibers {
            require_cube_fiber(ctx, fiber_points)?;
        }
        next_points.push(ctx.pow(&fiber_points[0], 3));

        let p0 = &fiber_points[0];
        let p1 = &fiber_points[1];
        let p2 = &fiber_points[2];
        let d01 = ctx.sub(p0, p1);
        let d02 = ctx.sub(p0, p2);
        let d10 = ctx.sub(p1, p0);
        let d12 = ctx.sub(p1, p2);
        let d20 = ctx.sub(p2, p0);
        let d21 = ctx.sub(p2, p1);
        denominator_products.push(ctx.mul(&d01, &d02));
        denominator_products.push(ctx.mul(&d10, &d12));
        denominator_products.push(ctx.mul(&d20, &d21));
    }
    let denominator_inverses = ctx.batch_inv(&denominator_products)?;

    let mut next_values = Vec::with_capacity(next_size);
    for group in 0..next_size {
        let base = group * 3;
        let p0 = &points[base];
        let p1 = &points[base + 1];
        let p2 = &points[base + 2];
        let a0 = ctx.sub(alpha, p0);
        let a1 = ctx.sub(alpha, p1);
        let a2 = ctx.sub(alpha, p2);
        let numerators = [ctx.mul(&a1, &a2), ctx.mul(&a0, &a2), ctx.mul(&a0, &a1)];

        let mut out = ctx.zero();
        for index in 0..3 {
            let basis_value = ctx.mul(&numerators[index], &denominator_inverses[base + index]);
            out = ctx.add(&out, &ctx.mul(&values[base + index], &basis_value));
        }
        next_values.push(out);
    }
    Ok((next_points, next_values))
}

fn fold_ordered_ternary_level_fast(
    ctx: &GrContext,
    points: &[GrElem],
    values: &[GrElem],
    alpha: &GrElem,
) -> Result<(Vec<GrElem>, Vec<GrElem>)> {
    let next_size = points.len() / 3;
    let mut next_points = Vec::with_capacity(next_size);
    let mut base_points = Vec::with_capacity(next_size);
    for group in 0..next_size {
        let base = group * 3;
        next_points.push(ctx.pow(&points[base], 3));
        base_points.push(points[base].clone());
    }
    let base_inverses = ctx.batch_inv(&base_points)?;

    let omega = ctx.mul(&points[1], &base_inverses[0]);
    let omega_squared = ctx.square(&omega);
    let one = ctx.one();
    let denominator_constants = [
        ctx.mul(&ctx.sub(&one, &omega), &ctx.sub(&one, &omega_squared)),
        ctx.mul(&ctx.sub(&omega, &one), &ctx.sub(&omega, &omega_squared)),
        ctx.mul(
            &ctx.sub(&omega_squared, &one),
            &ctx.sub(&omega_squared, &omega),
        ),
    ];
    let denominator_constant_inverses = ctx.batch_inv(&denominator_constants)?;

    let mut scratch = GrScratch::with_elements(ctx, TERNARY_FOLD_SCRATCH_ELEMS);
    let mut next_values = Vec::with_capacity(next_size);
    for (group, base_inverse) in base_inverses.iter().enumerate() {
        let base = group * 3;
        next_values.push(fold_ordered_ternary_value_with_scratch(
            ctx,
            TernaryFoldValueInput {
                points: [&points[base], &points[base + 1], &points[base + 2]],
                values: [&values[base], &values[base + 1], &values[base + 2]],
                base_inverse,
                denominator_constant_inverses: &denominator_constant_inverses,
                alpha,
            },
            &mut scratch,
        ));
    }
    Ok((next_points, next_values))
}

fn fold_ordered_ternary_level_fast_batches(
    ctx: &GrContext,
    point_batches: &[Vec<GrElem>],
    value_batches: &[Vec<GrElem>],
    alpha: &GrElem,
) -> Result<(FoldQueryBatches, FoldQueryBatches)> {
    if point_batches.len() != value_batches.len() {
        return Err(GrError::InvalidDomain(
            "fold_ordered_ternary_level_fast_batches requires equal batch counts",
        ));
    }
    if point_batches.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }

    let level_size = point_batches[0].len();
    if level_size == 0 || !level_size.is_multiple_of(3) {
        return Err(GrError::InvalidDomain(
            "fold_ordered_ternary_level_fast_batches saw non-ternary level",
        ));
    }
    for (points, values) in point_batches.iter().zip(value_batches) {
        if points.len() != level_size || values.len() != level_size {
            return Err(GrError::InvalidDomain(
                "fold_ordered_ternary_level_fast_batches requires uniform levels",
            ));
        }
    }

    let next_size = level_size / 3;
    let total_groups = next_size
        .checked_mul(point_batches.len())
        .ok_or(GrError::ArithmeticOverflow("batched fold group count"))?;
    let mut base_points = Vec::with_capacity(total_groups);
    for points in point_batches {
        for group in 0..next_size {
            base_points.push(points[group * 3].clone());
        }
    }
    let base_inverses = ctx.batch_inv(&base_points)?;

    let omega = ctx.mul(&point_batches[0][1], &base_inverses[0]);
    let omega_squared = ctx.square(&omega);
    let one = ctx.one();
    let denominator_constants = [
        ctx.mul(&ctx.sub(&one, &omega), &ctx.sub(&one, &omega_squared)),
        ctx.mul(&ctx.sub(&omega, &one), &ctx.sub(&omega, &omega_squared)),
        ctx.mul(
            &ctx.sub(&omega_squared, &one),
            &ctx.sub(&omega_squared, &omega),
        ),
    ];
    let denominator_constant_inverses = ctx.batch_inv(&denominator_constants)?;

    #[cfg(feature = "parallel")]
    {
        if should_parallel_fold(total_groups) {
            return fold_ordered_ternary_level_fast_batches_parallel(
                ctx,
                point_batches,
                value_batches,
                alpha,
                next_size,
                &base_inverses,
                &denominator_constant_inverses,
            );
        }
    }

    let mut inverse_index = 0;
    let mut scratch = GrScratch::with_elements(ctx, TERNARY_FOLD_SCRATCH_ELEMS);
    let mut next_point_batches = Vec::with_capacity(point_batches.len());
    let mut next_value_batches = Vec::with_capacity(value_batches.len());
    for (points, values) in point_batches.iter().zip(value_batches) {
        let mut next_points = Vec::with_capacity(next_size);
        let mut next_values = Vec::with_capacity(next_size);
        for group in 0..next_size {
            let base = group * 3;
            next_points.push(cube_with_scratch(ctx, &points[base], &mut scratch));

            let base_inverse = &base_inverses[inverse_index];
            inverse_index += 1;
            next_values.push(fold_ordered_ternary_value_with_scratch(
                ctx,
                TernaryFoldValueInput {
                    points: [&points[base], &points[base + 1], &points[base + 2]],
                    values: [&values[base], &values[base + 1], &values[base + 2]],
                    base_inverse,
                    denominator_constant_inverses: &denominator_constant_inverses,
                    alpha,
                },
                &mut scratch,
            ));
        }
        next_point_batches.push(next_points);
        next_value_batches.push(next_values);
    }

    Ok((next_point_batches, next_value_batches))
}

#[cfg(feature = "parallel")]
fn fold_ordered_ternary_level_fast_batches_parallel(
    ctx: &GrContext,
    point_batches: &[Vec<GrElem>],
    value_batches: &[Vec<GrElem>],
    alpha: &GrElem,
    next_size: usize,
    base_inverses: &[GrElem],
    denominator_constant_inverses: &[GrElem],
) -> Result<(FoldQueryBatches, FoldQueryBatches)> {
    let folded = point_batches
        .par_iter()
        .zip(value_batches)
        .enumerate()
        .map(|(batch_index, (points, values))| {
            #[cfg(feature = "parallel")]
            {
                if should_parallel_fold(next_size) {
                    let folded = (0..next_size)
                        .into_par_iter()
                        .map(|group| {
                            let mut scratch =
                                GrScratch::with_elements(ctx, TERNARY_FOLD_SCRATCH_ELEMS);
                            let base = group * 3;
                            (
                                cube_with_scratch(ctx, &points[base], &mut scratch),
                                fold_ordered_ternary_value_with_scratch(
                                    ctx,
                                    TernaryFoldValueInput {
                                        points: [
                                            &points[base],
                                            &points[base + 1],
                                            &points[base + 2],
                                        ],
                                        values: [
                                            &values[base],
                                            &values[base + 1],
                                            &values[base + 2],
                                        ],
                                        base_inverse: &base_inverses
                                            [batch_index * next_size + group],
                                        denominator_constant_inverses,
                                        alpha,
                                    },
                                    &mut scratch,
                                ),
                            )
                        })
                        .collect::<Vec<_>>();
                    let mut next_points = Vec::with_capacity(next_size);
                    let mut next_values = Vec::with_capacity(next_size);
                    for (next_point, next_value) in folded {
                        next_points.push(next_point);
                        next_values.push(next_value);
                    }
                    return Ok((next_points, next_values));
                }
            }

            let mut scratch = GrScratch::with_elements(ctx, TERNARY_FOLD_SCRATCH_ELEMS);
            let mut next_points = Vec::with_capacity(next_size);
            let mut next_values = Vec::with_capacity(next_size);
            for group in 0..next_size {
                let base = group * 3;
                next_points.push(cube_with_scratch(ctx, &points[base], &mut scratch));
                next_values.push(fold_ordered_ternary_value_with_scratch(
                    ctx,
                    TernaryFoldValueInput {
                        points: [&points[base], &points[base + 1], &points[base + 2]],
                        values: [&values[base], &values[base + 1], &values[base + 2]],
                        base_inverse: &base_inverses[batch_index * next_size + group],
                        denominator_constant_inverses,
                        alpha,
                    },
                    &mut scratch,
                ));
            }
            Ok((next_points, next_values))
        })
        .collect::<Result<Vec<_>>>()?;

    let mut next_point_batches = Vec::with_capacity(folded.len());
    let mut next_value_batches = Vec::with_capacity(folded.len());
    for (next_points, next_values) in folded {
        next_point_batches.push(next_points);
        next_value_batches.push(next_values);
    }
    Ok((next_point_batches, next_value_batches))
}

#[derive(Clone, Copy)]
struct TernaryFoldTableInput<'a> {
    base_point: &'a GrElem,
    base_inverse: &'a GrElem,
    evals: &'a [GrElem],
    base: usize,
    folded_size: usize,
    omega: &'a GrElem,
    omega_squared: &'a GrElem,
    denominator_constant_inverses: &'a [GrElem],
    alpha: &'a GrElem,
}

#[derive(Clone, Copy)]
struct TernaryFoldValueInput<'a> {
    points: [&'a GrElem; 3],
    values: [&'a GrElem; 3],
    base_inverse: &'a GrElem,
    denominator_constant_inverses: &'a [GrElem],
    alpha: &'a GrElem,
}

fn fold_ordered_ternary_table_value_with_scratch(
    ctx: &GrContext,
    input: TernaryFoldTableInput<'_>,
    scratch: &mut GrScratch,
) -> GrElem {
    let p1 = ctx.mul(input.base_point, input.omega);
    let p2 = ctx.mul(input.base_point, input.omega_squared);
    fold_ordered_ternary_value_with_scratch(
        ctx,
        TernaryFoldValueInput {
            points: [input.base_point, &p1, &p2],
            values: [
                &input.evals[input.base],
                &input.evals[input.base + input.folded_size],
                &input.evals[input.base + 2 * input.folded_size],
            ],
            base_inverse: input.base_inverse,
            denominator_constant_inverses: input.denominator_constant_inverses,
            alpha: input.alpha,
        },
        scratch,
    )
}

fn cube_with_scratch(ctx: &GrContext, value: &GrElem, scratch: &mut GrScratch) -> GrElem {
    let (elements, mul_scratch) = scratch.parts_mut();
    let [square, cube, ..] = elements else {
        unreachable!("ternary fold scratch has at least two elements");
    };
    ctx.square_into(square, value, mul_scratch);
    ctx.mul_into(cube, square, value, mul_scratch);
    cube.clone()
}

fn fold_ordered_ternary_value_with_scratch(
    ctx: &GrContext,
    input: TernaryFoldValueInput<'_>,
    scratch: &mut GrScratch,
) -> GrElem {
    let (elements, mul_scratch) = scratch.parts_mut();
    let [a0, a1, a2, numerator, inv_base_squared, denominator_inverse, basis_value, product, sum, sum_next, ..] =
        elements
    else {
        unreachable!("ternary fold scratch has enough elements");
    };

    ctx.square_into(inv_base_squared, input.base_inverse, mul_scratch);
    ctx.sub_into(a0, input.alpha, input.points[0]);
    ctx.sub_into(a1, input.alpha, input.points[1]);
    ctx.sub_into(a2, input.alpha, input.points[2]);

    clear_elem(sum);
    for (index, denominator_constant_inverse) in input
        .denominator_constant_inverses
        .iter()
        .enumerate()
        .take(3)
    {
        match index {
            0 => ctx.mul_into(numerator, a1, a2, mul_scratch),
            1 => ctx.mul_into(numerator, a0, a2, mul_scratch),
            _ => ctx.mul_into(numerator, a0, a1, mul_scratch),
        }
        ctx.mul_into(
            denominator_inverse,
            inv_base_squared,
            denominator_constant_inverse,
            mul_scratch,
        );
        ctx.mul_into(basis_value, numerator, denominator_inverse, mul_scratch);
        ctx.mul_into(product, input.values[index], basis_value, mul_scratch);
        ctx.add_into(sum_next, sum, product);
        std::mem::swap(sum, sum_next);
    }
    sum.clone()
}

#[cfg(feature = "parallel")]
fn should_parallel_fold(total_groups: usize) -> bool {
    total_groups >= PARALLEL_FOLD_GROUP_THRESHOLD && rayon::current_num_threads() > 1
}

pub fn evaluate_virtual_fold_query_from_leaf_payloads(
    domain: &Domain,
    b: u64,
    child_index: u64,
    leaf_payloads: &[Vec<u8>],
    alphas: &[GrElem],
) -> Result<GrElem> {
    if alphas.len() != checked_size(b, "fold depth")? {
        return Err(GrError::InvalidDomain(
            "evaluate_virtual_fold_query_from_leaf_payloads requires one alpha per level",
        ));
    }

    let parent_indices = virtual_fold_query_indices(domain.size(), b, child_index)?;
    if leaf_payloads.len() != parent_indices.len() {
        return Err(GrError::InvalidDomain(
            "evaluate_virtual_fold_query_from_leaf_payloads requires one payload per parent index",
        ));
    }

    let ctx = domain.context();
    let mut points = Vec::with_capacity(parent_indices.len());
    let mut values = Vec::with_capacity(parent_indices.len());
    for (index, payload) in parent_indices.iter().zip(leaf_payloads) {
        points.push(domain.element(*index)?);
        values.push(ctx.deserialize(payload)?);
    }

    evaluate_repeated_ternary_fold_from_values(ctx, &points, &values, alphas)
}

fn require_cube_fiber(ctx: &GrContext, fiber_points: &[GrElem]) -> Result<()> {
    if fiber_points.len() != 3 {
        return Err(GrError::InvalidDomain(
            "ternary fold requires exactly three fiber points",
        ));
    }
    for left in 0..fiber_points.len() {
        for right in left + 1..fiber_points.len() {
            if fiber_points[left] == fiber_points[right] {
                return Err(GrError::InvalidDomain(
                    "ternary fold requires distinct fiber points",
                ));
            }
            if !ctx.is_unit(&ctx.sub(&fiber_points[left], &fiber_points[right])) {
                return Err(GrError::InvalidDomain(
                    "ternary fold requires unit fiber point differences",
                ));
            }
        }
    }

    let mapped_point = ctx.pow(&fiber_points[0], 3);
    if fiber_points
        .iter()
        .skip(1)
        .any(|point| ctx.pow(point, 3) != mapped_point)
    {
        return Err(GrError::InvalidDomain("ternary fold requires cube fibers"));
    }
    Ok(())
}

fn build_virtual_fold_query_indices(
    domain_size: u64,
    b: u64,
    child_index: u64,
) -> Result<Vec<u64>> {
    if b == 0 {
        return Ok(vec![child_index]);
    }

    let next_domain_size = domain_size / 3;
    let next_indices = build_virtual_fold_query_indices(next_domain_size, b - 1, child_index)?;
    let mut out = Vec::with_capacity(next_indices.len() * 3);
    for next_index in next_indices {
        for offset in 0..3 {
            out.push(next_index + offset * next_domain_size);
        }
    }
    Ok(out)
}

fn build_virtual_fold_query_points(
    domain: &Domain,
    b: u64,
    child_index: u64,
) -> Result<Vec<GrElem>> {
    if b == 0 {
        return Ok(vec![domain.element(child_index)?]);
    }

    let next_domain_size = domain.size() / 3;
    let next_indices = build_virtual_fold_query_indices(next_domain_size, b - 1, child_index)?;
    let ctx = domain.context();
    let omega = ctx.pow(domain.root(), next_domain_size.into());
    let omega_squared = ctx.square(&omega);
    let mut out = Vec::with_capacity(next_indices.len() * 3);
    for next_index in next_indices {
        let base = domain.element(next_index)?;
        out.push(base.clone());
        out.push(ctx.mul(&base, &omega));
        out.push(ctx.mul(&base, &omega_squared));
    }
    Ok(out)
}

fn pow3_checked(exponent: u64, label: &'static str) -> Result<u64> {
    let mut out = 1u64;
    for _ in 0..exponent {
        out = out
            .checked_mul(3)
            .ok_or(GrError::ArithmeticOverflow(label))?;
    }
    Ok(out)
}

fn checked_size(value: u64, label: &'static str) -> Result<usize> {
    usize::try_from(value).map_err(|_| GrError::ArithmeticOverflow(label))
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::{
        algebra::galois_ring::{Domain, GrConfig, GrContext, GrError},
        protocols::whir_gr::folding::{
            evaluate_ordered_repeated_ternary_fold_batch_from_values,
            evaluate_ordered_repeated_ternary_fold_from_values,
            evaluate_repeated_ternary_fold_from_values,
            evaluate_virtual_fold_query_from_leaf_payloads, fold_eval, repeated_ternary_fold_table,
            virtual_fold_query_indices, virtual_fold_query_points,
        },
    };

    fn sample_context(r: usize) -> Arc<GrContext> {
        Arc::new(GrContext::new(GrConfig { p: 2, k_exp: 16, r }).unwrap())
    }

    fn sample_polynomial_coefficients(
        ctx: &GrContext,
        domain: &Domain,
    ) -> Vec<crate::algebra::galois_ring::GrElem> {
        let mut coefficients = Vec::with_capacity(domain.size() as usize);
        let mut root_power = ctx.one();
        let mut offset_power = ctx.one();
        let twist = domain.element(1).unwrap();
        let mut twist_power = ctx.one();
        for index in 0..domain.size() {
            let mut value = ctx.add(&root_power, &offset_power);
            if index % 3 == 2 {
                value = ctx.add(&value, &twist_power);
            }
            coefficients.push(value);
            root_power = ctx.mul(&root_power, domain.root());
            offset_power = ctx.mul(&offset_power, domain.offset());
            twist_power = ctx.mul(&twist_power, &twist);
        }
        let last = coefficients.last_mut().unwrap();
        *last = ctx.add(last, domain.root());
        coefficients
    }

    fn evaluate_polynomial(
        ctx: &GrContext,
        coefficients: &[crate::algebra::galois_ring::GrElem],
        point: &crate::algebra::galois_ring::GrElem,
    ) -> crate::algebra::galois_ring::GrElem {
        let mut out = ctx.zero();
        for coefficient in coefficients.iter().rev() {
            out = ctx.add(&ctx.mul(&out, point), coefficient);
        }
        out
    }

    fn sample_evals(ctx: &GrContext, domain: &Domain) -> Vec<crate::algebra::galois_ring::GrElem> {
        let coefficients = sample_polynomial_coefficients(ctx, domain);
        (0..domain.size())
            .map(|index| evaluate_polynomial(ctx, &coefficients, &domain.element(index).unwrap()))
            .collect()
    }

    fn sparse_points(domain: &Domain, indices: &[u64]) -> Vec<crate::algebra::galois_ring::GrElem> {
        indices
            .iter()
            .map(|&index| domain.element(index).unwrap())
            .collect()
    }

    fn sparse_values(
        evals: &[crate::algebra::galois_ring::GrElem],
        indices: &[u64],
    ) -> Vec<crate::algebra::galois_ring::GrElem> {
        indices
            .iter()
            .map(|&index| evals[index as usize].clone())
            .collect()
    }

    fn sparse_payloads(
        ctx: &GrContext,
        evals: &[crate::algebra::galois_ring::GrElem],
        indices: &[u64],
    ) -> Vec<Vec<u8>> {
        indices
            .iter()
            .map(|&index| ctx.serialize(&evals[index as usize]))
            .collect()
    }

    fn check_sparse_matches_full_table(
        domain: &Domain,
        alphas: &[crate::algebra::galois_ring::GrElem],
    ) {
        let ctx = domain.context();
        let evals = sample_evals(ctx, domain);
        let full = repeated_ternary_fold_table(domain, &evals, alphas).unwrap();
        let mut batch_points = Vec::new();
        let mut batch_values = Vec::new();
        let mut batch_expected = Vec::new();

        for child in 0..full.len() as u64 {
            let indices =
                virtual_fold_query_indices(domain.size(), alphas.len() as u64, child).unwrap();
            let points = sparse_points(domain, &indices);
            assert_eq!(
                virtual_fold_query_points(domain, alphas.len() as u64, child).unwrap(),
                points
            );
            let values = sparse_values(&evals, &indices);
            batch_points.push(points.clone());
            batch_values.push(values.clone());
            batch_expected.push(full[child as usize].clone());

            let sparse =
                evaluate_repeated_ternary_fold_from_values(ctx, &points, &values, alphas).unwrap();
            assert_eq!(sparse, full[child as usize]);
            let ordered_fast =
                evaluate_ordered_repeated_ternary_fold_from_values(ctx, &points, &values, alphas)
                    .unwrap();
            assert_eq!(ordered_fast, sparse);

            let payloads = sparse_payloads(ctx, &evals, &indices);
            let from_payloads = evaluate_virtual_fold_query_from_leaf_payloads(
                domain,
                alphas.len() as u64,
                child,
                &payloads,
                alphas,
            )
            .unwrap();
            assert_eq!(from_payloads, full[child as usize]);
        }

        let batched = evaluate_ordered_repeated_ternary_fold_batch_from_values(
            ctx,
            batch_points,
            batch_values,
            alphas,
        )
        .unwrap();
        assert_eq!(batched, batch_expected);
    }

    #[test]
    fn repeated_ternary_sparse_b1_should_match_full_table() {
        let ctx = sample_context(6);
        let domain = Domain::teichmuller_subgroup(Arc::clone(&ctx), 9).unwrap();
        let alphas = vec![ctx.add(&ctx.one(), &domain.element(2).unwrap())];

        check_sparse_matches_full_table(&domain, &alphas);
    }

    #[test]
    fn repeated_ternary_sparse_b2_should_match_full_table() {
        let ctx = sample_context(6);
        let domain = Domain::teichmuller_subgroup(Arc::clone(&ctx), 9).unwrap();
        let alphas = vec![
            ctx.add(&ctx.one(), &domain.element(2).unwrap()),
            ctx.add(&ctx.one(), &domain.element(4).unwrap()),
        ];

        check_sparse_matches_full_table(&domain, &alphas);
    }

    #[test]
    fn repeated_ternary_sparse_b3_should_match_full_table() {
        let ctx = sample_context(18);
        let domain = Domain::teichmuller_subgroup(Arc::clone(&ctx), 27).unwrap();
        let alphas = vec![
            ctx.add(&ctx.one(), &domain.element(2).unwrap()),
            ctx.add(&ctx.one(), &domain.element(4).unwrap()),
            ctx.add(&ctx.one(), &domain.element(8).unwrap()),
        ];

        check_sparse_matches_full_table(&domain, &alphas);
    }

    #[test]
    fn virtual_parent_index_shape_should_match_recursive_fibers() {
        let ctx = sample_context(18);
        let domain = Domain::teichmuller_subgroup(Arc::clone(&ctx), 27).unwrap();
        let folded = domain.pow_map(9).unwrap();

        let indices = virtual_fold_query_indices(domain.size(), 2, 1).unwrap();
        let expected = vec![1, 10, 19, 4, 13, 22, 7, 16, 25];

        assert_eq!(indices, expected);
        let folded_point = folded.element(1).unwrap();
        for index in indices {
            assert_eq!(ctx.pow(&domain.element(index).unwrap(), 9), folded_point);
        }
    }

    #[test]
    fn fold_eval_should_interpolate_at_alpha() {
        let ctx = sample_context(6);
        let domain = Domain::teichmuller_subgroup(Arc::clone(&ctx), 9).unwrap();
        let points = vec![
            domain.element(0).unwrap(),
            domain.element(3).unwrap(),
            domain.element(6).unwrap(),
        ];
        let values = points
            .iter()
            .map(|point| ctx.add(&ctx.from_u64(7), &ctx.mul(&ctx.from_u64(5), point)))
            .collect::<Vec<_>>();
        let alpha = ctx.add(&ctx.one(), &domain.element(2).unwrap());

        let folded = fold_eval(&ctx, &points, &values, &alpha).unwrap();

        assert_eq!(
            folded,
            ctx.add(&ctx.from_u64(7), &ctx.mul(&ctx.from_u64(5), &alpha))
        );
    }

    #[test]
    fn invalid_inputs_should_reject() {
        let ctx = sample_context(6);
        let domain = Domain::teichmuller_subgroup(Arc::clone(&ctx), 9).unwrap();
        let alphas = vec![ctx.add(&ctx.one(), &domain.element(1).unwrap())];
        let short_evals = vec![ctx.one(); 8];
        let short_points = vec![domain.element(0).unwrap(), domain.element(3).unwrap()];
        let short_values = vec![ctx.one(), domain.root().clone()];
        let duplicate_points = vec![ctx.one(), ctx.one(), domain.root().clone()];
        let duplicate_values = vec![ctx.one(), domain.root().clone(), ctx.one()];

        assert!(matches!(
            repeated_ternary_fold_table(&domain, &short_evals, &alphas),
            Err(GrError::InvalidDomain(_))
        ));
        assert!(matches!(
            virtual_fold_query_indices(10, 2, 0),
            Err(GrError::InvalidDomain(_))
        ));
        assert!(matches!(
            virtual_fold_query_indices(9, 2, 1),
            Err(GrError::IndexOutOfRange { .. })
        ));
        assert!(matches!(
            evaluate_repeated_ternary_fold_from_values(&ctx, &short_points, &short_values, &alphas,),
            Err(GrError::InvalidDomain(_))
        ));
        assert!(matches!(
            evaluate_repeated_ternary_fold_from_values(
                &ctx,
                &duplicate_points,
                &duplicate_values,
                &alphas,
            ),
            Err(GrError::InvalidDomain(_))
        ));
        assert!(matches!(
            evaluate_virtual_fold_query_from_leaf_payloads(&domain, 1, 0, &[], &alphas),
            Err(GrError::InvalidDomain(_))
        ));
    }
}
