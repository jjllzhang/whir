use crate::algebra::galois_ring::{Domain, GrContext, GrElem, GrError, Result};

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
    let mut out = Vec::with_capacity(checked_size(folded_size, "folded size")?);
    for base in 0..folded_size {
        let mut points = Vec::with_capacity(3);
        let mut values = Vec::with_capacity(3);
        for offset in 0..3 {
            let index = base + offset * folded_size;
            points.push(domain.element(index)?);
            values.push(evals[checked_size(index, "fold index")?].clone());
        }
        out.push(fold_eval(ctx, &points, &values, alpha)?);
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

pub fn evaluate_repeated_ternary_fold_from_values(
    ctx: &GrContext,
    points: &[GrElem],
    values: &[GrElem],
    alphas: &[GrElem],
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
        if !current_points.len().is_multiple_of(3) {
            return Err(GrError::InvalidDomain(
                "evaluate_repeated_ternary_fold_from_values saw non-ternary level",
            ));
        }

        let next_size = current_points.len() / 3;
        let mut next_points = Vec::with_capacity(next_size);
        let mut next_values = Vec::with_capacity(next_size);
        for group in 0..next_size {
            let base = group * 3;
            let fiber_points = &current_points[base..base + 3];
            let fiber_values = &current_values[base..base + 3];

            require_cube_fiber(ctx, fiber_points)?;
            next_points.push(ctx.pow(&fiber_points[0], 3));
            next_values.push(fold_eval(ctx, fiber_points, fiber_values, alpha)?);
        }

        current_points = next_points;
        current_values = next_values;
    }

    current_values
        .into_iter()
        .next()
        .ok_or(GrError::InvalidDomain("fold result is unexpectedly empty"))
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
            evaluate_repeated_ternary_fold_from_values,
            evaluate_virtual_fold_query_from_leaf_payloads, fold_eval, repeated_ternary_fold_table,
            virtual_fold_query_indices,
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

        for child in 0..full.len() as u64 {
            let indices =
                virtual_fold_query_indices(domain.size(), alphas.len() as u64, child).unwrap();
            let points = sparse_points(domain, &indices);
            let values = sparse_values(&evals, &indices);

            let sparse =
                evaluate_repeated_ternary_fold_from_values(ctx, &points, &values, alphas).unwrap();
            assert_eq!(sparse, full[child as usize]);

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
