use crate::{
    algebra::galois_ring::{teichmuller_subgroup_generator, GrContext, GrElem, GrError, Result},
    protocols::whir_gr::{
        common::WhirGrSumcheckPolynomial, multiquadratic::MultiQuadraticPolynomial,
    },
};

pub type TernaryGrid = [GrElem; 3];

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EqTerm {
    pub weight: GrElem,
    pub point: Vec<GrElem>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WhirConstraint {
    grid: TernaryGrid,
    terms: Vec<EqTerm>,
}

impl WhirConstraint {
    pub const fn new(grid: TernaryGrid) -> Self {
        Self {
            grid,
            terms: Vec::new(),
        }
    }

    pub fn with_terms(grid: TernaryGrid, terms: Vec<EqTerm>) -> Result<Self> {
        let mut constraint = Self::new(grid);
        for term in terms {
            constraint.add_shift_term(term.weight, term.point)?;
        }
        Ok(constraint)
    }

    pub const fn grid(&self) -> &TernaryGrid {
        &self.grid
    }

    pub fn terms(&self) -> &[EqTerm] {
        &self.terms
    }

    pub const fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn variable_count(&self) -> u64 {
        self.terms.first().map_or(0, |term| term.point.len() as u64)
    }

    pub fn add_shift_term(&mut self, weight: GrElem, point: Vec<GrElem>) -> Result<()> {
        if let Some(first) = self.terms.first() {
            if first.point.len() != point.len() {
                return Err(GrError::InvalidPolynomial(
                    "constraint terms must have the same arity",
                ));
            }
        }
        self.terms.push(EqTerm { weight, point });
        Ok(())
    }

    pub fn evaluate_a(&self, ctx: &GrContext, x: &[GrElem]) -> Result<GrElem> {
        require_valid_grid(ctx, &self.grid)?;

        let mut out = ctx.zero();
        for term in &self.terms {
            if term.point.len() != x.len() {
                return Err(GrError::InvalidPolynomial(
                    "constraint point length mismatch",
                ));
            }
            let weighted = ctx.mul(&term.weight, &eq_b(ctx, &self.grid, &term.point, x)?);
            out = ctx.add(&out, &weighted);
        }
        Ok(out)
    }

    pub fn evaluate_w(&self, ctx: &GrContext, z: &GrElem, x: &[GrElem]) -> Result<GrElem> {
        Ok(ctx.mul(z, &self.evaluate_a(ctx, x)?))
    }

    pub fn restrict_prefix(&self, ctx: &GrContext, alphas: &[GrElem]) -> Result<Self> {
        require_valid_grid(ctx, &self.grid)?;

        let mut restricted_terms = Vec::with_capacity(self.terms.len());
        for term in &self.terms {
            if alphas.len() > term.point.len() {
                return Err(GrError::InvalidPolynomial(
                    "constraint prefix fixes too many variables",
                ));
            }

            let mut restricted_weight = term.weight.clone();
            for (point, alpha) in term.point.iter().zip(alphas) {
                restricted_weight = ctx.mul(
                    &restricted_weight,
                    &eq_b_coordinate(ctx, &self.grid, point, alpha)?,
                );
            }
            restricted_terms.push(EqTerm {
                weight: restricted_weight,
                point: term.point[alphas.len()..].to_vec(),
            });
        }

        Self::with_terms(self.grid.clone(), restricted_terms)
    }
}

pub fn ternary_grid(ctx: &GrContext, omega: &GrElem) -> Result<TernaryGrid> {
    if !ctx.is_unit(omega) || *omega == ctx.one() || ctx.pow(omega, 3) != ctx.one() {
        return Err(GrError::InvalidDomain(
            "ternary grid requires omega of exact order 3",
        ));
    }

    let grid = [ctx.one(), omega.clone(), ctx.square(omega)];
    require_valid_grid(ctx, &grid)?;
    Ok(grid)
}

pub fn points_have_pairwise_unit_differences(ctx: &GrContext, points: &[GrElem]) -> bool {
    for i in 0..points.len() {
        for j in i + 1..points.len() {
            if !ctx.is_unit(&ctx.sub(&points[i], &points[j])) {
                return false;
            }
        }
    }
    true
}

pub fn lagrange_basis_on_ternary_grid(
    ctx: &GrContext,
    grid: &TernaryGrid,
    basis_index: usize,
    x: &GrElem,
) -> Result<GrElem> {
    if basis_index >= grid.len() {
        return Err(GrError::IndexOutOfRange {
            index: basis_index as u64,
            size: grid.len() as u64,
        });
    }
    require_valid_grid(ctx, grid)?;

    let basis_point = &grid[basis_index];
    let mut numerator = ctx.one();
    let mut denominator = ctx.one();
    for (index, point) in grid.iter().enumerate() {
        if index == basis_index {
            continue;
        }
        numerator = ctx.mul(&numerator, &ctx.sub(x, point));
        denominator = ctx.mul(&denominator, &ctx.sub(basis_point, point));
    }

    Ok(ctx.mul(&numerator, &ctx.inv(&denominator)?))
}

pub fn eq_b_coordinate(
    ctx: &GrContext,
    grid: &TernaryGrid,
    z: &GrElem,
    x: &GrElem,
) -> Result<GrElem> {
    require_valid_grid(ctx, grid)?;

    let mut out = ctx.zero();
    for index in 0..grid.len() {
        let term = ctx.mul(
            &lagrange_basis_on_ternary_grid(ctx, grid, index, z)?,
            &lagrange_basis_on_ternary_grid(ctx, grid, index, x)?,
        );
        out = ctx.add(&out, &term);
    }
    Ok(out)
}

pub fn eq_b(ctx: &GrContext, grid: &TernaryGrid, z: &[GrElem], x: &[GrElem]) -> Result<GrElem> {
    if z.len() != x.len() {
        return Err(GrError::InvalidPolynomial(
            "eq_b requires equal-length points",
        ));
    }

    let mut out = ctx.one();
    for (z_coordinate, x_coordinate) in z.iter().zip(x) {
        out = ctx.mul(
            &out,
            &eq_b_coordinate(ctx, grid, z_coordinate, x_coordinate)?,
        );
    }
    Ok(out)
}

pub fn sumcheck_interpolation_points(ctx: &GrContext) -> Result<Vec<GrElem>> {
    let root = teichmuller_subgroup_generator(ctx, 7)?;
    let mut points = Vec::with_capacity(5);
    points.push(ctx.zero());
    let mut current = ctx.one();
    for _ in 0..4 {
        points.push(current.clone());
        current = ctx.mul(&current, &root);
    }
    if !points_have_pairwise_unit_differences(ctx, &points) {
        return Err(GrError::InvalidDomain(
            "sumcheck interpolation points require pairwise unit differences",
        ));
    }
    Ok(points)
}

pub fn honest_sumcheck_polynomial(
    ctx: &GrContext,
    polynomial: &MultiQuadraticPolynomial,
    constraint: &WhirConstraint,
    prefix: &[GrElem],
) -> Result<WhirGrSumcheckPolynomial> {
    honest_sumcheck_polynomial_with_evaluator(
        ctx,
        polynomial.variable_count(),
        constraint,
        prefix,
        |point| polynomial.evaluate(ctx, point),
    )
}

pub fn honest_sumcheck_polynomial_with_evaluator<F>(
    ctx: &GrContext,
    variable_count: u64,
    constraint: &WhirConstraint,
    prefix: &[GrElem],
    evaluate_f: F,
) -> Result<WhirGrSumcheckPolynomial>
where
    F: Fn(&[GrElem]) -> Result<GrElem>,
{
    let variable_count_usize = checked_size(variable_count, "variable count")?;
    if prefix.len() >= variable_count_usize {
        return Err(GrError::InvalidPolynomial(
            "honest sumcheck requires one live variable",
        ));
    }
    if !constraint.is_empty() && constraint.variable_count() != variable_count {
        return Err(GrError::InvalidPolynomial(
            "sumcheck constraint arity mismatch",
        ));
    }

    let interpolation_points = sumcheck_interpolation_points(ctx)?;
    let remaining_count = variable_count - prefix.len() as u64 - 1;
    let assignment_count = pow3_small(remaining_count)?;
    let mut full_point = vec![ctx.zero(); variable_count_usize];
    full_point[..prefix.len()].clone_from_slice(prefix);
    let variable_index = prefix.len();
    let suffix_begin = variable_index + 1;

    let mut values = Vec::with_capacity(interpolation_points.len());
    for t in &interpolation_points {
        full_point[variable_index] = t.clone();
        let mut h_t = ctx.zero();
        for assignment in 0..assignment_count {
            append_grid_assignment(
                constraint.grid(),
                assignment,
                &mut full_point[suffix_begin..],
            );
            let product = ctx.mul(
                &evaluate_f(&full_point)?,
                &constraint.evaluate_a(ctx, &full_point)?,
            );
            h_t = ctx.add(&h_t, &product);
        }
        values.push(h_t);
    }

    Ok(WhirGrSumcheckPolynomial {
        coefficients: interpolate(ctx, &interpolation_points, &values)?,
    })
}

pub const fn sumcheck_declared_degree(polynomial: &WhirGrSumcheckPolynomial) -> usize {
    polynomial.coefficients.len().saturating_sub(1)
}

pub fn evaluate_sumcheck_polynomial(
    ctx: &GrContext,
    polynomial: &WhirGrSumcheckPolynomial,
    x: &GrElem,
) -> GrElem {
    let mut out = ctx.zero();
    for coefficient in polynomial.coefficients.iter().rev() {
        out = ctx.add(&ctx.mul(&out, x), coefficient);
    }
    out
}

pub const fn check_sumcheck_degree(
    polynomial: &WhirGrSumcheckPolynomial,
    degree_bound: u64,
) -> bool {
    sumcheck_declared_degree(polynomial) <= degree_bound as usize
}

pub fn sumcheck_grid_sum(
    ctx: &GrContext,
    grid: &TernaryGrid,
    polynomial: &WhirGrSumcheckPolynomial,
) -> Result<GrElem> {
    require_valid_grid(ctx, grid)?;

    let mut out = ctx.zero();
    for point in grid {
        out = ctx.add(&out, &evaluate_sumcheck_polynomial(ctx, polynomial, point));
    }
    Ok(out)
}

pub fn check_sumcheck_identity(
    ctx: &GrContext,
    grid: &TernaryGrid,
    polynomial: &WhirGrSumcheckPolynomial,
    current_sigma: &GrElem,
    degree_bound: u64,
) -> Result<bool> {
    if !check_sumcheck_degree(polynomial, degree_bound) {
        return Ok(false);
    }
    Ok(sumcheck_grid_sum(ctx, grid, polynomial)? == *current_sigma)
}

pub fn sumcheck_next_sigma(
    ctx: &GrContext,
    polynomial: &WhirGrSumcheckPolynomial,
    alpha: &GrElem,
) -> GrElem {
    evaluate_sumcheck_polynomial(ctx, polynomial, alpha)
}

fn require_valid_grid(ctx: &GrContext, grid: &TernaryGrid) -> Result<()> {
    if points_have_pairwise_unit_differences(ctx, grid) {
        Ok(())
    } else {
        Err(GrError::InvalidDomain(
            "ternary grid requires pairwise unit differences",
        ))
    }
}

fn append_grid_assignment(grid: &TernaryGrid, mut assignment: u64, suffix: &mut [GrElem]) {
    for coordinate in suffix {
        *coordinate = grid[(assignment % 3) as usize].clone();
        assignment /= 3;
    }
}

fn interpolate(ctx: &GrContext, points: &[GrElem], values: &[GrElem]) -> Result<Vec<GrElem>> {
    if points.len() != values.len() {
        return Err(GrError::InvalidPolynomial(
            "interpolation point/value length mismatch",
        ));
    }
    if !points_have_pairwise_unit_differences(ctx, points) {
        return Err(GrError::InvalidDomain(
            "interpolation points require pairwise unit differences",
        ));
    }

    let mut out = vec![ctx.zero(); points.len()];
    for (basis_index, point) in points.iter().enumerate() {
        let mut basis = vec![ctx.one()];
        let mut denominator = ctx.one();
        for (other_index, other) in points.iter().enumerate() {
            if basis_index == other_index {
                continue;
            }
            basis = multiply_by_linear_factor(ctx, &basis, &ctx.neg(other));
            denominator = ctx.mul(&denominator, &ctx.sub(point, other));
        }
        let scale = ctx.mul(&values[basis_index], &ctx.inv(&denominator)?);
        for (index, coefficient) in basis.iter().enumerate() {
            out[index] = ctx.add(&out[index], &ctx.mul(coefficient, &scale));
        }
    }
    Ok(trim_trailing_zeros(out))
}

fn multiply_by_linear_factor(
    ctx: &GrContext,
    coefficients: &[GrElem],
    constant: &GrElem,
) -> Vec<GrElem> {
    let mut out = vec![ctx.zero(); coefficients.len() + 1];
    for (index, coefficient) in coefficients.iter().enumerate() {
        out[index] = ctx.add(&out[index], &ctx.mul(coefficient, constant));
        out[index + 1] = ctx.add(&out[index + 1], coefficient);
    }
    out
}

fn trim_trailing_zeros(mut coefficients: Vec<GrElem>) -> Vec<GrElem> {
    while coefficients.last().is_some_and(GrElem::is_zero) {
        coefficients.pop();
    }
    coefficients
}

fn pow3_small(exponent: u64) -> Result<u64> {
    let mut out = 1u64;
    for _ in 0..exponent {
        out = out
            .checked_mul(3)
            .ok_or(GrError::ArithmeticOverflow("constraint pow3"))?;
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
        algebra::galois_ring::{Domain, GrConfig, GrContext, GrElem, GrError},
        protocols::whir_gr::{
            common::WhirGrSumcheckPolynomial,
            constraint::{
                check_sumcheck_degree, check_sumcheck_identity, evaluate_sumcheck_polynomial,
                honest_sumcheck_polynomial, lagrange_basis_on_ternary_grid,
                points_have_pairwise_unit_differences, sumcheck_interpolation_points,
                sumcheck_next_sigma, ternary_grid, EqTerm, TernaryGrid, WhirConstraint,
            },
            multiquadratic::MultiQuadraticPolynomial,
        },
    };

    fn sample_context() -> GrContext {
        GrContext::new(GrConfig {
            p: 2,
            k_exp: 16,
            r: 6,
        })
        .unwrap()
    }

    fn make_grid(ctx: &GrContext) -> TernaryGrid {
        let domain = Domain::teichmuller_subgroup(Arc::new(ctx.clone()), 3).unwrap();
        ternary_grid(ctx, domain.root()).unwrap()
    }

    fn sample_coefficients(ctx: &GrContext, variable_count: u64) -> Vec<GrElem> {
        let count =
            crate::protocols::whir_gr::multiquadratic::pow3_checked(variable_count).unwrap();
        (0..count).map(|i| ctx.from_u64((5 * i + 2) % 13)).collect()
    }

    fn sample_point(ctx: &GrContext, grid: &TernaryGrid, variable_count: u64) -> Vec<GrElem> {
        (0..variable_count)
            .map(|i| ctx.add(&grid[((i + 1) % 3) as usize], &ctx.from_u64(i + 2)))
            .collect()
    }

    fn grid_point_from_index(
        grid: &TernaryGrid,
        variable_count: u64,
        mut index: u64,
    ) -> Vec<GrElem> {
        let mut point = Vec::with_capacity(variable_count as usize);
        for _ in 0..variable_count {
            point.push(grid[(index % 3) as usize].clone());
            index /= 3;
        }
        point
    }

    fn sum_over_grid(
        ctx: &GrContext,
        grid: &TernaryGrid,
        polynomial: &MultiQuadraticPolynomial,
        constraint: &WhirConstraint,
    ) -> GrElem {
        let mut sum = ctx.zero();
        let count =
            crate::protocols::whir_gr::multiquadratic::pow3_checked(polynomial.variable_count())
                .unwrap();
        for index in 0..count {
            let point = grid_point_from_index(grid, polynomial.variable_count(), index);
            let product = ctx.mul(
                &polynomial.evaluate(ctx, &point).unwrap(),
                &constraint.evaluate_a(ctx, &point).unwrap(),
            );
            sum = ctx.add(&sum, &product);
        }
        sum
    }

    #[test]
    fn ternary_grid_and_lagrange_basis_should_match_kronecker_delta() {
        let ctx = sample_context();
        let grid = make_grid(&ctx);

        assert!(points_have_pairwise_unit_differences(&ctx, &grid));
        for i in 0..grid.len() {
            for j in 0..grid.len() {
                let actual = lagrange_basis_on_ternary_grid(&ctx, &grid, i, &grid[j]).unwrap();
                assert_eq!(actual, if i == j { ctx.one() } else { ctx.zero() });
            }
        }
        assert_eq!(sumcheck_interpolation_points(&ctx).unwrap().len(), 5);
    }

    #[test]
    fn equality_kernel_should_reproduce_multiquadratic_on_grid_sum() {
        let ctx = sample_context();
        let grid = make_grid(&ctx);
        let polynomial = MultiQuadraticPolynomial::new(3, sample_coefficients(&ctx, 3)).unwrap();
        let z = sample_point(&ctx, &grid, 3);
        let constraint = WhirConstraint::with_terms(
            grid.clone(),
            vec![EqTerm {
                weight: ctx.one(),
                point: z.clone(),
            }],
        )
        .unwrap();

        assert_eq!(
            sum_over_grid(&ctx, &grid, &polynomial, &constraint),
            polynomial.evaluate(&ctx, &z).unwrap()
        );
    }

    #[test]
    fn constraint_restriction_should_match_direct_evaluation() {
        let ctx = sample_context();
        let grid = make_grid(&ctx);
        let first_point = sample_point(&ctx, &grid, 3);
        let second_point = vec![
            grid[2].clone(),
            ctx.add(&grid[0], &ctx.from_u64(3)),
            grid[1].clone(),
        ];
        let prefix = vec![ctx.add(&grid[1], &ctx.from_u64(5))];
        let tail = vec![ctx.add(&grid[0], &ctx.from_u64(7)), grid[2].clone()];
        let constraint = WhirConstraint::with_terms(
            grid.clone(),
            vec![
                EqTerm {
                    weight: ctx.from_u64(2),
                    point: first_point,
                },
                EqTerm {
                    weight: ctx.add(&grid[1], &ctx.one()),
                    point: second_point,
                },
            ],
        )
        .unwrap();

        let restricted = constraint.restrict_prefix(&ctx, &prefix).unwrap();
        let mut full_point = prefix;
        full_point.extend_from_slice(&tail);

        assert_eq!(
            restricted.evaluate_a(&ctx, &tail).unwrap(),
            constraint.evaluate_a(&ctx, &full_point).unwrap()
        );
    }

    #[test]
    fn honest_sumcheck_should_satisfy_verifier_identity() {
        let ctx = sample_context();
        let grid = make_grid(&ctx);
        let polynomial = MultiQuadraticPolynomial::new(3, sample_coefficients(&ctx, 3)).unwrap();
        let z = sample_point(&ctx, &grid, 3);
        let constraint = WhirConstraint::with_terms(
            grid.clone(),
            vec![EqTerm {
                weight: ctx.one(),
                point: z.clone(),
            }],
        )
        .unwrap();

        let mut current_sigma = polynomial.evaluate(&ctx, &z).unwrap();
        let mut prefix = Vec::new();
        for round in 0..polynomial.variable_count() {
            let h = honest_sumcheck_polynomial(&ctx, &polynomial, &constraint, &prefix).unwrap();
            assert!(check_sumcheck_degree(&h, 4));
            assert!(check_sumcheck_identity(&ctx, &grid, &h, &current_sigma, 4).unwrap());
            let alpha = ctx.add(&grid[(round % 3) as usize], &ctx.from_u64(round + 4));
            current_sigma = sumcheck_next_sigma(&ctx, &h, &alpha);
            prefix.push(alpha);
        }
    }

    #[test]
    fn tampering_and_declared_degree_should_fail() {
        let ctx = sample_context();
        let grid = make_grid(&ctx);
        let polynomial = MultiQuadraticPolynomial::new(2, sample_coefficients(&ctx, 2)).unwrap();
        let z = sample_point(&ctx, &grid, 2);
        let constraint = WhirConstraint::with_terms(
            grid.clone(),
            vec![EqTerm {
                weight: ctx.one(),
                point: z.clone(),
            }],
        )
        .unwrap();
        let honest = honest_sumcheck_polynomial(&ctx, &polynomial, &constraint, &[]).unwrap();
        let current_sigma = polynomial.evaluate(&ctx, &z).unwrap();
        assert!(check_sumcheck_identity(&ctx, &grid, &honest, &current_sigma, 4).unwrap());

        let mut tampered = honest.clone();
        if tampered.coefficients.is_empty() {
            tampered.coefficients.push(ctx.one());
        } else {
            tampered.coefficients[0] = ctx.add(&tampered.coefficients[0], &ctx.one());
        }
        assert!(!check_sumcheck_identity(&ctx, &grid, &tampered, &current_sigma, 4).unwrap());

        let mut declared_too_large = honest;
        declared_too_large.coefficients.resize(6, ctx.zero());
        assert!(!check_sumcheck_degree(&declared_too_large, 4));
    }

    #[test]
    fn invalid_shapes_should_reject() {
        let ctx = sample_context();
        let grid = make_grid(&ctx);
        let mut constraint = WhirConstraint::new(grid.clone());
        constraint
            .add_shift_term(ctx.one(), sample_point(&ctx, &grid, 2))
            .unwrap();

        assert!(matches!(
            constraint.add_shift_term(ctx.one(), sample_point(&ctx, &grid, 3)),
            Err(GrError::InvalidPolynomial(_))
        ));

        let polynomial = MultiQuadraticPolynomial::new(2, sample_coefficients(&ctx, 2)).unwrap();
        let full_prefix = sample_point(&ctx, &grid, 2);
        assert!(matches!(
            honest_sumcheck_polynomial(&ctx, &polynomial, &constraint, &full_prefix),
            Err(GrError::InvalidPolynomial(_))
        ));
    }

    #[test]
    fn sumcheck_polynomial_evaluation_should_use_horner_order() {
        let ctx = sample_context();
        let polynomial = WhirGrSumcheckPolynomial {
            coefficients: vec![ctx.from_u64(2), ctx.from_u64(3), ctx.from_u64(5)],
        };
        let x = ctx.from_u64(7);
        let expected = ctx.add(
            &ctx.from_u64(2),
            &ctx.add(
                &ctx.mul(&ctx.from_u64(3), &x),
                &ctx.mul(&ctx.from_u64(5), &ctx.square(&x)),
            ),
        );

        assert_eq!(
            evaluate_sumcheck_polynomial(&ctx, &polynomial, &x),
            expected
        );
    }
}
