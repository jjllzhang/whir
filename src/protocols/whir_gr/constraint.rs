use std::time::Instant;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::{
    algebra::galois_ring::{
        clear_elem, teichmuller_subgroup_generator, GrContext, GrElem, GrError, GrScratch, Result,
    },
    protocols::whir_gr::{
        common::WhirGrSumcheckPolynomial, multiquadratic::MultiQuadraticPolynomial,
    },
};

pub type TernaryGrid = [GrElem; 3];

const SUMCHECK_ROW_SCRATCH_ELEMS: usize = 8;
const SUMCHECK_ACCUMULATE_SCRATCH_ELEMS: usize = 4;
const PARALLEL_SUMCHECK_ROW_THRESHOLD: usize = 4096;
const PARALLEL_SUMCHECK_ACCUMULATE_THRESHOLD: usize = 4096;
const PARALLEL_TERNARY_TRANSFORM_THRESHOLD: usize = 4096;

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct HonestSumcheckTimings {
    pub constraint_plan: f64,
    pub poly_eval: f64,
    pub accumulate: f64,
    pub interpolate: f64,
}

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
        for term in &self.terms {
            if term.point.len() != x.len() {
                return Err(GrError::InvalidPolynomial(
                    "constraint point length mismatch",
                ));
            }
        }

        let basis_table = TernaryBasisTable::new(ctx, &self.grid)?;
        let x_basis = x
            .iter()
            .map(|coordinate| basis_table.values(ctx, &self.grid, coordinate))
            .collect::<Vec<_>>();

        let mut out = ctx.zero();
        for term in &self.terms {
            let eq_value =
                eq_b_against_precomputed_rhs(ctx, &self.grid, &basis_table, &term.point, &x_basis);
            let weighted = ctx.mul(&term.weight, &eq_value);
            out = ctx.add(&out, &weighted);
        }
        Ok(out)
    }

    pub fn evaluate_w(&self, ctx: &GrContext, z: &GrElem, x: &[GrElem]) -> Result<GrElem> {
        Ok(ctx.mul(z, &self.evaluate_a(ctx, x)?))
    }

    pub fn restrict_prefix(&self, ctx: &GrContext, alphas: &[GrElem]) -> Result<Self> {
        require_valid_grid(ctx, &self.grid)?;
        let basis_table = TernaryBasisTable::new(ctx, &self.grid)?;
        let alpha_basis = alphas
            .iter()
            .map(|alpha| basis_table.values(ctx, &self.grid, alpha))
            .collect::<Vec<_>>();

        let mut restricted_terms = Vec::with_capacity(self.terms.len());
        for term in &self.terms {
            if alphas.len() > term.point.len() {
                return Err(GrError::InvalidPolynomial(
                    "constraint prefix fixes too many variables",
                ));
            }

            let mut restricted_weight = term.weight.clone();
            for (point, alpha_basis) in term.point.iter().zip(&alpha_basis) {
                let point_basis = basis_table.values(ctx, &self.grid, point);
                restricted_weight = ctx.mul(
                    &restricted_weight,
                    &TernaryBasisTable::eq_from_values(ctx, &point_basis, alpha_basis),
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

    let basis_table = TernaryBasisTable::new(ctx, grid)?;
    let z_basis = basis_table.values(ctx, grid, z);
    let x_basis = basis_table.values(ctx, grid, x);
    Ok(TernaryBasisTable::eq_from_values(ctx, &z_basis, &x_basis))
}

pub fn eq_b(ctx: &GrContext, grid: &TernaryGrid, z: &[GrElem], x: &[GrElem]) -> Result<GrElem> {
    if z.len() != x.len() {
        return Err(GrError::InvalidPolynomial(
            "eq_b requires equal-length points",
        ));
    }
    require_valid_grid(ctx, grid)?;

    let basis_table = TernaryBasisTable::new(ctx, grid)?;
    let x_basis = x
        .iter()
        .map(|coordinate| basis_table.values(ctx, grid, coordinate))
        .collect::<Vec<_>>();
    Ok(eq_b_against_precomputed_rhs(
        ctx,
        grid,
        &basis_table,
        z,
        &x_basis,
    ))
}

fn eq_b_against_precomputed_rhs(
    ctx: &GrContext,
    grid: &TernaryGrid,
    basis_table: &TernaryBasisTable,
    z: &[GrElem],
    x_basis: &[[GrElem; 3]],
) -> GrElem {
    let mut out = ctx.one();
    for (z_coordinate, x_basis_coordinate) in z.iter().zip(x_basis) {
        let z_basis = basis_table.values(ctx, grid, z_coordinate);
        out = ctx.mul(
            &out,
            &TernaryBasisTable::eq_from_values(ctx, &z_basis, x_basis_coordinate),
        );
    }
    out
}

struct TernaryBasisTable {
    denominator_inverses: [GrElem; 3],
}

impl TernaryBasisTable {
    fn new(ctx: &GrContext, grid: &TernaryGrid) -> Result<Self> {
        require_valid_grid(ctx, grid)?;
        let mut denominators = Vec::with_capacity(grid.len());
        for (basis_index, basis_point) in grid.iter().enumerate() {
            let mut denominator = ctx.one();
            for (point_index, point) in grid.iter().enumerate() {
                if basis_index != point_index {
                    denominator = ctx.mul(&denominator, &ctx.sub(basis_point, point));
                }
            }
            denominators.push(denominator);
        }
        let denominator_inverses = ctx.batch_inv(&denominators)?.try_into().map_err(|_| {
            GrError::InvalidDomain("ternary basis table requires three grid points")
        })?;
        Ok(Self {
            denominator_inverses,
        })
    }

    fn values(&self, ctx: &GrContext, grid: &TernaryGrid, x: &GrElem) -> [GrElem; 3] {
        let mut out = [ctx.zero(), ctx.zero(), ctx.zero()];
        for (basis_index, out_value) in out.iter_mut().enumerate().take(grid.len()) {
            let mut numerator = ctx.one();
            for (point_index, point) in grid.iter().enumerate() {
                if basis_index != point_index {
                    numerator = ctx.mul(&numerator, &ctx.sub(x, point));
                }
            }
            *out_value = ctx.mul(&numerator, &self.denominator_inverses[basis_index]);
        }
        out
    }

    fn eq_from_values(
        ctx: &GrContext,
        lhs_values: &[GrElem; 3],
        rhs_values: &[GrElem; 3],
    ) -> GrElem {
        let mut out = ctx.zero();
        for (lhs, rhs) in lhs_values.iter().zip(rhs_values) {
            out = ctx.add(&out, &ctx.mul(lhs, rhs));
        }
        out
    }
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
    let restricted;
    let restricted_polynomial = if prefix.is_empty() {
        polynomial
    } else {
        restricted = polynomial.restrict_prefix(ctx, prefix)?;
        &restricted
    };
    Ok(honest_sumcheck_polynomial_for_restricted(
        ctx,
        restricted_polynomial,
        constraint,
        prefix,
        false,
    )?
    .0)
}

pub(crate) fn honest_sumcheck_polynomial_for_restricted(
    ctx: &GrContext,
    polynomial: &MultiQuadraticPolynomial,
    constraint: &WhirConstraint,
    prefix: &[GrElem],
    capture_timings: bool,
) -> Result<(WhirGrSumcheckPolynomial, HonestSumcheckTimings)> {
    let mut timings = HonestSumcheckTimings::default();
    let restricted_variable_count = polynomial.variable_count();
    if restricted_variable_count == 0 {
        return Err(GrError::InvalidPolynomial(
            "honest sumcheck requires one live variable",
        ));
    }
    let full_variable_count = if constraint.is_empty() {
        prefix.len() as u64 + restricted_variable_count
    } else {
        constraint.variable_count()
    };
    if prefix.len() as u64 > full_variable_count
        || restricted_variable_count != full_variable_count - prefix.len() as u64
    {
        return Err(GrError::InvalidPolynomial(
            "restricted sumcheck polynomial arity mismatch",
        ));
    }

    let interpolation_points = sumcheck_interpolation_points(ctx)?;
    let remaining_count = restricted_variable_count - 1;
    let assignment_count = pow3_small(remaining_count)?;
    let assignment_count_usize = checked_size(assignment_count, "sumcheck assignment count")?;

    let constraint_plan_start = capture_timings.then(Instant::now);
    let constraint_plan = SumcheckConstraintPlan::new(
        ctx,
        constraint,
        prefix,
        &interpolation_points,
        remaining_count,
        assignment_count_usize,
    )?;
    record_elapsed(&mut timings.constraint_plan, constraint_plan_start);

    let poly_eval_start = capture_timings.then(Instant::now);
    let polynomial_values = evaluate_restricted_polynomial_rows(
        ctx,
        polynomial,
        constraint.grid(),
        &interpolation_points,
    )?;
    record_elapsed(&mut timings.poly_eval, poly_eval_start);

    let accumulate_start = capture_timings.then(Instant::now);
    let values = accumulate_sumcheck_values(
        ctx,
        &polynomial_values,
        &constraint_plan,
        interpolation_points.len(),
        assignment_count_usize,
    );
    record_elapsed(&mut timings.accumulate, accumulate_start);

    let interpolate_start = capture_timings.then(Instant::now);
    let coefficients = interpolate(ctx, &interpolation_points, &values)?;
    record_elapsed(&mut timings.interpolate, interpolate_start);

    Ok((WhirGrSumcheckPolynomial { coefficients }, timings))
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
    let assignment_count_usize = checked_size(assignment_count, "sumcheck assignment count")?;
    let constraint_plan = SumcheckConstraintPlan::new(
        ctx,
        constraint,
        prefix,
        &interpolation_points,
        remaining_count,
        assignment_count_usize,
    )?;
    let mut full_point = vec![ctx.zero(); variable_count_usize];
    full_point[..prefix.len()].clone_from_slice(prefix);
    let variable_index = prefix.len();
    let suffix_begin = variable_index + 1;

    let mut values = Vec::with_capacity(interpolation_points.len());
    for (t_index, t) in interpolation_points.iter().enumerate() {
        full_point[variable_index] = t.clone();
        let mut h_t = ctx.zero();
        for assignment in 0..assignment_count_usize {
            append_grid_assignment(
                constraint.grid(),
                assignment as u64,
                &mut full_point[suffix_begin..],
            );
            let product = ctx.mul(
                &evaluate_f(&full_point)?,
                constraint_plan.value(t_index, assignment),
            );
            h_t = ctx.add(&h_t, &product);
        }
        values.push(h_t);
    }

    Ok(WhirGrSumcheckPolynomial {
        coefficients: interpolate(ctx, &interpolation_points, &values)?,
    })
}

fn evaluate_restricted_polynomial_rows(
    ctx: &GrContext,
    polynomial: &MultiQuadraticPolynomial,
    grid: &TernaryGrid,
    interpolation_points: &[GrElem],
) -> Result<Vec<GrElem>> {
    require_valid_grid(ctx, grid)?;
    let variable_count = checked_size(polynomial.variable_count(), "restricted variable count")?;
    if variable_count == 0 {
        return Err(GrError::InvalidPolynomial(
            "restricted sumcheck polynomial requires one live variable",
        ));
    }

    let coefficient_count = checked_size(
        pow3_small(polynomial.variable_count())?,
        "restricted coefficient count",
    )?;
    let mut work = vec![ctx.zero(); coefficient_count];
    for (target, coefficient) in work.iter_mut().zip(polynomial.coefficients()) {
        target.clone_from(coefficient);
    }

    for dimension in 1..variable_count {
        apply_ternary_grid_transform(ctx, &mut work, dimension, grid)?;
    }

    let assignment_count = checked_size(
        pow3_small(polynomial.variable_count() - 1)?,
        "sumcheck assignment count",
    )?;
    #[cfg(feature = "parallel")]
    {
        let row_count = interpolation_points.len() * assignment_count;
        if should_parallel_sumcheck_rows(row_count) {
            let row_blocks = interpolation_points
                .par_iter()
                .map(|t| evaluate_restricted_polynomial_row(ctx, &work, t, assignment_count))
                .collect::<Vec<_>>();
            let mut rows = Vec::with_capacity(row_count);
            for row in row_blocks {
                rows.extend(row);
            }
            return Ok(rows);
        }
    }

    let mut rows = Vec::with_capacity(interpolation_points.len() * assignment_count);
    for t in interpolation_points {
        rows.extend(evaluate_restricted_polynomial_row(
            ctx,
            &work,
            t,
            assignment_count,
        ));
    }
    Ok(rows)
}

fn evaluate_restricted_polynomial_row(
    ctx: &GrContext,
    work: &[GrElem],
    t: &GrElem,
    assignment_count: usize,
) -> Vec<GrElem> {
    let mut scratch = GrScratch::with_elements(ctx, SUMCHECK_ROW_SCRATCH_ELEMS);
    let mut t_squared = ctx.zero();
    {
        let (_, mul_scratch) = scratch.parts_mut();
        ctx.square_into(&mut t_squared, t, mul_scratch);
    }

    let mut row = Vec::with_capacity(assignment_count);
    let (elements, mul_scratch) = scratch.parts_mut();
    let [linear, quadratic, linear_plus_quadratic, value, ..] = elements else {
        unreachable!("sumcheck row scratch has enough elements");
    };
    for assignment in 0..assignment_count {
        let base = 3 * assignment;
        ctx.mul_into(linear, &work[base + 1], t, mul_scratch);
        ctx.mul_into(quadratic, &work[base + 2], &t_squared, mul_scratch);
        ctx.add_into(linear_plus_quadratic, linear, quadratic);
        ctx.add_into(value, &work[base], linear_plus_quadratic);
        row.push(value.clone());
    }
    row
}

fn apply_ternary_grid_transform(
    ctx: &GrContext,
    values: &mut [GrElem],
    dimension: usize,
    grid: &TernaryGrid,
) -> Result<()> {
    let stride = checked_size(pow3_small(dimension as u64)?, "transform stride")?;
    let block = stride
        .checked_mul(3)
        .ok_or(GrError::ArithmeticOverflow("transform block"))?;
    if block == 0 || !values.len().is_multiple_of(block) {
        return Err(GrError::InvalidPolynomial(
            "ternary transform input length mismatch",
        ));
    }

    let grid_squares = [
        ctx.square(&grid[0]),
        ctx.square(&grid[1]),
        ctx.square(&grid[2]),
    ];

    #[cfg(feature = "parallel")]
    {
        if should_parallel_ternary_transform(values.len()) {
            values.par_chunks_mut(block).for_each(|block_values| {
                let mut scratch = GrScratch::with_elements(ctx, SUMCHECK_ROW_SCRATCH_ELEMS);
                for offset in 0..stride {
                    apply_ternary_grid_transform_slot(
                        ctx,
                        block_values,
                        offset,
                        stride,
                        grid,
                        &grid_squares,
                        &mut scratch,
                    );
                }
            });
            return Ok(());
        }
    }

    let mut scratch = GrScratch::with_elements(ctx, SUMCHECK_ROW_SCRATCH_ELEMS);
    for block_begin in (0..values.len()).step_by(block) {
        for offset in 0..stride {
            apply_ternary_grid_transform_slot(
                ctx,
                &mut values[block_begin..block_begin + block],
                offset,
                stride,
                grid,
                &grid_squares,
                &mut scratch,
            );
        }
    }
    Ok(())
}

fn apply_ternary_grid_transform_slot(
    ctx: &GrContext,
    block_values: &mut [GrElem],
    offset: usize,
    stride: usize,
    grid: &TernaryGrid,
    grid_squares: &[GrElem; 3],
    scratch: &mut GrScratch,
) {
    let c0 = block_values[offset].clone();
    let c1 = block_values[stride + offset].clone();
    let c2 = block_values[2 * stride + offset].clone();
    let (elements, mul_scratch) = scratch.parts_mut();
    let [linear, quadratic, linear_plus_quadratic, value, ..] = elements else {
        unreachable!("ternary transform scratch has enough elements");
    };
    for (grid_index, point) in grid.iter().enumerate() {
        ctx.mul_into(linear, &c1, point, mul_scratch);
        ctx.mul_into(quadratic, &c2, &grid_squares[grid_index], mul_scratch);
        ctx.add_into(linear_plus_quadratic, linear, quadratic);
        ctx.add_into(value, &c0, linear_plus_quadratic);
        block_values[grid_index * stride + offset].clone_from(value);
    }
}

struct SumcheckConstraintPlan {
    assignment_count: usize,
    values: Vec<GrElem>,
}

impl SumcheckConstraintPlan {
    fn new(
        ctx: &GrContext,
        constraint: &WhirConstraint,
        prefix: &[GrElem],
        interpolation_points: &[GrElem],
        remaining_count: u64,
        assignment_count: usize,
    ) -> Result<Self> {
        require_valid_grid(ctx, constraint.grid())?;
        let basis_table = TernaryBasisTable::new(ctx, constraint.grid())?;
        let interpolation_basis = interpolation_points
            .iter()
            .map(|point| basis_table.values(ctx, constraint.grid(), point))
            .collect::<Vec<_>>();
        let mut values = vec![ctx.zero(); interpolation_points.len() * assignment_count];
        if constraint.is_empty() {
            return Ok(Self {
                assignment_count,
                values,
            });
        }

        let live_index = prefix.len();
        let suffix_begin = live_index + 1;
        for term in constraint.terms() {
            if term.point.len() < suffix_begin {
                return Err(GrError::InvalidPolynomial(
                    "sumcheck constraint term arity mismatch",
                ));
            }

            let mut fixed_weight = term.weight.clone();
            for (point, alpha) in term.point.iter().zip(prefix) {
                let point_basis = basis_table.values(ctx, constraint.grid(), point);
                let alpha_basis = basis_table.values(ctx, constraint.grid(), alpha);
                fixed_weight = ctx.mul(
                    &fixed_weight,
                    &TernaryBasisTable::eq_from_values(ctx, &point_basis, &alpha_basis),
                );
            }

            let live_point_basis =
                basis_table.values(ctx, constraint.grid(), &term.point[live_index]);
            let live_weights = interpolation_points
                .iter()
                .enumerate()
                .map(|(t_index, _)| {
                    ctx.mul(
                        &fixed_weight,
                        &TernaryBasisTable::eq_from_values(
                            ctx,
                            &live_point_basis,
                            &interpolation_basis[t_index],
                        ),
                    )
                })
                .collect::<Vec<_>>();
            let suffix_basis = term.point[suffix_begin..]
                .iter()
                .map(|point| basis_table.values(ctx, constraint.grid(), point))
                .collect::<Vec<_>>();
            if suffix_basis.len() != checked_size(remaining_count, "sumcheck suffix count")? {
                return Err(GrError::InvalidPolynomial(
                    "sumcheck suffix basis length mismatch",
                ));
            }

            let suffix_products =
                suffix_products(ctx, &suffix_basis, assignment_count, remaining_count)?;
            for (t_index, live_weight) in live_weights.iter().enumerate() {
                let row_offset = t_index * assignment_count;
                add_live_weighted_suffix_products(
                    ctx,
                    &mut values[row_offset..row_offset + assignment_count],
                    live_weight,
                    &suffix_products,
                );
            }
        }

        Ok(Self {
            assignment_count,
            values,
        })
    }

    fn value(&self, t_index: usize, assignment: usize) -> &GrElem {
        &self.values[t_index * self.assignment_count + assignment]
    }
}

fn add_live_weighted_suffix_products(
    ctx: &GrContext,
    row_values: &mut [GrElem],
    live_weight: &GrElem,
    suffix_products: &[GrElem],
) {
    debug_assert_eq!(row_values.len(), suffix_products.len());
    #[cfg(feature = "parallel")]
    {
        if should_parallel_sumcheck_accumulate(row_values.len()) {
            let target_chunks = rayon::current_num_threads().saturating_mul(4).max(1);
            let chunk_size = row_values.len().div_ceil(target_chunks).max(1);
            row_values
                .par_chunks_mut(chunk_size)
                .zip(suffix_products.par_chunks(chunk_size))
                .for_each(|(row_chunk, suffix_chunk)| {
                    let mut scratch =
                        GrScratch::with_elements(ctx, SUMCHECK_ACCUMULATE_SCRATCH_ELEMS);
                    let (elements, mul_scratch) = scratch.parts_mut();
                    let [contribution, ..] = elements else {
                        unreachable!("sumcheck contribution scratch has at least one element");
                    };
                    for (slot, suffix_product) in row_chunk.iter_mut().zip(suffix_chunk) {
                        ctx.mul_into(contribution, live_weight, suffix_product, mul_scratch);
                        ctx.add_assign(slot, contribution);
                    }
                });
            return;
        }
    }

    let mut scratch = GrScratch::with_elements(ctx, SUMCHECK_ACCUMULATE_SCRATCH_ELEMS);
    let (elements, mul_scratch) = scratch.parts_mut();
    let [contribution, ..] = elements else {
        unreachable!("sumcheck contribution scratch has at least one element");
    };
    for (slot, suffix_product) in row_values.iter_mut().zip(suffix_products) {
        ctx.mul_into(contribution, live_weight, suffix_product, mul_scratch);
        ctx.add_assign(slot, contribution);
    }
}

fn accumulate_sumcheck_values(
    ctx: &GrContext,
    polynomial_values: &[GrElem],
    constraint_plan: &SumcheckConstraintPlan,
    interpolation_count: usize,
    assignment_count: usize,
) -> Vec<GrElem> {
    (0..interpolation_count)
        .map(|t_index| {
            accumulate_sumcheck_value(
                ctx,
                polynomial_values,
                constraint_plan,
                t_index,
                assignment_count,
            )
        })
        .collect()
}

fn accumulate_sumcheck_value(
    ctx: &GrContext,
    polynomial_values: &[GrElem],
    constraint_plan: &SumcheckConstraintPlan,
    t_index: usize,
    assignment_count: usize,
) -> GrElem {
    #[cfg(feature = "parallel")]
    {
        if should_parallel_sumcheck_accumulate(assignment_count) {
            let target_chunks = rayon::current_num_threads().saturating_mul(4).max(1);
            let chunk_size = assignment_count.div_ceil(target_chunks).max(1);
            let starts = (0..assignment_count)
                .step_by(chunk_size)
                .collect::<Vec<_>>();
            let partials = starts
                .par_iter()
                .map(|&begin| {
                    let end = (begin + chunk_size).min(assignment_count);
                    accumulate_sumcheck_chunk(
                        ctx,
                        polynomial_values,
                        constraint_plan,
                        t_index,
                        begin,
                        end,
                    )
                })
                .collect::<Vec<_>>();
            return add_sumcheck_partials(ctx, &partials);
        }
    }

    accumulate_sumcheck_chunk(
        ctx,
        polynomial_values,
        constraint_plan,
        t_index,
        0,
        assignment_count,
    )
}

fn accumulate_sumcheck_chunk(
    ctx: &GrContext,
    polynomial_values: &[GrElem],
    constraint_plan: &SumcheckConstraintPlan,
    t_index: usize,
    begin: usize,
    end: usize,
) -> GrElem {
    let row_offset = t_index * constraint_plan.assignment_count;
    let mut scratch = GrScratch::with_elements(ctx, SUMCHECK_ACCUMULATE_SCRATCH_ELEMS);
    let (elements, mul_scratch) = scratch.parts_mut();
    let [product, sum, sum_next, ..] = elements else {
        unreachable!("sumcheck accumulation scratch has enough elements");
    };
    clear_elem(sum);
    for assignment in begin..end {
        ctx.mul_into(
            product,
            &polynomial_values[row_offset + assignment],
            constraint_plan.value(t_index, assignment),
            mul_scratch,
        );
        ctx.add_into(sum_next, sum, product);
        std::mem::swap(sum, sum_next);
    }
    sum.clone()
}

fn add_sumcheck_partials(ctx: &GrContext, partials: &[GrElem]) -> GrElem {
    let mut scratch = GrScratch::with_elements(ctx, 2);
    let (elements, _) = scratch.parts_mut();
    let [sum, sum_next] = elements else {
        unreachable!("partial sum scratch has exactly two elements");
    };
    clear_elem(sum);
    for partial in partials {
        ctx.add_into(sum_next, sum, partial);
        std::mem::swap(sum, sum_next);
    }
    sum.clone()
}

fn suffix_products(
    ctx: &GrContext,
    suffix_basis: &[[GrElem; 3]],
    assignment_count: usize,
    remaining_count: u64,
) -> Result<Vec<GrElem>> {
    let suffix_count = checked_size(remaining_count, "suffix count")?;
    let mut products = vec![ctx.one()];
    for basis in suffix_basis.iter().take(suffix_count) {
        products = expand_suffix_products(ctx, &products, basis);
    }
    if products.len() != assignment_count {
        return Err(GrError::InvalidPolynomial("suffix product count mismatch"));
    }
    Ok(products)
}

fn expand_suffix_products(
    ctx: &GrContext,
    products: &[GrElem],
    basis: &[GrElem; 3],
) -> Vec<GrElem> {
    let mut out = Vec::with_capacity(products.len() * 3);
    for basis_value in basis {
        out.extend(multiply_suffix_products(ctx, products, basis_value));
    }
    out
}

fn multiply_suffix_products(
    ctx: &GrContext,
    products: &[GrElem],
    basis_value: &GrElem,
) -> Vec<GrElem> {
    #[cfg(feature = "parallel")]
    {
        if should_parallel_sumcheck_accumulate(products.len()) {
            let target_chunks = rayon::current_num_threads().saturating_mul(4).max(1);
            let chunk_size = products.len().div_ceil(target_chunks).max(1);
            let chunks = products
                .par_chunks(chunk_size)
                .map(|chunk| multiply_suffix_product_chunk(ctx, chunk, basis_value))
                .collect::<Vec<_>>();
            let mut out = Vec::with_capacity(products.len());
            for chunk in chunks {
                out.extend(chunk);
            }
            return out;
        }
    }
    multiply_suffix_product_chunk(ctx, products, basis_value)
}

fn multiply_suffix_product_chunk(
    ctx: &GrContext,
    products: &[GrElem],
    basis_value: &GrElem,
) -> Vec<GrElem> {
    let mut scratch = GrScratch::with_elements(ctx, 1);
    let (elements, mul_scratch) = scratch.parts_mut();
    let [product, ..] = elements else {
        unreachable!("suffix product scratch has at least one element");
    };
    let mut out = Vec::with_capacity(products.len());
    for value in products {
        ctx.mul_into(product, value, basis_value, mul_scratch);
        out.push(product.clone());
    }
    out
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

#[cfg(feature = "parallel")]
fn should_parallel_sumcheck_rows(row_count: usize) -> bool {
    row_count >= PARALLEL_SUMCHECK_ROW_THRESHOLD && rayon::current_num_threads() > 1
}

#[cfg(feature = "parallel")]
fn should_parallel_sumcheck_accumulate(assignment_count: usize) -> bool {
    assignment_count >= PARALLEL_SUMCHECK_ACCUMULATE_THRESHOLD && rayon::current_num_threads() > 1
}

#[cfg(feature = "parallel")]
fn should_parallel_ternary_transform(value_count: usize) -> bool {
    value_count >= PARALLEL_TERNARY_TRANSFORM_THRESHOLD && rayon::current_num_threads() > 1
}

fn record_elapsed(slot: &mut f64, start: Option<Instant>) {
    if let Some(start) = start {
        *slot += start.elapsed().as_secs_f64() * 1000.0;
    }
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
                honest_sumcheck_polynomial, honest_sumcheck_polynomial_with_evaluator,
                lagrange_basis_on_ternary_grid, points_have_pairwise_unit_differences,
                sumcheck_interpolation_points, sumcheck_next_sigma, ternary_grid, EqTerm,
                TernaryGrid, WhirConstraint,
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
    fn sumcheck_constraint_plan_should_match_direct_evaluation() {
        let ctx = sample_context();
        let grid = make_grid(&ctx);
        let first_point = sample_point(&ctx, &grid, 3);
        let second_point = vec![
            grid[2].clone(),
            ctx.add(&grid[0], &ctx.from_u64(3)),
            grid[1].clone(),
        ];
        let prefix = vec![ctx.add(&grid[1], &ctx.from_u64(5))];
        let interpolation_points = sumcheck_interpolation_points(&ctx).unwrap();
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
        let plan = super::SumcheckConstraintPlan::new(
            &ctx,
            &constraint,
            &prefix,
            &interpolation_points,
            1,
            3,
        )
        .unwrap();

        for (t_index, t) in interpolation_points.iter().enumerate() {
            for (assignment, grid_point) in grid.iter().enumerate() {
                let mut point = prefix.clone();
                point.push(t.clone());
                point.push(grid_point.clone());

                assert_eq!(
                    plan.value(t_index, assignment),
                    &constraint.evaluate_a(&ctx, &point).unwrap()
                );
            }
        }
    }

    #[test]
    fn restricted_polynomial_rows_should_match_direct_evaluation() {
        let ctx = sample_context();
        let grid = make_grid(&ctx);
        let interpolation_points = sumcheck_interpolation_points(&ctx).unwrap();

        for variable_count in 1..=4 {
            let polynomial = MultiQuadraticPolynomial::new(
                variable_count,
                sample_coefficients(&ctx, variable_count),
            )
            .unwrap();
            let rows = super::evaluate_restricted_polynomial_rows(
                &ctx,
                &polynomial,
                &grid,
                &interpolation_points,
            )
            .unwrap();
            let assignment_count = crate::protocols::whir_gr::multiquadratic::pow3_checked(
                variable_count - 1,
            )
            .unwrap() as usize;

            for (t_index, t) in interpolation_points.iter().enumerate() {
                for assignment in 0..assignment_count {
                    let mut point = vec![t.clone()];
                    point.extend(grid_point_from_index(
                        &grid,
                        variable_count - 1,
                        assignment as u64,
                    ));
                    assert_eq!(
                        rows[t_index * assignment_count + assignment],
                        polynomial.evaluate(&ctx, &point).unwrap()
                    );
                }
            }
        }
    }

    #[test]
    fn batched_sumcheck_should_match_pointwise_evaluator() {
        let ctx = sample_context();
        let grid = make_grid(&ctx);
        let polynomial = MultiQuadraticPolynomial::new(4, sample_coefficients(&ctx, 4)).unwrap();
        let first_point = sample_point(&ctx, &grid, 4);
        let second_point = vec![
            grid[2].clone(),
            ctx.add(&grid[0], &ctx.from_u64(3)),
            grid[1].clone(),
            ctx.add(&grid[2], &ctx.from_u64(7)),
        ];
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

        for prefix_len in 0..polynomial.variable_count() as usize {
            let prefix = sample_point(&ctx, &grid, prefix_len as u64);
            let batched =
                honest_sumcheck_polynomial(&ctx, &polynomial, &constraint, &prefix).unwrap();
            let pointwise = honest_sumcheck_polynomial_with_evaluator(
                &ctx,
                polynomial.variable_count(),
                &constraint,
                &prefix,
                |point| polynomial.evaluate(&ctx, point),
            )
            .unwrap();

            assert_eq!(batched, pointwise);
        }
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
