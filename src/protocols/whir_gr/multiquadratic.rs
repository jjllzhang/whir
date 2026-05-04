use crate::algebra::galois_ring::{GrContext, GrElem, GrError, Result};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultiQuadraticPolynomial {
    variable_count: u64,
    coefficients: Vec<GrElem>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultilinearPolynomial {
    variable_count: u64,
    coefficients: Vec<GrElem>,
}

impl MultiQuadraticPolynomial {
    pub fn new(variable_count: u64, coefficients: Vec<GrElem>) -> Result<Self> {
        let max_coefficients = checked_size(pow3_checked(variable_count)?, "coefficient bound")?;
        if coefficients.len() > max_coefficients {
            return Err(GrError::InvalidPolynomial(
                "multi-quadratic coefficient length exceeds 3^m",
            ));
        }

        Ok(Self {
            variable_count,
            coefficients: trim_trailing_zeros(coefficients),
        })
    }

    pub const fn variable_count(&self) -> u64 {
        self.variable_count
    }

    pub fn coefficients(&self) -> &[GrElem] {
        &self.coefficients
    }

    pub fn evaluate(&self, ctx: &GrContext, point: &[GrElem]) -> Result<GrElem> {
        if point.len() != checked_size(self.variable_count, "variable count")? {
            return Err(GrError::InvalidPolynomial(
                "multi-quadratic point length mismatch",
            ));
        }

        let mut acc = ctx.zero();
        for (index, coefficient) in self.coefficients.iter().enumerate() {
            if coefficient.is_zero() {
                continue;
            }

            let mut digits = index as u64;
            let mut monomial = ctx.one();
            for coordinate in point.iter().take(self.variable_count as usize) {
                match digits % 3 {
                    0 => {}
                    1 => {
                        monomial = ctx.mul(&monomial, coordinate);
                    }
                    2 => {
                        monomial = ctx.mul(&monomial, &ctx.square(coordinate));
                    }
                    _ => unreachable!("digit is reduced modulo 3"),
                }
                digits /= 3;
            }
            acc = ctx.add(&acc, &ctx.mul(coefficient, &monomial));
        }
        Ok(acc)
    }

    pub fn evaluate_pow(&self, ctx: &GrContext, x: &GrElem) -> GrElem {
        let mut acc = ctx.zero();
        for coefficient in self.coefficients.iter().rev() {
            acc = ctx.add(&ctx.mul(&acc, x), coefficient);
        }
        acc
    }

    pub fn restrict_prefix(&self, ctx: &GrContext, alphas: &[GrElem]) -> Result<Self> {
        if alphas.len() > checked_size(self.variable_count, "variable count")? {
            return Err(GrError::InvalidPolynomial(
                "multi-quadratic prefix fixes too many variables",
            ));
        }

        let fixed_count = alphas.len() as u64;
        let remaining_count = self.variable_count - fixed_count;
        let prefix_size = pow3_checked(fixed_count)?;
        let tail_bound = pow3_checked(remaining_count)?;

        if self.coefficients.is_empty() {
            return Self::new(remaining_count, Vec::new());
        }

        let required_tail_terms = (self.coefficients.len() as u64).div_ceil(prefix_size);
        let output_terms = required_tail_terms.min(tail_bound);
        let mut restricted = vec![ctx.zero(); checked_size(output_terms, "restricted terms")?];

        for (index, coefficient) in self.coefficients.iter().enumerate() {
            let flat_index = index as u64;
            let prefix_index = flat_index % prefix_size;
            let tail_index = checked_size(flat_index / prefix_size, "tail index")?;
            let weighted = ctx.mul(coefficient, &prefix_weight(ctx, prefix_index, alphas));
            restricted[tail_index] = ctx.add(&restricted[tail_index], &weighted);
        }

        Self::new(remaining_count, restricted)
    }

    pub fn to_univariate_pow_coefficients(&self) -> &[GrElem] {
        &self.coefficients
    }
}

impl MultilinearPolynomial {
    pub fn new(variable_count: u64, coefficients: Vec<GrElem>) -> Result<Self> {
        let max_coefficients = checked_size(
            pow2_checked(variable_count)?,
            "multilinear coefficient bound",
        )?;
        if coefficients.len() > max_coefficients {
            return Err(GrError::InvalidPolynomial(
                "multilinear coefficient length exceeds 2^m",
            ));
        }

        Ok(Self {
            variable_count,
            coefficients: trim_trailing_zeros(coefficients),
        })
    }

    pub const fn variable_count(&self) -> u64 {
        self.variable_count
    }

    pub fn coefficients(&self) -> &[GrElem] {
        &self.coefficients
    }

    pub fn evaluate(&self, ctx: &GrContext, point: &[GrElem]) -> Result<GrElem> {
        if point.len() != checked_size(self.variable_count, "variable count")? {
            return Err(GrError::InvalidPolynomial(
                "multilinear point length mismatch",
            ));
        }

        let mut acc = ctx.zero();
        for (index, coefficient) in self.coefficients.iter().enumerate() {
            if coefficient.is_zero() {
                continue;
            }

            let mut bits = index as u64;
            let mut monomial = ctx.one();
            for coordinate in point.iter().take(self.variable_count as usize) {
                if bits & 1 == 1 {
                    monomial = ctx.mul(&monomial, coordinate);
                }
                bits >>= 1;
            }
            acc = ctx.add(&acc, &ctx.mul(coefficient, &monomial));
        }
        Ok(acc)
    }

    pub fn to_multi_quadratic(&self, ctx: &GrContext) -> Result<MultiQuadraticPolynomial> {
        let embedded_len = checked_size(
            pow3_checked(self.variable_count)?,
            "embedded multilinear coefficient count",
        )?;
        let mut embedded = vec![ctx.zero(); embedded_len];
        for (index, coefficient) in self.coefficients.iter().enumerate() {
            let ternary_index = checked_size(
                binary_index_to_ternary_index(index as u64, self.variable_count)?,
                "embedded multilinear coefficient index",
            )?;
            embedded[ternary_index] = coefficient.clone();
        }
        MultiQuadraticPolynomial::new(self.variable_count, embedded)
    }

    pub fn evaluate_pow(&self, ctx: &GrContext, x: &GrElem) -> Result<GrElem> {
        Ok(self.to_multi_quadratic(ctx)?.evaluate_pow(ctx, x))
    }
}

pub fn pow3_checked(exponent: u64) -> Result<u64> {
    let mut out = 1u64;
    for _ in 0..exponent {
        out = out
            .checked_mul(3)
            .ok_or(GrError::ArithmeticOverflow("pow3_checked"))?;
    }
    Ok(out)
}

pub fn pow2_checked(exponent: u64) -> Result<u64> {
    let mut out = 1u64;
    for _ in 0..exponent {
        out = out
            .checked_mul(2)
            .ok_or(GrError::ArithmeticOverflow("pow2_checked"))?;
    }
    Ok(out)
}

pub fn encode_base3_index(digits: &[u8]) -> Result<u64> {
    let mut index = 0u64;
    let mut place = 1u64;
    for (position, &digit) in digits.iter().enumerate() {
        if digit > 2 {
            return Err(GrError::InvalidPolynomial(
                "base-3 digits must be in [0, 2]",
            ));
        }
        let term = place
            .checked_mul(u64::from(digit))
            .ok_or(GrError::ArithmeticOverflow("encode_base3_index"))?;
        index = index
            .checked_add(term)
            .ok_or(GrError::ArithmeticOverflow("encode_base3_index"))?;
        if position + 1 < digits.len() {
            place = place
                .checked_mul(3)
                .ok_or(GrError::ArithmeticOverflow("encode_base3_index"))?;
        }
    }
    Ok(index)
}

pub fn decode_base3_index(mut index: u64, digit_count: u64) -> Result<Vec<u8>> {
    let bound = pow3_checked(digit_count)?;
    if index >= bound {
        return Err(GrError::IndexOutOfRange { index, size: bound });
    }

    let mut digits = Vec::with_capacity(checked_size(digit_count, "digit count")?);
    for _ in 0..digit_count {
        digits.push((index % 3) as u8);
        index /= 3;
    }
    Ok(digits)
}

pub fn pow_m(ctx: &GrContext, x: &GrElem, variable_count: u64) -> Result<Vec<GrElem>> {
    let mut powers = Vec::with_capacity(checked_size(variable_count, "variable count")?);
    let mut current = x.clone();
    for _ in 0..variable_count {
        powers.push(current.clone());
        current = ctx.mul(&ctx.square(&current), &current);
    }
    Ok(powers)
}

fn prefix_weight(ctx: &GrContext, mut prefix_index: u64, alphas: &[GrElem]) -> GrElem {
    let mut weight = ctx.one();
    for alpha in alphas {
        match prefix_index % 3 {
            0 => {}
            1 => {
                weight = ctx.mul(&weight, alpha);
            }
            2 => {
                weight = ctx.mul(&weight, &ctx.square(alpha));
            }
            _ => unreachable!("digit is reduced modulo 3"),
        }
        prefix_index /= 3;
    }
    weight
}

fn binary_index_to_ternary_index(mut binary_index: u64, variable_count: u64) -> Result<u64> {
    let mut ternary_index = 0u64;
    let mut ternary_place = 1u64;
    for variable in 0..variable_count {
        if binary_index & 1 == 1 {
            ternary_index = ternary_index
                .checked_add(ternary_place)
                .ok_or(GrError::ArithmeticOverflow("binary_index_to_ternary_index"))?;
        }
        binary_index >>= 1;
        if variable + 1 < variable_count {
            ternary_place = ternary_place
                .checked_mul(3)
                .ok_or(GrError::ArithmeticOverflow("binary_index_to_ternary_index"))?;
        }
    }
    Ok(ternary_index)
}

fn trim_trailing_zeros(mut coefficients: Vec<GrElem>) -> Vec<GrElem> {
    while coefficients.last().is_some_and(GrElem::is_zero) {
        coefficients.pop();
    }
    coefficients
}

fn checked_size(value: u64, label: &'static str) -> Result<usize> {
    usize::try_from(value).map_err(|_| GrError::ArithmeticOverflow(label))
}

#[cfg(test)]
mod tests {
    use crate::{
        algebra::galois_ring::{Domain, GrConfig, GrContext, GrError},
        protocols::whir_gr::multiquadratic::{
            decode_base3_index, encode_base3_index, pow2_checked, pow3_checked, pow_m,
            MultiQuadraticPolynomial, MultilinearPolynomial,
        },
    };

    fn sample_context() -> GrContext {
        GrContext::new(GrConfig {
            p: 2,
            k_exp: 16,
            r: 2,
        })
        .unwrap()
    }

    fn sample_coefficients(
        ctx: &GrContext,
        count: u64,
    ) -> Vec<crate::algebra::galois_ring::GrElem> {
        (0..count).map(|i| ctx.from_u64((i * 7 + 3) % 11)).collect()
    }

    #[test]
    fn base3_helpers_should_use_little_endian_digits() {
        let digits = vec![2, 0, 1, 2];

        let index = encode_base3_index(&digits).unwrap();
        let decoded = decode_base3_index(index, 4).unwrap();

        assert_eq!(pow3_checked(5).unwrap(), 243);
        assert_eq!(index, 65);
        assert_eq!(decoded, digits);
    }

    #[test]
    fn pow_m_should_return_repeated_cubes() {
        let ctx = sample_context();
        let x = ctx.from_u64(5);

        let powers = pow_m(&ctx, &x, 4).unwrap();

        assert_eq!(powers.len(), 4);
        assert_eq!(powers[0], x);
        assert_eq!(powers[1], ctx.pow(&x, 3));
        assert_eq!(powers[2], ctx.pow(&x, 9));
        assert_eq!(powers[3], ctx.pow(&x, 27));
    }

    #[test]
    fn evaluate_pow_should_match_evaluation_at_pow_m() {
        let ctx = sample_context();
        let poly = MultiQuadraticPolynomial::new(3, sample_coefficients(&ctx, 27)).unwrap();
        let x = ctx.from_u64(7);
        let pow_point = pow_m(&ctx, &x, poly.variable_count()).unwrap();

        assert_eq!(
            poly.evaluate(&ctx, &pow_point).unwrap(),
            poly.evaluate_pow(&ctx, &x)
        );
    }

    #[test]
    fn restrict_prefix_should_match_original_evaluation() {
        let ctx = sample_context();
        let poly = MultiQuadraticPolynomial::new(3, sample_coefficients(&ctx, 27)).unwrap();
        let alphas = vec![ctx.from_u64(2), ctx.from_u64(5)];
        let tail = vec![ctx.from_u64(9)];

        let restricted = poly.restrict_prefix(&ctx, &alphas).unwrap();
        let mut full_point = alphas;
        full_point.extend_from_slice(&tail);

        assert_eq!(
            restricted.evaluate(&ctx, &tail).unwrap(),
            poly.evaluate(&ctx, &full_point).unwrap()
        );
    }

    #[test]
    fn multilinear_embedding_should_insert_zero_quadratic_coefficients() {
        let ctx = sample_context();
        let coefficients = vec![
            ctx.from_u64(2),
            ctx.from_u64(3),
            ctx.from_u64(5),
            ctx.from_u64(7),
        ];
        let multilinear = MultilinearPolynomial::new(2, coefficients.clone()).unwrap();

        let embedded = multilinear.to_multi_quadratic(&ctx).unwrap();

        assert_eq!(pow2_checked(4).unwrap(), 16);
        assert_eq!(embedded.variable_count(), 2);
        assert_eq!(embedded.coefficients().len(), 5);
        assert_eq!(embedded.coefficients()[0], coefficients[0]);
        assert_eq!(embedded.coefficients()[1], coefficients[1]);
        assert_eq!(embedded.coefficients()[2], ctx.zero());
        assert_eq!(embedded.coefficients()[3], coefficients[2]);
        assert_eq!(embedded.coefficients()[4], coefficients[3]);
    }

    #[test]
    fn multilinear_evaluation_should_match_embedded_evaluation() {
        let ctx = sample_context();
        let coefficients = vec![
            ctx.from_u64(2),
            ctx.from_u64(3),
            ctx.from_u64(5),
            ctx.from_u64(7),
        ];
        let multilinear = MultilinearPolynomial::new(2, coefficients).unwrap();
        let embedded = multilinear.to_multi_quadratic(&ctx).unwrap();
        let point = vec![ctx.from_u64(5), ctx.from_u64(7)];

        assert_eq!(
            multilinear.evaluate(&ctx, &point).unwrap(),
            embedded.evaluate(&ctx, &point).unwrap()
        );
    }

    #[test]
    fn evaluate_pow_should_work_on_ternary_domain_points() {
        let ctx = GrContext::new(GrConfig {
            p: 2,
            k_exp: 16,
            r: 18,
        })
        .unwrap();
        let domain = Domain::teichmuller_subgroup(std::sync::Arc::new(ctx.clone()), 27).unwrap();
        let poly = MultiQuadraticPolynomial::new(3, sample_coefficients(&ctx, 27)).unwrap();

        for index in 0..domain.size() {
            let point = domain.element(index).unwrap();
            assert_eq!(
                poly.evaluate(&ctx, &pow_m(&ctx, &point, 3).unwrap())
                    .unwrap(),
                poly.evaluate_pow(&ctx, &point)
            );
        }
    }

    #[test]
    fn invalid_inputs_should_reject() {
        let ctx = sample_context();

        assert!(matches!(
            MultiQuadraticPolynomial::new(2, vec![ctx.one(); 10]),
            Err(GrError::InvalidPolynomial(_))
        ));
        assert!(matches!(
            MultiQuadraticPolynomial::new(2, sample_coefficients(&ctx, 9))
                .unwrap()
                .evaluate(&ctx, &[ctx.one()]),
            Err(GrError::InvalidPolynomial(_))
        ));
        assert!(matches!(
            encode_base3_index(&[0, 3]),
            Err(GrError::InvalidPolynomial(_))
        ));
        assert!(matches!(
            decode_base3_index(9, 2),
            Err(GrError::IndexOutOfRange { .. })
        ));
        assert!(matches!(
            pow3_checked(41),
            Err(GrError::ArithmeticOverflow(_))
        ));
    }
}
