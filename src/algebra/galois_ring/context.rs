use core::fmt;

use ark_std::rand::Rng;
use serde::{Deserialize, Serialize};

use super::{poly_f2::BinaryPolynomial, GrElem};

pub type Result<T> = core::result::Result<T, GrError>;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct GrConfig {
    pub p: u64,
    pub k_exp: u32,
    pub r: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GrError {
    UnsupportedPrime(u64),
    ZeroPrecision,
    UnsupportedPrecision(u32),
    ZeroDegree,
    CoefficientLength { expected: usize, actual: usize },
    DeserializeSize { expected: usize, actual: usize },
    NonUnit,
    InvalidSubgroupSize { size: u64 },
    InvalidDomain(&'static str),
    InvalidPolynomial(&'static str),
    ArithmeticOverflow(&'static str),
    IndexOutOfRange { index: u64, size: u64 },
    DifferentRings,
    NoIrreduciblePolynomial { degree: usize, attempts: usize },
    InvalidDefiningPolynomial(&'static str),
}

impl fmt::Display for GrError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedPrime(p) => {
                write!(formatter, "only p = 2 is supported for WHIR_GR P1, got p = {p}")
            }
            Self::ZeroPrecision => formatter.write_str("k_exp must be greater than zero"),
            Self::UnsupportedPrecision(k_exp) => {
                write!(formatter, "k_exp must be at most 64, got {k_exp}")
            }
            Self::ZeroDegree => formatter.write_str("extension degree r must be greater than zero"),
            Self::CoefficientLength { expected, actual } => write!(
                formatter,
                "coefficient length mismatch: expected {expected}, got {actual}"
            ),
            Self::DeserializeSize { expected, actual } => write!(
                formatter,
                "serialized element size mismatch: expected {expected}, got {actual}"
            ),
            Self::NonUnit => formatter.write_str("operation requires a unit element"),
            Self::InvalidSubgroupSize { size } => {
                write!(formatter, "invalid Teichmuller subgroup size {size}")
            }
            Self::InvalidDomain(message) => write!(formatter, "invalid domain: {message}"),
            Self::InvalidPolynomial(message) => {
                write!(formatter, "invalid polynomial: {message}")
            }
            Self::ArithmeticOverflow(message) => write!(formatter, "arithmetic overflow: {message}"),
            Self::IndexOutOfRange { index, size } => {
                write!(formatter, "domain index {index} is out of range for size {size}")
            }
            Self::DifferentRings => formatter.write_str("operation requires matching ring contexts"),
            Self::NoIrreduciblePolynomial { degree, attempts } => write!(
                formatter,
                "failed to find irreducible binary polynomial of degree {degree} after {attempts} attempts"
            ),
            Self::InvalidDefiningPolynomial(message) => {
                write!(formatter, "invalid defining polynomial: {message}")
            }
        }
    }
}

impl std::error::Error for GrError {}

#[derive(Clone, Debug)]
pub struct GrContext {
    config: GrConfig,
    coeff_bytes: usize,
    mask: u64,
    defining_polynomial: BinaryPolynomial,
    reduction_coefficients: Vec<u8>,
}

impl GrContext {
    pub fn new(config: GrConfig) -> Result<Self> {
        if config.p != 2 {
            return Err(GrError::UnsupportedPrime(config.p));
        }
        if config.k_exp == 0 {
            return Err(GrError::ZeroPrecision);
        }
        if config.k_exp > u64::BITS {
            return Err(GrError::UnsupportedPrecision(config.k_exp));
        }
        if config.r == 0 {
            return Err(GrError::ZeroDegree);
        }

        let defining_polynomial = BinaryPolynomial::irreducible_for_degree(config.r)?;
        let reduction_coefficients = defining_polynomial.low_coefficients(config.r);
        let coeff_bytes = config.k_exp.div_ceil(8) as usize;
        let mask = mask_for_precision(config.k_exp);

        Ok(Self {
            config,
            coeff_bytes,
            mask,
            defining_polynomial,
            reduction_coefficients,
        })
    }

    pub const fn config(&self) -> &GrConfig {
        &self.config
    }

    pub const fn coeff_bytes(&self) -> usize {
        self.coeff_bytes
    }

    pub const fn elem_bytes(&self) -> usize {
        self.coeff_bytes * self.config.r
    }

    pub const fn mul_scratch_len(&self) -> usize {
        self.config.r.saturating_mul(2).saturating_sub(1)
    }

    pub const fn modulus(&self) -> u128 {
        1u128 << self.config.k_exp
    }

    pub fn defining_polynomial_low_coefficients(&self) -> &[u8] {
        &self.reduction_coefficients
    }

    pub fn zero(&self) -> GrElem {
        GrElem::new_unchecked(vec![0; self.config.r])
    }

    pub fn one(&self) -> GrElem {
        self.from_u64(1)
    }

    pub fn x(&self) -> GrElem {
        let mut coefficients = vec![0; self.config.r];
        if self.config.r == 1 {
            coefficients[0] = self.neg_u64(self.reduction_coefficients[0].into());
        } else {
            coefficients[1] = 1;
        }
        GrElem::new_unchecked(coefficients)
    }

    pub fn from_u64(&self, value: u64) -> GrElem {
        let mut coefficients = vec![0; self.config.r];
        coefficients[0] = self.normalize(value);
        GrElem::new_unchecked(coefficients)
    }

    pub fn from_coefficients(&self, coefficients: &[u64]) -> Result<GrElem> {
        if coefficients.len() != self.config.r {
            return Err(GrError::CoefficientLength {
                expected: self.config.r,
                actual: coefficients.len(),
            });
        }
        Ok(GrElem::new_unchecked(
            coefficients
                .iter()
                .map(|&coefficient| self.normalize(coefficient))
                .collect(),
        ))
    }

    pub fn random_element(&self, rng: &mut impl Rng) -> GrElem {
        GrElem::new_unchecked(
            (0..self.config.r)
                .map(|_| self.normalize(rng.gen()))
                .collect(),
        )
    }

    pub fn add(&self, lhs: &GrElem, rhs: &GrElem) -> GrElem {
        self.debug_assert_element(lhs);
        self.debug_assert_element(rhs);
        GrElem::new_unchecked(
            lhs.coefficients()
                .iter()
                .zip(rhs.coefficients())
                .map(|(&lhs, &rhs)| self.normalize(lhs.wrapping_add(rhs)))
                .collect(),
        )
    }

    pub fn add_into(&self, out: &mut GrElem, lhs: &GrElem, rhs: &GrElem) {
        self.debug_assert_element(out);
        self.debug_assert_element(lhs);
        self.debug_assert_element(rhs);
        for ((out, &lhs), &rhs) in out
            .coefficients_mut()
            .iter_mut()
            .zip(lhs.coefficients())
            .zip(rhs.coefficients())
        {
            *out = self.normalize(lhs.wrapping_add(rhs));
        }
    }

    pub fn add_assign(&self, lhs: &mut GrElem, rhs: &GrElem) {
        self.debug_assert_element(lhs);
        self.debug_assert_element(rhs);
        for (lhs, &rhs) in lhs.coefficients_mut().iter_mut().zip(rhs.coefficients()) {
            *lhs = self.normalize(lhs.wrapping_add(rhs));
        }
    }

    pub fn sub(&self, lhs: &GrElem, rhs: &GrElem) -> GrElem {
        self.debug_assert_element(lhs);
        self.debug_assert_element(rhs);
        GrElem::new_unchecked(
            lhs.coefficients()
                .iter()
                .zip(rhs.coefficients())
                .map(|(&lhs, &rhs)| self.normalize(lhs.wrapping_sub(rhs)))
                .collect(),
        )
    }

    pub fn sub_into(&self, out: &mut GrElem, lhs: &GrElem, rhs: &GrElem) {
        self.debug_assert_element(out);
        self.debug_assert_element(lhs);
        self.debug_assert_element(rhs);
        for ((out, &lhs), &rhs) in out
            .coefficients_mut()
            .iter_mut()
            .zip(lhs.coefficients())
            .zip(rhs.coefficients())
        {
            *out = self.normalize(lhs.wrapping_sub(rhs));
        }
    }

    pub fn neg(&self, value: &GrElem) -> GrElem {
        self.debug_assert_element(value);
        GrElem::new_unchecked(
            value
                .coefficients()
                .iter()
                .map(|&coefficient| self.neg_u64(coefficient))
                .collect(),
        )
    }

    pub fn mul(&self, lhs: &GrElem, rhs: &GrElem) -> GrElem {
        self.debug_assert_element(lhs);
        self.debug_assert_element(rhs);

        let mut coefficients: Vec<u64> = vec![0; self.config.r.saturating_mul(2).saturating_sub(1)];
        for (lhs_index, &lhs_coefficient) in lhs.coefficients().iter().enumerate() {
            for (rhs_index, &rhs_coefficient) in rhs.coefficients().iter().enumerate() {
                let product = self.mul_coeff(lhs_coefficient, rhs_coefficient);
                let index = lhs_index + rhs_index;
                coefficients[index] = self.normalize(coefficients[index].wrapping_add(product));
            }
        }

        self.reduce_coefficients(coefficients)
    }

    pub fn mul_into(&self, out: &mut GrElem, lhs: &GrElem, rhs: &GrElem, scratch: &mut [u64]) {
        self.debug_assert_element(out);
        self.debug_assert_element(lhs);
        self.debug_assert_element(rhs);
        debug_assert!(scratch.len() >= self.mul_scratch_len());

        let scratch_len = self.mul_scratch_len();
        let scratch = &mut scratch[..scratch_len];
        scratch.fill(0);
        for (lhs_index, &lhs_coefficient) in lhs.coefficients().iter().enumerate() {
            for (rhs_index, &rhs_coefficient) in rhs.coefficients().iter().enumerate() {
                let product = self.mul_coeff(lhs_coefficient, rhs_coefficient);
                let index = lhs_index + rhs_index;
                scratch[index] = self.normalize(scratch[index].wrapping_add(product));
            }
        }

        for degree in (self.config.r..scratch.len()).rev() {
            let coefficient = scratch[degree];
            if coefficient == 0 {
                continue;
            }
            scratch[degree] = 0;
            let offset = degree - self.config.r;
            for (index, &defining_coefficient) in self.reduction_coefficients.iter().enumerate() {
                if defining_coefficient == 1 {
                    let target = offset + index;
                    scratch[target] = self.normalize(scratch[target].wrapping_sub(coefficient));
                }
            }
        }

        for (out, &coefficient) in out
            .coefficients_mut()
            .iter_mut()
            .zip(scratch.iter().take(self.config.r))
        {
            *out = self.normalize(coefficient);
        }
    }

    pub fn mul_base_scalar_into(&self, out: &mut GrElem, value: &GrElem, scalar: u64) {
        self.debug_assert_element(out);
        self.debug_assert_element(value);
        let scalar = self.normalize(scalar);
        for (out, &coefficient) in out.coefficients_mut().iter_mut().zip(value.coefficients()) {
            *out = self.mul_coeff(coefficient, scalar);
        }
    }

    pub fn square(&self, value: &GrElem) -> GrElem {
        self.mul(value, value)
    }

    pub fn square_into(&self, out: &mut GrElem, value: &GrElem, scratch: &mut [u64]) {
        self.mul_into(out, value, value, scratch);
    }

    pub fn pow(&self, base: &GrElem, mut exponent: u128) -> GrElem {
        let mut result = self.one();
        let mut power = base.clone();
        while exponent != 0 {
            if exponent & 1 == 1 {
                result = self.mul(&result, &power);
            }
            exponent >>= 1;
            if exponent != 0 {
                power = self.square(&power);
            }
        }
        result
    }

    pub fn pow_words_le(&self, base: &GrElem, exponent_words: &[u64]) -> GrElem {
        let mut result = self.one();
        let mut power = base.clone();
        for &word in exponent_words {
            let mut remaining = word;
            for _ in 0..u64::BITS {
                if remaining & 1 == 1 {
                    result = self.mul(&result, &power);
                }
                remaining >>= 1;
                power = self.square(&power);
            }
        }
        result
    }

    pub fn is_unit(&self, value: &GrElem) -> bool {
        self.debug_assert_element(value);
        value
            .coefficients()
            .iter()
            .any(|&coefficient| coefficient & 1 == 1)
    }

    pub fn inv(&self, value: &GrElem) -> Result<GrElem> {
        self.debug_assert_element(value);
        if !self.is_unit(value) {
            return Err(GrError::NonUnit);
        }

        let value_mod2 = BinaryPolynomial::from_coefficients_mod2(value.coefficients());
        let inverse_mod2 = value_mod2.inverse_mod(&self.defining_polynomial)?;
        let inverse_coefficients = (0..self.config.r)
            .map(|index| u64::from(inverse_mod2.coefficient(index)))
            .collect::<Vec<_>>();
        let mut inverse = GrElem::new_unchecked(inverse_coefficients);

        let two = self.from_u64(2);
        let mut precision = 1;
        while precision < self.config.k_exp {
            let value_times_inverse = self.mul(value, &inverse);
            let correction = self.sub(&two, &value_times_inverse);
            inverse = self.mul(&inverse, &correction);
            precision = precision.saturating_mul(2);
        }

        Ok(inverse)
    }

    pub fn batch_inv(&self, values: &[GrElem]) -> Result<Vec<GrElem>> {
        if values.is_empty() {
            return Ok(Vec::new());
        }

        let mut prefix_products = Vec::with_capacity(values.len());
        let mut running_product = self.one();
        for value in values {
            if !self.is_unit(value) {
                return Err(GrError::NonUnit);
            }
            prefix_products.push(running_product.clone());
            running_product = self.mul(&running_product, value);
        }

        let mut suffix_inverse = self.inv(&running_product)?;
        let mut inverses = vec![self.zero(); values.len()];
        for index in (0..values.len()).rev() {
            inverses[index] = self.mul(&prefix_products[index], &suffix_inverse);
            suffix_inverse = self.mul(&suffix_inverse, &values[index]);
        }

        Ok(inverses)
    }

    pub fn serialize(&self, value: &GrElem) -> Vec<u8> {
        self.debug_assert_element(value);
        let mut out = vec![0; self.elem_bytes()];
        self.serialize_into(&mut out, value);
        out
    }

    pub fn serialize_into(&self, out: &mut [u8], value: &GrElem) {
        self.debug_assert_element(value);
        debug_assert_eq!(out.len(), self.elem_bytes());
        for (coefficient_index, &coefficient) in value.coefficients().iter().enumerate() {
            let offset = coefficient_index * self.coeff_bytes;
            out[offset..offset + self.coeff_bytes]
                .copy_from_slice(&coefficient.to_le_bytes()[..self.coeff_bytes]);
        }
    }

    pub fn deserialize(&self, bytes: &[u8]) -> Result<GrElem> {
        if bytes.len() != self.elem_bytes() {
            return Err(GrError::DeserializeSize {
                expected: self.elem_bytes(),
                actual: bytes.len(),
            });
        }

        let coefficients = bytes
            .chunks_exact(self.coeff_bytes)
            .map(|chunk| {
                let mut buffer = [0; 8];
                buffer[..chunk.len()].copy_from_slice(chunk);
                self.normalize(u64::from_le_bytes(buffer))
            })
            .collect();
        Ok(GrElem::new_unchecked(coefficients))
    }

    fn reduce_coefficients(&self, mut coefficients: Vec<u64>) -> GrElem {
        coefficients.resize(self.config.r.saturating_mul(2).saturating_sub(1), 0);
        for degree in (self.config.r..coefficients.len()).rev() {
            let coefficient = coefficients[degree];
            if coefficient == 0 {
                continue;
            }
            coefficients[degree] = 0;
            let offset = degree - self.config.r;
            for (index, &defining_coefficient) in self.reduction_coefficients.iter().enumerate() {
                if defining_coefficient == 1 {
                    let target = offset + index;
                    coefficients[target] =
                        self.normalize(coefficients[target].wrapping_sub(coefficient));
                }
            }
        }
        coefficients.truncate(self.config.r);
        for coefficient in &mut coefficients {
            *coefficient = self.normalize(*coefficient);
        }
        GrElem::new_unchecked(coefficients)
    }

    const fn normalize(&self, value: u64) -> u64 {
        value & self.mask
    }

    const fn neg_u64(&self, value: u64) -> u64 {
        self.normalize(0u64.wrapping_sub(value))
    }

    fn mul_coeff(&self, lhs: u64, rhs: u64) -> u64 {
        if self.config.k_exp == u64::BITS {
            lhs.wrapping_mul(rhs)
        } else {
            ((u128::from(lhs) * u128::from(rhs)) & u128::from(self.mask)) as u64
        }
    }

    fn debug_assert_element(&self, value: &GrElem) {
        debug_assert_eq!(value.coefficients().len(), self.config.r);
    }
}

const fn mask_for_precision(k_exp: u32) -> u64 {
    if k_exp == u64::BITS {
        u64::MAX
    } else {
        (1u64 << k_exp) - 1
    }
}

#[cfg(test)]
mod tests {
    use ark_std::{rand::SeedableRng, test_rng};

    use super::{GrConfig, GrContext, GrError};

    fn small_context() -> GrContext {
        GrContext::new(GrConfig {
            p: 2,
            k_exp: 16,
            r: 2,
        })
        .unwrap()
    }

    #[test]
    fn metadata_should_match_expected_byte_widths() {
        let ctx = GrContext::new(GrConfig {
            p: 2,
            k_exp: 16,
            r: 54,
        })
        .unwrap();

        assert_eq!(ctx.coeff_bytes(), 2);
        assert_eq!(ctx.elem_bytes(), 108);
    }

    #[test]
    fn constructor_should_reject_unsupported_parameters() {
        assert!(matches!(
            GrContext::new(GrConfig {
                p: 3,
                k_exp: 16,
                r: 2,
            }),
            Err(GrError::UnsupportedPrime(3))
        ));
        assert!(matches!(
            GrContext::new(GrConfig {
                p: 2,
                k_exp: 0,
                r: 2,
            }),
            Err(GrError::ZeroPrecision)
        ));
        assert!(matches!(
            GrContext::new(GrConfig {
                p: 2,
                k_exp: 16,
                r: 0,
            }),
            Err(GrError::ZeroDegree)
        ));
    }

    #[test]
    fn serialization_should_round_trip() {
        let ctx = small_context();
        let value = ctx.from_coefficients(&[0x1234, 0xFEDC]).unwrap();

        let encoded = ctx.serialize(&value);
        let decoded = ctx.deserialize(&encoded).unwrap();

        assert_eq!(decoded, value);
    }

    #[test]
    fn deserialize_should_reject_wrong_size() {
        let ctx = small_context();

        let result = ctx.deserialize(&[1, 2, 3]);

        assert!(matches!(result, Err(GrError::DeserializeSize { .. })));
    }

    #[test]
    fn multiplication_should_reduce_by_defining_polynomial() {
        let ctx = small_context();
        let x = ctx.x();

        let x_squared = ctx.square(&x);

        assert_eq!(
            x_squared,
            ctx.from_coefficients(&[u16::MAX.into(), u16::MAX.into()])
                .unwrap()
        );
    }

    #[test]
    fn in_place_arithmetic_should_match_allocating_apis() {
        let ctx = small_context();
        let lhs = ctx.from_coefficients(&[0x1234, 0xFEDC]).unwrap();
        let rhs = ctx.from_coefficients(&[0x0102, 0x0304]).unwrap();
        let mut out = ctx.zero();
        let mut scratch = vec![0; ctx.mul_scratch_len()];

        ctx.add_into(&mut out, &lhs, &rhs);
        assert_eq!(out, ctx.add(&lhs, &rhs));

        let mut assigned = lhs.clone();
        ctx.add_assign(&mut assigned, &rhs);
        assert_eq!(assigned, ctx.add(&lhs, &rhs));

        ctx.sub_into(&mut out, &lhs, &rhs);
        assert_eq!(out, ctx.sub(&lhs, &rhs));

        ctx.mul_into(&mut out, &lhs, &rhs, &mut scratch);
        assert_eq!(out, ctx.mul(&lhs, &rhs));

        ctx.mul_base_scalar_into(&mut out, &lhs, 17);
        assert_eq!(out, ctx.mul(&lhs, &ctx.from_u64(17)));

        ctx.square_into(&mut out, &lhs, &mut scratch);
        assert_eq!(out, ctx.square(&lhs));

        let mut serialized = vec![0; ctx.elem_bytes()];
        ctx.serialize_into(&mut serialized, &lhs);
        assert_eq!(serialized, vec![0x34, 0x12, 0xdc, 0xfe]);
    }

    #[test]
    fn inverse_should_multiply_unit_to_one() {
        let ctx = small_context();
        let unit = ctx.from_coefficients(&[1, 1]).unwrap();

        let inverse = ctx.inv(&unit).unwrap();
        let product = ctx.mul(&unit, &inverse);

        assert_eq!(product, ctx.one());
    }

    #[test]
    fn inverse_should_reject_non_unit() {
        let ctx = small_context();
        let non_unit = ctx.from_coefficients(&[2, 0]).unwrap();

        let result = ctx.inv(&non_unit);

        assert_eq!(result, Err(GrError::NonUnit));
    }

    #[test]
    fn batch_inverse_should_multiply_units_to_one() {
        let ctx = small_context();
        let values = vec![
            ctx.from_coefficients(&[1, 1]).unwrap(),
            ctx.from_coefficients(&[3, 2]).unwrap(),
            ctx.from_coefficients(&[5, 7]).unwrap(),
        ];

        let inverses = ctx.batch_inv(&values).unwrap();

        assert_eq!(inverses.len(), values.len());
        for (value, inverse) in values.iter().zip(&inverses) {
            assert_eq!(ctx.mul(value, inverse), ctx.one());
        }
    }

    #[test]
    fn ring_axioms_should_hold_for_sampled_elements() {
        let ctx = GrContext::new(GrConfig {
            p: 2,
            k_exp: 8,
            r: 3,
        })
        .unwrap();
        let mut rng = test_rng();
        let samples = (0..8)
            .map(|_| ctx.random_element(&mut rng))
            .collect::<Vec<_>>();

        for a in &samples {
            assert_eq!(ctx.add(a, &ctx.zero()), *a);
            assert_eq!(ctx.mul(a, &ctx.one()), *a);
            assert_eq!(ctx.add(a, &ctx.neg(a)), ctx.zero());
            for b in &samples {
                assert_eq!(ctx.add(a, b), ctx.add(b, a));
                assert_eq!(ctx.mul(a, b), ctx.mul(b, a));
                for c in &samples {
                    assert_eq!(ctx.add(&ctx.add(a, b), c), ctx.add(a, &ctx.add(b, c)));
                    assert_eq!(ctx.mul(&ctx.mul(a, b), c), ctx.mul(a, &ctx.mul(b, c)));
                    assert_eq!(
                        ctx.mul(a, &ctx.add(b, c)),
                        ctx.add(&ctx.mul(a, b), &ctx.mul(a, c))
                    );
                }
            }
        }
    }

    #[test]
    fn random_element_should_be_deterministic_for_seeded_rng() {
        let ctx = small_context();
        let mut lhs_rng = ark_std::rand::rngs::StdRng::seed_from_u64(7);
        let mut rhs_rng = ark_std::rand::rngs::StdRng::seed_from_u64(7);

        let lhs = ctx.random_element(&mut lhs_rng);
        let rhs = ctx.random_element(&mut rhs_rng);

        assert_eq!(lhs, rhs);
    }
}
