use super::{GrError, Result};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct BinaryPolynomial {
    bits: Vec<u64>,
}

impl BinaryPolynomial {
    pub(super) fn from_exponents(exponents: &[usize]) -> Self {
        let mut polynomial = Self::zero();
        for &exponent in exponents {
            polynomial.set_bit(exponent);
        }
        polynomial
    }

    pub(super) fn from_coefficients_mod2(coefficients: &[u64]) -> Self {
        let mut polynomial = Self::zero();
        for (index, &coefficient) in coefficients.iter().enumerate() {
            if coefficient & 1 == 1 {
                polynomial.set_bit(index);
            }
        }
        polynomial
    }

    pub(super) fn irreducible_for_degree(degree: usize) -> Result<Self> {
        if degree == 0 {
            return Err(GrError::ZeroDegree);
        }
        if degree == 1 {
            return Ok(Self::from_exponents(&[0, 1]));
        }

        for middle in 1..degree {
            let candidate = Self::from_exponents(&[0, middle, degree]);
            if candidate.is_irreducible() {
                return Ok(candidate);
            }
        }

        let mut state = deterministic_seed(degree);
        for _ in 0..8192 {
            let candidate = random_monic_candidate(degree, &mut state);
            if candidate.is_irreducible() {
                return Ok(candidate);
            }
        }

        Err(GrError::NoIrreduciblePolynomial {
            degree,
            attempts: 8192 + degree.saturating_sub(1),
        })
    }

    pub(super) fn degree(&self) -> Option<usize> {
        for (word_index, &word) in self.bits.iter().enumerate().rev() {
            if word != 0 {
                let high_bit = u64::BITS - 1 - word.leading_zeros();
                return Some(word_index * u64::BITS as usize + high_bit as usize);
            }
        }
        None
    }

    pub(super) fn coefficient(&self, degree: usize) -> bool {
        let word_index = degree / u64::BITS as usize;
        let bit_index = degree % u64::BITS as usize;
        self.bits
            .get(word_index)
            .is_some_and(|word| (word >> bit_index) & 1 == 1)
    }

    pub(super) fn low_coefficients(&self, degree: usize) -> Vec<u8> {
        (0..degree)
            .map(|index| u8::from(self.coefficient(index)))
            .collect()
    }

    pub(super) fn inverse_mod(&self, modulus: &Self) -> Result<Self> {
        let modulus_degree = modulus.degree().ok_or(GrError::InvalidDefiningPolynomial(
            "modulus must be nonzero",
        ))?;
        if !modulus.coefficient(modulus_degree) {
            return Err(GrError::InvalidDefiningPolynomial("modulus must be monic"));
        }

        let mut r0 = modulus.clone();
        let mut r1 = self.rem(modulus);
        let mut t0 = Self::zero();
        let mut t1 = Self::one();

        while !r1.is_zero() {
            let (quotient, remainder) = r0.div_rem(&r1)?;
            r0 = r1;
            r1 = remainder;
            let product = quotient.mul(&t1);
            let next_t = t0.xor(&product);
            t0 = t1;
            t1 = next_t;
        }

        if !r0.is_one() {
            return Err(GrError::NonUnit);
        }

        Ok(t0.rem(modulus))
    }

    const fn zero() -> Self {
        Self { bits: Vec::new() }
    }

    fn one() -> Self {
        Self::from_exponents(&[0])
    }

    fn x() -> Self {
        Self::from_exponents(&[1])
    }

    fn is_zero(&self) -> bool {
        self.degree().is_none()
    }

    fn is_one(&self) -> bool {
        self.degree() == Some(0) && self.coefficient(0)
    }

    fn set_bit(&mut self, degree: usize) {
        let word_index = degree / u64::BITS as usize;
        let bit_index = degree % u64::BITS as usize;
        if self.bits.len() <= word_index {
            self.bits.resize(word_index + 1, 0);
        }
        self.bits[word_index] |= 1u64 << bit_index;
    }

    fn xor_assign_shifted(&mut self, rhs: &Self, shift: usize) {
        if rhs.is_zero() {
            return;
        }

        let word_shift = shift / u64::BITS as usize;
        let bit_shift = shift % u64::BITS as usize;
        let required_words = rhs.bits.len() + word_shift + usize::from(bit_shift != 0);
        if self.bits.len() < required_words {
            self.bits.resize(required_words, 0);
        }

        for (index, &word) in rhs.bits.iter().enumerate() {
            self.bits[index + word_shift] ^= word << bit_shift;
            if bit_shift != 0 {
                self.bits[index + word_shift + 1] ^= word >> (u64::BITS as usize - bit_shift);
            }
        }
        self.trim();
    }

    fn xor(&self, rhs: &Self) -> Self {
        let mut out = self.clone();
        if out.bits.len() < rhs.bits.len() {
            out.bits.resize(rhs.bits.len(), 0);
        }
        for (index, &word) in rhs.bits.iter().enumerate() {
            out.bits[index] ^= word;
        }
        out.trim();
        out
    }

    fn mul(&self, rhs: &Self) -> Self {
        let mut out = Self::zero();
        for (word_index, &word) in self.bits.iter().enumerate() {
            let mut remaining = word;
            while remaining != 0 {
                let bit_index = remaining.trailing_zeros() as usize;
                out.xor_assign_shifted(rhs, word_index * u64::BITS as usize + bit_index);
                remaining &= remaining - 1;
            }
        }
        out
    }

    fn mul_mod(&self, rhs: &Self, modulus: &Self) -> Self {
        self.mul(rhs).rem(modulus)
    }

    fn square_mod(&self, modulus: &Self) -> Self {
        self.mul_mod(self, modulus)
    }

    fn rem(&self, modulus: &Self) -> Self {
        let Some(modulus_degree) = modulus.degree() else {
            return self.clone();
        };
        let mut out = self.clone();
        while let Some(out_degree) = out.degree() {
            if out_degree < modulus_degree {
                break;
            }
            out.xor_assign_shifted(modulus, out_degree - modulus_degree);
        }
        out
    }

    fn div_rem(&self, divisor: &Self) -> Result<(Self, Self)> {
        let divisor_degree = divisor
            .degree()
            .ok_or(GrError::InvalidDefiningPolynomial("division by zero"))?;
        let mut quotient = Self::zero();
        let mut remainder = self.clone();

        while let Some(remainder_degree) = remainder.degree() {
            if remainder_degree < divisor_degree {
                break;
            }
            let shift = remainder_degree - divisor_degree;
            quotient.set_bit(shift);
            remainder.xor_assign_shifted(divisor, shift);
        }

        Ok((quotient, remainder))
    }

    fn gcd(mut lhs: Self, mut rhs: Self) -> Self {
        while !rhs.is_zero() {
            let remainder = lhs.rem(&rhs);
            lhs = rhs;
            rhs = remainder;
        }
        lhs
    }

    fn frobenius_power_x(exponent: usize, modulus: &Self) -> Self {
        let mut value = Self::x();
        for _ in 0..exponent {
            value = value.square_mod(modulus);
        }
        value
    }

    fn is_irreducible(&self) -> bool {
        let Some(degree) = self.degree() else {
            return false;
        };
        if degree == 0 || !self.coefficient(0) {
            return false;
        }

        for prime_divisor in prime_divisors(degree) {
            let exponent = degree / prime_divisor;
            let powered = Self::frobenius_power_x(exponent, self);
            let difference = powered.xor(&Self::x());
            if !Self::gcd(difference, self.clone()).is_one() {
                return false;
            }
        }

        Self::frobenius_power_x(degree, self) == Self::x()
    }

    fn trim(&mut self) {
        while self.bits.last().is_some_and(|word| *word == 0) {
            self.bits.pop();
        }
    }
}

fn random_monic_candidate(degree: usize, state: &mut u64) -> BinaryPolynomial {
    let mut candidate = BinaryPolynomial::from_exponents(&[0, degree]);
    for index in 1..degree {
        if splitmix64(state) & 1 == 1 {
            candidate.set_bit(index);
        }
    }
    candidate
}

const fn deterministic_seed(degree: usize) -> u64 {
    0x5758_4749_5252_4544u64 ^ degree as u64
}

const fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn prime_divisors(mut value: usize) -> Vec<usize> {
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
    use super::BinaryPolynomial;

    #[test]
    fn irreducibility_should_accept_x2_plus_x_plus_one() {
        let polynomial = BinaryPolynomial::from_exponents(&[0, 1, 2]);

        assert!(polynomial.is_irreducible());
    }

    #[test]
    fn inverse_mod_should_invert_nonzero_field_element() {
        let modulus = BinaryPolynomial::from_exponents(&[0, 1, 2]);
        let element = BinaryPolynomial::from_exponents(&[0, 1]);

        let inverse = element.inverse_mod(&modulus).unwrap();
        let product = element.mul_mod(&inverse, &modulus);

        assert_eq!(product, BinaryPolynomial::one());
    }
}
