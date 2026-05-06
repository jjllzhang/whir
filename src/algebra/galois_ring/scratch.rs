use super::{GrContext, GrElem};

#[derive(Debug)]
pub struct GrScratch {
    elements: Vec<GrElem>,
    mul_scratch: Vec<u64>,
}

impl GrScratch {
    pub fn with_elements(ctx: &GrContext, element_count: usize) -> Self {
        Self {
            elements: vec![ctx.zero(); element_count],
            mul_scratch: vec![0; ctx.mul_scratch_len()],
        }
    }

    pub fn parts_mut(&mut self) -> (&mut [GrElem], &mut [u64]) {
        (&mut self.elements, &mut self.mul_scratch)
    }
}

pub fn clear_elem(value: &mut GrElem) {
    value.coefficients_mut().fill(0);
}

#[cfg(test)]
mod tests {
    use super::{clear_elem, GrScratch};
    use crate::algebra::galois_ring::{GrConfig, GrContext};

    #[test]
    fn scratch_parts_should_support_in_place_arithmetic() {
        let ctx = GrContext::new(GrConfig {
            p: 2,
            k_exp: 16,
            r: 3,
        })
        .unwrap();
        let lhs = ctx.from_coefficients(&[0x1234, 0x4321, 0x2222]).unwrap();
        let rhs = ctx.from_coefficients(&[0x0102, 0x0304, 0x0506]).unwrap();
        let mut scratch = GrScratch::with_elements(&ctx, 2);
        let (elements, mul_scratch) = scratch.parts_mut();
        let [product, square] = elements else {
            unreachable!("scratch was created with exactly two elements");
        };

        ctx.mul_into(product, &lhs, &rhs, mul_scratch);
        assert_eq!(*product, ctx.mul(&lhs, &rhs));

        ctx.square_into(square, &lhs, mul_scratch);
        assert_eq!(*square, ctx.square(&lhs));

        clear_elem(product);
        assert_eq!(*product, ctx.zero());
    }
}
