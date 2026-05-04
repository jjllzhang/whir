use std::sync::Arc;

use super::{
    has_exact_multiplicative_order, is_teichmuller_element, teichmuller_subgroup_generator,
    GrContext, GrElem, GrError, Result,
};

#[derive(Clone, Debug)]
pub struct Domain {
    ctx: Arc<GrContext>,
    offset: GrElem,
    root: GrElem,
    size: u64,
}

impl Domain {
    pub fn teichmuller_subgroup(ctx: Arc<GrContext>, size: u64) -> Result<Self> {
        let offset = ctx.one();
        let root = teichmuller_subgroup_generator(&ctx, size)?;
        Self::new(ctx, offset, root, size)
    }

    pub fn teichmuller_coset(ctx: Arc<GrContext>, offset: GrElem, size: u64) -> Result<Self> {
        let root = teichmuller_subgroup_generator(&ctx, size)?;
        Self::new(ctx, offset, root, size)
    }

    pub fn new(ctx: Arc<GrContext>, offset: GrElem, root: GrElem, size: u64) -> Result<Self> {
        if size == 0 {
            return Err(GrError::InvalidDomain("size must be greater than zero"));
        }
        if !ctx.is_unit(&offset) {
            return Err(GrError::InvalidDomain("offset must be a unit"));
        }
        if !ctx.is_unit(&root) {
            return Err(GrError::InvalidDomain("root must be a unit"));
        }
        if !has_exact_multiplicative_order(&ctx, &root, size) {
            return Err(GrError::InvalidDomain(
                "root must have exact multiplicative order equal to size",
            ));
        }

        Ok(Self {
            ctx,
            offset,
            root,
            size,
        })
    }

    pub const fn context(&self) -> &Arc<GrContext> {
        &self.ctx
    }

    pub const fn size(&self) -> u64 {
        self.size
    }

    pub const fn offset(&self) -> &GrElem {
        &self.offset
    }

    pub const fn root(&self) -> &GrElem {
        &self.root
    }

    pub fn element(&self, index: u64) -> Result<GrElem> {
        if index >= self.size {
            return Err(GrError::IndexOutOfRange {
                index,
                size: self.size,
            });
        }

        Ok(self
            .ctx
            .mul(&self.offset, &self.ctx.pow(&self.root, index.into())))
    }

    pub fn elements(&self) -> Vec<GrElem> {
        let mut values = Vec::with_capacity(self.size as usize);
        let mut current = self.offset.clone();
        for _ in 0..self.size {
            values.push(current.clone());
            current = self.ctx.mul(&current, &self.root);
        }
        values
    }

    pub fn contains(&self, value: &GrElem) -> bool {
        let mut current = self.offset.clone();
        for _ in 0..self.size {
            if current == *value {
                return true;
            }
            current = self.ctx.mul(&current, &self.root);
        }
        false
    }

    pub fn is_teichmuller_subset(&self) -> bool {
        is_teichmuller_element(&self.ctx, &self.offset)
            && is_teichmuller_element(&self.ctx, &self.root)
    }

    pub fn scale(&self, power_factor: u64) -> Result<Self> {
        if power_factor == 0 || !self.size.is_multiple_of(power_factor) {
            return Err(GrError::InvalidDomain("scale requires power dividing size"));
        }

        Self::new(
            Arc::clone(&self.ctx),
            self.offset.clone(),
            self.ctx.pow(&self.root, power_factor.into()),
            self.size / power_factor,
        )
    }

    pub fn scale_offset(&self, power_factor: u64) -> Result<Self> {
        if power_factor == 0 || !self.size.is_multiple_of(power_factor) {
            return Err(GrError::InvalidDomain(
                "scale_offset requires power dividing size",
            ));
        }

        Self::new(
            Arc::clone(&self.ctx),
            self.ctx.mul(&self.offset, &self.root),
            self.ctx.pow(&self.root, power_factor.into()),
            self.size / power_factor,
        )
    }

    pub fn pow_map(&self, exponent: u64) -> Result<Self> {
        if exponent == 0 {
            return Err(GrError::InvalidDomain("pow_map requires exponent > 0"));
        }

        let common = gcd(self.size, exponent);
        Self::new(
            Arc::clone(&self.ctx),
            self.ctx.pow(&self.offset, exponent.into()),
            self.ctx.pow(&self.root, exponent.into()),
            self.size / common,
        )
    }

    pub fn disjoint_with(&self, other: &Self) -> Result<bool> {
        if self.ctx.config() != other.ctx.config() {
            return Err(GrError::DifferentRings);
        }

        for lhs in self.elements() {
            for rhs in other.elements() {
                if lhs == rhs {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }
}

const fn gcd(mut lhs: u64, mut rhs: u64) -> u64 {
    while rhs != 0 {
        let remainder = lhs % rhs;
        lhs = rhs;
        rhs = remainder;
    }
    lhs
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::Domain;
    use crate::algebra::galois_ring::{
        is_teichmuller_element, teichmuller_generator, teichmuller_subgroup_size_supported,
        GrConfig, GrContext, GrError,
    };

    fn ctx_r6() -> Arc<GrContext> {
        Arc::new(
            GrContext::new(GrConfig {
                p: 2,
                k_exp: 16,
                r: 6,
            })
            .unwrap(),
        )
    }

    #[test]
    fn domain_should_support_element_access_and_derived_subdomains() {
        let ctx = ctx_r6();
        let domain = Domain::teichmuller_subgroup(Arc::clone(&ctx), 9).unwrap();
        let offset = teichmuller_generator(&ctx).unwrap();
        let coset = Domain::teichmuller_coset(Arc::clone(&ctx), offset.clone(), 9).unwrap();

        assert_eq!(domain.size(), 9);
        assert_eq!(domain.element(0).unwrap(), ctx.one());
        assert_eq!(domain.scale(3).unwrap().size(), 3);
        assert_eq!(domain.scale_offset(3).unwrap().size(), 3);
        assert_eq!(domain.pow_map(3).unwrap().size(), 3);
        assert!(domain
            .pow_map(3)
            .unwrap()
            .disjoint_with(&domain.scale_offset(3).unwrap())
            .unwrap());
        assert_eq!(coset.size(), domain.size());
        assert_eq!(coset.root(), domain.root());
        assert_eq!(coset.offset(), &offset);
    }

    #[test]
    fn domain_root_order_and_pow_map_should_match() {
        let ctx = ctx_r6();
        let domain = Domain::teichmuller_subgroup(Arc::clone(&ctx), 9).unwrap();

        assert_eq!(ctx.pow(domain.root(), 9), ctx.one());
        assert_ne!(domain.root(), &ctx.one());
        assert_eq!(ctx.pow(domain.root(), 3), *domain.scale(3).unwrap().root());
        assert_eq!(domain.pow_map(9).unwrap().element(0).unwrap(), ctx.one());
    }

    #[test]
    fn domain_contains_should_accept_only_domain_points() {
        let ctx = ctx_r6();
        let subgroup = Domain::teichmuller_subgroup(Arc::clone(&ctx), 9).unwrap();
        let coset =
            Domain::teichmuller_coset(Arc::clone(&ctx), teichmuller_generator(&ctx).unwrap(), 9)
                .unwrap();

        assert!(subgroup.contains(&ctx.one()));
        assert!(subgroup.contains(&subgroup.element(7).unwrap()));
        assert!(!subgroup.contains(coset.offset()));
        assert!(!subgroup.contains(&ctx.zero()));
        assert!(!subgroup.contains(&ctx.from_u64(2)));
    }

    #[test]
    fn domain_teichmuller_subset_should_reject_generic_unit_coset() {
        let ctx = ctx_r6();
        let subgroup = Domain::teichmuller_subgroup(Arc::clone(&ctx), 9).unwrap();
        let teich_coset =
            Domain::teichmuller_coset(Arc::clone(&ctx), teichmuller_generator(&ctx).unwrap(), 9)
                .unwrap();
        let non_teich_unit = ctx.from_u64(3);
        let non_teich_coset =
            Domain::teichmuller_coset(Arc::clone(&ctx), non_teich_unit, 9).unwrap();

        assert!(subgroup.is_teichmuller_subset());
        assert!(teich_coset.is_teichmuller_subset());
        assert!(!non_teich_coset.is_teichmuller_subset());
    }

    #[test]
    fn invalid_domain_inputs_should_reject() {
        let ctx = ctx_r6();

        assert!(matches!(
            Domain::teichmuller_subgroup(Arc::clone(&ctx), 10),
            Err(GrError::InvalidSubgroupSize { .. })
        ));
        assert!(matches!(
            Domain::teichmuller_coset(Arc::clone(&ctx), ctx.zero(), 9),
            Err(GrError::InvalidDomain(_))
        ));
    }

    #[test]
    fn mid_extension_subgroups_should_be_supported() {
        let ctx = Arc::new(
            GrContext::new(GrConfig {
                p: 2,
                k_exp: 16,
                r: 18,
            })
            .unwrap(),
        );

        assert!(teichmuller_subgroup_size_supported(&ctx, 27));
        let subgroup = Domain::teichmuller_subgroup(Arc::clone(&ctx), 27).unwrap();
        assert_eq!(subgroup.size(), 27);
        assert!(is_teichmuller_element(&ctx, subgroup.root()));
    }
}
