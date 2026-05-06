use core::fmt;

use serde::{Deserialize, Serialize};

#[derive(Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GrElem {
    coefficients: Vec<u64>,
}

impl GrElem {
    pub(crate) const fn new_unchecked(coefficients: Vec<u64>) -> Self {
        Self { coefficients }
    }

    pub fn coefficients(&self) -> &[u64] {
        &self.coefficients
    }

    pub(crate) fn coefficients_mut(&mut self) -> &mut [u64] {
        &mut self.coefficients
    }

    pub fn into_coefficients(self) -> Vec<u64> {
        self.coefficients
    }

    pub fn is_zero(&self) -> bool {
        self.coefficients
            .iter()
            .all(|&coefficient| coefficient == 0)
    }
}

impl fmt::Debug for GrElem {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_tuple("GrElem")
            .field(&self.coefficients)
            .finish()
    }
}
