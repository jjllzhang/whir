use std::sync::Arc;

use crate::{
    algebra::galois_ring::{Domain, GrContext, GrElem},
    engines::EngineId,
    hash::{Hash, BLAKE3},
};

#[derive(Clone, Debug)]
pub struct WhirGrPublicParameters {
    pub ctx: Arc<GrContext>,
    pub initial_domain: Domain,
    pub variable_count: u64,
    pub layer_widths: Vec<u64>,
    pub shift_repetitions: Vec<u64>,
    pub final_repetitions: u64,
    pub degree_bounds: Vec<u64>,
    pub omega: GrElem,
    pub ternary_grid: [GrElem; 3],
    pub lambda_target: u64,
    pub hash_id: EngineId,
}

impl WhirGrPublicParameters {
    pub const fn new(
        ctx: Arc<GrContext>,
        initial_domain: Domain,
        variable_count: u64,
        omega: GrElem,
        ternary_grid: [GrElem; 3],
    ) -> Self {
        Self {
            ctx,
            initial_domain,
            variable_count,
            layer_widths: Vec::new(),
            shift_repetitions: Vec::new(),
            final_repetitions: 0,
            degree_bounds: Vec::new(),
            omega,
            ternary_grid,
            lambda_target: 128,
            hash_id: BLAKE3,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WhirGrCommitment {
    pub oracle_root: Hash,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WhirGrSumcheckPolynomial {
    pub coefficients: Vec<GrElem>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WhirGrRoundProof {
    pub sumcheck_polynomials: Vec<WhirGrSumcheckPolynomial>,
    pub g_root: Hash,
    pub virtual_fold_openings: crate::protocols::whir_gr::merkle::MerkleProof,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WhirGrProof {
    pub rounds: Vec<WhirGrRoundProof>,
    pub final_constant: GrElem,
    pub final_openings: crate::protocols::whir_gr::merkle::MerkleProof,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WhirGrOpening {
    pub value: GrElem,
    pub proof: WhirGrProof,
}
