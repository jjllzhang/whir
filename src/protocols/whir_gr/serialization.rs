use crate::{
    algebra::galois_ring::{Domain, GrContext, GrElem},
    hash::Hash,
    protocols::whir_gr::{
        common::{
            WhirGrOpening, WhirGrProof, WhirGrPublicParameters, WhirGrRoundProof,
            WhirGrSumcheckPolynomial,
        },
        merkle::MerkleProof,
    },
};

#[derive(Clone, Debug, Default)]
pub struct ByteWriter {
    bytes: Vec<u8>,
}

impl ByteWriter {
    pub const fn new() -> Self {
        Self { bytes: Vec::new() }
    }

    pub fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }

    pub fn write_u64(&mut self, value: u64) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    pub fn write_bytes(&mut self, bytes: &[u8]) {
        self.write_u64(bytes.len() as u64);
        self.bytes.extend_from_slice(bytes);
    }

    pub fn write_raw_bytes(&mut self, bytes: &[u8]) {
        self.bytes.extend_from_slice(bytes);
    }

    pub fn write_hash(&mut self, hash: &Hash) {
        self.bytes.extend_from_slice(&hash.0);
    }

    pub fn write_ring_element(&mut self, ctx: &GrContext, value: &GrElem) {
        self.write_raw_bytes(&ctx.serialize(value));
    }

    pub fn write_ring_vector(&mut self, ctx: &GrContext, values: &[GrElem]) {
        self.write_u64(values.len() as u64);
        for value in values {
            self.write_ring_element(ctx, value);
        }
    }

    pub fn write_u64_vector(&mut self, values: &[u64]) {
        self.write_u64(values.len() as u64);
        for &value in values {
            self.write_u64(value);
        }
    }

    pub fn write_byte_vectors(&mut self, values: &[Vec<u8>]) {
        self.write_u64(values.len() as u64);
        for value in values {
            self.write_bytes(value);
        }
    }
}

pub fn serialize_ring_element(ctx: &GrContext, value: &GrElem) -> Vec<u8> {
    ctx.serialize(value)
}

pub fn serialize_ring_vector(ctx: &GrContext, values: &[GrElem]) -> Vec<u8> {
    let mut writer = ByteWriter::new();
    writer.write_ring_vector(ctx, values);
    writer.into_bytes()
}

pub fn serialize_domain(domain: &Domain) -> Vec<u8> {
    let ctx = domain.context();
    let mut writer = ByteWriter::new();
    writer.write_u64(ctx.config().p);
    writer.write_u64(u64::from(ctx.config().k_exp));
    writer.write_u64(ctx.config().r as u64);
    writer.write_u64(domain.size());
    writer.write_ring_element(ctx, domain.offset());
    writer.write_ring_element(ctx, domain.root());
    writer.into_bytes()
}

pub fn serialize_public_parameters(params: &WhirGrPublicParameters) -> Vec<u8> {
    let ctx = &params.ctx;
    let mut writer = ByteWriter::new();
    writer.write_u64(ctx.config().p);
    writer.write_u64(u64::from(ctx.config().k_exp));
    writer.write_u64(ctx.config().r as u64);
    writer.write_u64(params.variable_count);
    writer.write_bytes(&serialize_domain(&params.initial_domain));
    writer.write_ring_element(ctx, &params.omega);
    writer.write_u64(params.ternary_grid.len() as u64);
    for point in &params.ternary_grid {
        writer.write_ring_element(ctx, point);
    }
    writer.write_u64(params.lambda_target);
    writer.write_raw_bytes(params.hash_id.as_slice());
    writer.write_u64_vector(&params.layer_widths);
    writer.write_u64_vector(&params.shift_repetitions);
    writer.write_u64(params.final_repetitions);
    writer.write_u64_vector(&params.degree_bounds);
    writer.into_bytes()
}

pub fn serialize_sumcheck_polynomial(
    ctx: &GrContext,
    polynomial: &WhirGrSumcheckPolynomial,
) -> Vec<u8> {
    serialize_ring_vector(ctx, &polynomial.coefficients)
}

pub fn serialize_merkle_proof(proof: &MerkleProof) -> Vec<u8> {
    let mut writer = ByteWriter::new();
    writer.write_u64(proof.leaf_count as u64);
    writer.write_u64_vector(
        &proof
            .queried_indices
            .iter()
            .map(|&index| index as u64)
            .collect::<Vec<_>>(),
    );
    writer.write_byte_vectors(&proof.leaf_payloads);
    writer.write_u64(proof.sibling_hashes.len() as u64);
    for hash in &proof.sibling_hashes {
        writer.write_hash(hash);
    }
    writer.into_bytes()
}

pub fn serialize_round_proof(ctx: &GrContext, proof: &WhirGrRoundProof) -> Vec<u8> {
    let mut writer = ByteWriter::new();
    writer.write_u64(proof.sumcheck_polynomials.len() as u64);
    for polynomial in &proof.sumcheck_polynomials {
        writer.write_bytes(&serialize_sumcheck_polynomial(ctx, polynomial));
    }
    writer.write_hash(&proof.g_root);
    writer.write_bytes(&serialize_merkle_proof(&proof.virtual_fold_openings));
    writer.into_bytes()
}

pub fn serialize_proof(ctx: &GrContext, proof: &WhirGrProof) -> Vec<u8> {
    let mut writer = ByteWriter::new();
    writer.write_u64(proof.rounds.len() as u64);
    for round in &proof.rounds {
        writer.write_bytes(&serialize_round_proof(ctx, round));
    }
    writer.write_ring_element(ctx, &proof.final_constant);
    writer.write_bytes(&serialize_merkle_proof(&proof.final_openings));
    writer.into_bytes()
}

pub fn serialize_opening(ctx: &GrContext, opening: &WhirGrOpening) -> Vec<u8> {
    let mut writer = ByteWriter::new();
    writer.write_ring_element(ctx, &opening.value);
    writer.write_bytes(&serialize_proof(ctx, &opening.proof));
    writer.into_bytes()
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::{
        algebra::galois_ring::{teichmuller_subgroup_generator, Domain, GrConfig, GrContext},
        protocols::whir_gr::{
            common::WhirGrPublicParameters,
            serialization::{serialize_domain, serialize_public_parameters, serialize_ring_vector},
        },
    };

    fn sample_context() -> Arc<GrContext> {
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
    fn ring_vector_serialization_should_be_length_prefixed() {
        let ctx = sample_context();
        let values = vec![ctx.one(), ctx.from_u64(7)];

        let encoded = serialize_ring_vector(&ctx, &values);

        assert_eq!(encoded.len(), 8 + values.len() * ctx.elem_bytes());
    }

    #[test]
    fn domain_serialization_should_change_when_offset_changes() {
        let ctx = sample_context();
        let subgroup = Domain::teichmuller_subgroup(Arc::clone(&ctx), 9).unwrap();
        let coset = Domain::teichmuller_coset(
            Arc::clone(&ctx),
            teichmuller_subgroup_generator(&ctx, 7).unwrap(),
            9,
        )
        .unwrap();

        assert_ne!(serialize_domain(&subgroup), serialize_domain(&coset));
    }

    #[test]
    fn public_parameters_serialization_should_be_deterministic() {
        let ctx = sample_context();
        let domain = Domain::teichmuller_subgroup(Arc::clone(&ctx), 9).unwrap();
        let omega = domain.root().clone();
        let params = WhirGrPublicParameters::new(
            Arc::clone(&ctx),
            domain,
            2,
            omega.clone(),
            [ctx.one(), omega.clone(), ctx.square(&omega)],
        );

        assert_eq!(
            serialize_public_parameters(&params),
            serialize_public_parameters(&params)
        );
    }
}
