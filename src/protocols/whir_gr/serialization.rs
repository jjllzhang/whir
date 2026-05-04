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
    transcript::{Encoding, NargDeserialize, VerificationError, VerificationResult},
};

#[derive(Clone, Debug, Default)]
pub struct ByteWriter {
    bytes: Vec<u8>,
}

#[derive(Clone, Debug)]
pub struct ByteReader<'a> {
    bytes: &'a [u8],
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WhirGrOpeningPayload {
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

impl<'a> ByteReader<'a> {
    pub const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes }
    }

    pub const fn remaining(&self) -> &[u8] {
        self.bytes
    }

    pub const fn is_empty(&self) -> bool {
        self.bytes.is_empty()
    }

    pub fn read_u64(&mut self) -> crate::algebra::galois_ring::Result<u64> {
        if self.bytes.len() < 8 {
            return Err(crate::algebra::galois_ring::GrError::InvalidDomain(
                "WHIR_GR byte stream ended while reading u64",
            ));
        }
        let (head, tail) = self.bytes.split_at(8);
        self.bytes = tail;
        let value = u64::from_le_bytes(head.try_into().map_err(|_| {
            crate::algebra::galois_ring::GrError::InvalidDomain(
                "WHIR_GR byte stream contains malformed u64",
            )
        })?);
        Ok(value)
    }

    pub fn read_bytes(&mut self) -> crate::algebra::galois_ring::Result<Vec<u8>> {
        let len = self.read_u64()?;
        let len = usize::try_from(len).map_err(|_| {
            crate::algebra::galois_ring::GrError::ArithmeticOverflow("byte vector length")
        })?;
        self.read_raw_bytes(len)
    }

    pub fn read_raw_bytes(&mut self, len: usize) -> crate::algebra::galois_ring::Result<Vec<u8>> {
        if self.bytes.len() < len {
            return Err(crate::algebra::galois_ring::GrError::InvalidDomain(
                "WHIR_GR byte stream ended while reading bytes",
            ));
        }
        let (head, tail) = self.bytes.split_at(len);
        self.bytes = tail;
        Ok(head.to_vec())
    }

    pub fn read_hash(&mut self) -> crate::algebra::galois_ring::Result<Hash> {
        let bytes = self.read_raw_bytes(32)?;
        Ok(Hash(bytes.try_into().map_err(|_| {
            crate::algebra::galois_ring::GrError::InvalidDomain("malformed WHIR_GR hash")
        })?))
    }

    pub fn read_ring_element(
        &mut self,
        ctx: &GrContext,
    ) -> crate::algebra::galois_ring::Result<GrElem> {
        let bytes = self.read_raw_bytes(ctx.elem_bytes())?;
        ctx.deserialize(&bytes)
    }

    pub fn read_ring_vector(
        &mut self,
        ctx: &GrContext,
    ) -> crate::algebra::galois_ring::Result<Vec<GrElem>> {
        let len = self.read_u64()?;
        let len = usize::try_from(len).map_err(|_| {
            crate::algebra::galois_ring::GrError::ArithmeticOverflow("ring vector length")
        })?;
        (0..len).map(|_| self.read_ring_element(ctx)).collect()
    }

    pub fn read_u64_vector(&mut self) -> crate::algebra::galois_ring::Result<Vec<u64>> {
        let len = self.read_u64()?;
        let len = usize::try_from(len).map_err(|_| {
            crate::algebra::galois_ring::GrError::ArithmeticOverflow("u64 vector length")
        })?;
        (0..len).map(|_| self.read_u64()).collect()
    }

    pub fn read_byte_vectors(&mut self) -> crate::algebra::galois_ring::Result<Vec<Vec<u8>>> {
        let len = self.read_u64()?;
        let len = usize::try_from(len).map_err(|_| {
            crate::algebra::galois_ring::GrError::ArithmeticOverflow("byte vector count")
        })?;
        (0..len).map(|_| self.read_bytes()).collect()
    }
}

impl WhirGrOpeningPayload {
    pub const fn new(bytes: Vec<u8>) -> Self {
        Self { bytes }
    }

    pub fn from_opening(ctx: &GrContext, opening: &WhirGrOpening) -> Self {
        Self::new(serialize_opening(ctx, opening))
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.bytes
    }

    pub fn into_opening(
        self,
        ctx: &GrContext,
    ) -> crate::algebra::galois_ring::Result<WhirGrOpening> {
        deserialize_opening(ctx, &self.bytes)
    }
}

impl Encoding<[u8]> for WhirGrOpeningPayload {
    fn encode(&self) -> impl AsRef<[u8]> {
        let mut encoded = Vec::with_capacity(8 + self.bytes.len());
        encoded.extend_from_slice(&(self.bytes.len() as u64).to_le_bytes());
        encoded.extend_from_slice(&self.bytes);
        encoded
    }
}

impl NargDeserialize for WhirGrOpeningPayload {
    fn deserialize_from_narg(buf: &mut &[u8]) -> VerificationResult<Self> {
        if buf.len() < 8 {
            return Err(VerificationError);
        }
        let (len_bytes, tail) = buf.split_at(8);
        let len = u64::from_le_bytes(len_bytes.try_into().map_err(|_| VerificationError)?);
        let len = usize::try_from(len).map_err(|_| VerificationError)?;
        if tail.len() < len {
            return Err(VerificationError);
        }
        let (payload, remaining) = tail.split_at(len);
        *buf = remaining;
        Ok(Self::new(payload.to_vec()))
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

pub fn deserialize_ring_element(
    ctx: &GrContext,
    bytes: &[u8],
) -> crate::algebra::galois_ring::Result<GrElem> {
    ctx.deserialize(bytes)
}

pub fn deserialize_ring_vector(
    ctx: &GrContext,
    bytes: &[u8],
) -> crate::algebra::galois_ring::Result<Vec<GrElem>> {
    let mut reader = ByteReader::new(bytes);
    let values = reader.read_ring_vector(ctx)?;
    ensure_eof(&reader)?;
    Ok(values)
}

pub fn deserialize_sumcheck_polynomial(
    ctx: &GrContext,
    bytes: &[u8],
) -> crate::algebra::galois_ring::Result<WhirGrSumcheckPolynomial> {
    Ok(WhirGrSumcheckPolynomial {
        coefficients: deserialize_ring_vector(ctx, bytes)?,
    })
}

pub fn deserialize_merkle_proof(bytes: &[u8]) -> crate::algebra::galois_ring::Result<MerkleProof> {
    let mut reader = ByteReader::new(bytes);
    let proof = read_merkle_proof(&mut reader)?;
    ensure_eof(&reader)?;
    Ok(proof)
}

pub fn deserialize_round_proof(
    ctx: &GrContext,
    bytes: &[u8],
) -> crate::algebra::galois_ring::Result<WhirGrRoundProof> {
    let mut reader = ByteReader::new(bytes);
    let proof = read_round_proof(ctx, &mut reader)?;
    ensure_eof(&reader)?;
    Ok(proof)
}

pub fn deserialize_proof(
    ctx: &GrContext,
    bytes: &[u8],
) -> crate::algebra::galois_ring::Result<WhirGrProof> {
    let mut reader = ByteReader::new(bytes);
    let proof = read_proof(ctx, &mut reader)?;
    ensure_eof(&reader)?;
    Ok(proof)
}

pub fn deserialize_opening(
    ctx: &GrContext,
    bytes: &[u8],
) -> crate::algebra::galois_ring::Result<WhirGrOpening> {
    let mut reader = ByteReader::new(bytes);
    let value = reader.read_ring_element(ctx)?;
    let proof_bytes = reader.read_bytes()?;
    let proof = deserialize_proof(ctx, &proof_bytes)?;
    ensure_eof(&reader)?;
    Ok(WhirGrOpening { value, proof })
}

fn read_merkle_proof(
    reader: &mut ByteReader<'_>,
) -> crate::algebra::galois_ring::Result<MerkleProof> {
    let leaf_count = reader.read_u64()?;
    let leaf_count = usize::try_from(leaf_count).map_err(|_| {
        crate::algebra::galois_ring::GrError::ArithmeticOverflow("Merkle leaf count")
    })?;
    let queried_indices = reader
        .read_u64_vector()?
        .into_iter()
        .map(|index| {
            usize::try_from(index).map_err(|_| {
                crate::algebra::galois_ring::GrError::ArithmeticOverflow("Merkle query index")
            })
        })
        .collect::<crate::algebra::galois_ring::Result<Vec<_>>>()?;
    let leaf_payloads = reader.read_byte_vectors()?;
    let sibling_count = reader.read_u64()?;
    let sibling_count = usize::try_from(sibling_count).map_err(|_| {
        crate::algebra::galois_ring::GrError::ArithmeticOverflow("Merkle sibling count")
    })?;
    let sibling_hashes = (0..sibling_count)
        .map(|_| reader.read_hash())
        .collect::<crate::algebra::galois_ring::Result<Vec<_>>>()?;
    Ok(MerkleProof {
        leaf_count,
        queried_indices,
        leaf_payloads,
        sibling_hashes,
    })
}

fn read_round_proof(
    ctx: &GrContext,
    reader: &mut ByteReader<'_>,
) -> crate::algebra::galois_ring::Result<WhirGrRoundProof> {
    let polynomial_count = reader.read_u64()?;
    let polynomial_count = usize::try_from(polynomial_count).map_err(|_| {
        crate::algebra::galois_ring::GrError::ArithmeticOverflow("sumcheck polynomial count")
    })?;
    let mut sumcheck_polynomials = Vec::with_capacity(polynomial_count);
    for _ in 0..polynomial_count {
        let polynomial_bytes = reader.read_bytes()?;
        sumcheck_polynomials.push(deserialize_sumcheck_polynomial(ctx, &polynomial_bytes)?);
    }
    let g_root = reader.read_hash()?;
    let virtual_fold_bytes = reader.read_bytes()?;
    let virtual_fold_openings = deserialize_merkle_proof(&virtual_fold_bytes)?;
    Ok(WhirGrRoundProof {
        sumcheck_polynomials,
        g_root,
        virtual_fold_openings,
    })
}

fn read_proof(
    ctx: &GrContext,
    reader: &mut ByteReader<'_>,
) -> crate::algebra::galois_ring::Result<WhirGrProof> {
    let round_count = reader.read_u64()?;
    let round_count = usize::try_from(round_count)
        .map_err(|_| crate::algebra::galois_ring::GrError::ArithmeticOverflow("round count"))?;
    let mut rounds = Vec::with_capacity(round_count);
    for _ in 0..round_count {
        let round_bytes = reader.read_bytes()?;
        rounds.push(deserialize_round_proof(ctx, &round_bytes)?);
    }
    let final_constant = reader.read_ring_element(ctx)?;
    let final_opening_bytes = reader.read_bytes()?;
    let final_openings = deserialize_merkle_proof(&final_opening_bytes)?;
    Ok(WhirGrProof {
        rounds,
        final_constant,
        final_openings,
    })
}

const fn ensure_eof(reader: &ByteReader<'_>) -> crate::algebra::galois_ring::Result<()> {
    if reader.is_empty() {
        Ok(())
    } else {
        Err(crate::algebra::galois_ring::GrError::InvalidDomain(
            "WHIR_GR byte stream has trailing bytes",
        ))
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::{
        algebra::galois_ring::{teichmuller_subgroup_generator, Domain, GrConfig, GrContext},
        protocols::whir_gr::{
            common::WhirGrPublicParameters,
            merkle::build_oracle_tree,
            serialization::{
                deserialize_merkle_proof, deserialize_ring_vector, serialize_domain,
                serialize_merkle_proof, serialize_public_parameters, serialize_ring_vector,
                WhirGrOpeningPayload,
            },
        },
        transcript::{NargDeserialize, NargSerialize},
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
        assert_eq!(deserialize_ring_vector(&ctx, &encoded).unwrap(), values);
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

    #[test]
    fn merkle_proof_serialization_should_roundtrip() {
        let ctx = sample_context();
        let values = vec![ctx.one(), ctx.from_u64(2), ctx.from_u64(3), ctx.from_u64(4)];
        let tree = build_oracle_tree(crate::hash::BLAKE3, &ctx, &values).unwrap();
        let proof = tree.open(&[1, 3]).unwrap();

        let encoded = serialize_merkle_proof(&proof);
        let decoded = deserialize_merkle_proof(&encoded).unwrap();

        assert_eq!(decoded, proof);
    }

    #[test]
    fn opening_payload_should_roundtrip_through_narg() {
        let payload = WhirGrOpeningPayload::new(vec![1, 2, 3, 4]);
        let mut encoded = Vec::new();
        payload.serialize_into_narg(&mut encoded);

        let mut read = encoded.as_slice();
        let decoded = WhirGrOpeningPayload::deserialize_from_narg(&mut read).unwrap();

        assert_eq!(decoded.as_bytes(), payload.as_bytes());
        assert!(read.is_empty());
    }
}
