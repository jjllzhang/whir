use crate::{
    algebra::galois_ring::{GrContext, GrElem, GrError, Result},
    engines::EngineId,
    hash::{self, Hash, ENGINES},
    protocols::merkle_tree,
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ByteMerkleTree {
    hash_id: EngineId,
    leaf_count: usize,
    config: merkle_tree::Config,
    commitment: merkle_tree::Commitment,
    witness: merkle_tree::Witness,
    leaf_payloads: Vec<Vec<u8>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MerkleProof {
    pub leaf_count: usize,
    pub queried_indices: Vec<usize>,
    pub leaf_payloads: Vec<Vec<u8>>,
    pub sibling_hashes: Vec<Hash>,
}

impl ByteMerkleTree {
    pub fn commit(hash_id: EngineId, leaf_payloads: Vec<Vec<u8>>) -> Result<Self> {
        if leaf_payloads.is_empty() {
            return Err(GrError::InvalidDomain(
                "Merkle tree requires at least one leaf",
            ));
        }

        let leaf_count = leaf_payloads.len();
        let leaf_hashes = leaf_payloads
            .iter()
            .map(|payload| hash_leaf(hash_id, payload))
            .collect::<Result<Vec<_>>>()?;
        let config = merkle_tree::Config::with_hash(hash_id, leaf_count);
        let (commitment, witness) = config.build(leaf_hashes);

        Ok(Self {
            hash_id,
            leaf_count,
            config,
            commitment,
            witness,
            leaf_payloads,
        })
    }

    pub const fn root(&self) -> Hash {
        self.commitment.hash()
    }

    pub const fn leaf_count(&self) -> usize {
        self.leaf_count
    }

    pub const fn hash_id(&self) -> EngineId {
        self.hash_id
    }

    pub fn open(&self, indices: &[usize]) -> Result<MerkleProof> {
        if indices.iter().any(|&index| index >= self.leaf_count) {
            return Err(GrError::IndexOutOfRange {
                index: indices.iter().copied().max().unwrap_or(0) as u64,
                size: self.leaf_count as u64,
            });
        }

        let sibling_hashes = self
            .config
            .open_paths(&self.witness, indices)
            .map_err(|_| GrError::InvalidDomain("failed to open WHIR_GR Merkle paths"))?;

        Ok(MerkleProof {
            leaf_count: self.leaf_count,
            queried_indices: indices.to_vec(),
            leaf_payloads: indices
                .iter()
                .map(|&index| self.leaf_payloads[index].clone())
                .collect(),
            sibling_hashes,
        })
    }

    pub fn verify(hash_id: EngineId, root: Hash, proof: &MerkleProof) -> Result<bool> {
        if proof.leaf_count == 0 {
            return Err(GrError::InvalidDomain("Merkle proof leaf count is zero"));
        }
        if proof.queried_indices.len() != proof.leaf_payloads.len() {
            return Err(GrError::InvalidDomain(
                "Merkle proof index/payload length mismatch",
            ));
        }
        if proof
            .queried_indices
            .iter()
            .any(|&index| index >= proof.leaf_count)
        {
            return Err(GrError::IndexOutOfRange {
                index: proof.queried_indices.iter().copied().max().unwrap_or(0) as u64,
                size: proof.leaf_count as u64,
            });
        }

        let depth = proof.leaf_count.next_power_of_two().ilog2() as usize;
        if proof.sibling_hashes.len() != proof.queried_indices.len() * depth {
            return Err(GrError::InvalidDomain(
                "Merkle proof sibling count mismatch",
            ));
        }

        if ENGINES.retrieve(hash_id).is_none() {
            return Err(GrError::InvalidDomain(
                "WHIR_GR hash engine is not registered",
            ));
        }
        let config = merkle_tree::Config::with_hash(hash_id, proof.leaf_count);
        let commitment = merkle_tree::Commitment::new(root);
        let leaf_hashes = proof
            .leaf_payloads
            .iter()
            .map(|payload| hash_leaf(hash_id, payload))
            .collect::<Result<Vec<_>>>()?;
        Ok(config
            .verify_paths(
                &commitment,
                &proof.queried_indices,
                &leaf_hashes,
                &proof.sibling_hashes,
            )
            .is_ok())
    }
}

pub fn build_oracle_leaves(ctx: &GrContext, oracle_evals: &[GrElem]) -> Vec<Vec<u8>> {
    oracle_evals
        .iter()
        .map(|value| ctx.serialize(value))
        .collect()
}

pub fn build_oracle_tree(
    hash_id: EngineId,
    ctx: &GrContext,
    oracle_evals: &[GrElem],
) -> Result<ByteMerkleTree> {
    ByteMerkleTree::commit(hash_id, build_oracle_leaves(ctx, oracle_evals))
}

fn hash_leaf(hash_id: EngineId, payload: &[u8]) -> Result<Hash> {
    hash_framed(
        hash_id,
        &[
            b"whir-gr.merkle.leaf.v1",
            &(payload.len() as u64).to_le_bytes(),
            payload,
        ],
    )
}

fn hash_framed(hash_id: EngineId, parts: &[&[u8]]) -> Result<Hash> {
    let engine = ENGINES.retrieve(hash_id).ok_or(GrError::InvalidDomain(
        "WHIR_GR hash engine is not registered",
    ))?;
    let input_len = parts
        .iter()
        .try_fold(0usize, |acc, part| acc.checked_add(part.len()))
        .ok_or(GrError::ArithmeticOverflow("WHIR_GR hash input length"))?;
    if hash_id == hash::COPY && input_len > 32 {
        return Err(GrError::InvalidDomain(
            "WHIR_GR Merkle hashing cannot use Copy for framed inputs larger than 32 bytes",
        ));
    }

    let mut input = Vec::with_capacity(input_len);
    for part in parts {
        input.extend_from_slice(part);
    }
    let mut out = Hash::default();
    engine.hash_many(input.len(), &input, std::slice::from_mut(&mut out));
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::{build_oracle_tree, ByteMerkleTree};
    use crate::{
        algebra::galois_ring::{GrConfig, GrContext},
        hash::{BLAKE3, SHA2},
    };

    fn sample_context() -> GrContext {
        GrContext::new(GrConfig {
            p: 2,
            k_exp: 16,
            r: 6,
        })
        .unwrap()
    }

    #[test]
    fn merkle_root_should_be_deterministic() {
        let ctx = sample_context();
        let values = vec![ctx.one(), ctx.from_u64(2), ctx.from_u64(3)];

        let lhs = build_oracle_tree(BLAKE3, &ctx, &values).unwrap();
        let rhs = build_oracle_tree(BLAKE3, &ctx, &values).unwrap();

        assert_eq!(lhs.root(), rhs.root());
        assert_eq!(lhs.hash_id(), BLAKE3);
    }

    #[test]
    fn merkle_root_should_depend_on_hash_engine() {
        let ctx = sample_context();
        let values = vec![ctx.one(), ctx.from_u64(2), ctx.from_u64(3)];

        let blake3_tree = build_oracle_tree(BLAKE3, &ctx, &values).unwrap();
        let sha2_tree = build_oracle_tree(SHA2, &ctx, &values).unwrap();

        assert_ne!(blake3_tree.root(), sha2_tree.root());
    }

    #[test]
    fn merkle_opening_should_verify() {
        let ctx = sample_context();
        let values = vec![ctx.one(), ctx.from_u64(2), ctx.from_u64(3), ctx.from_u64(4)];
        let tree = build_oracle_tree(BLAKE3, &ctx, &values).unwrap();

        let proof = tree.open(&[1, 3]).unwrap();

        assert!(ByteMerkleTree::verify(BLAKE3, tree.root(), &proof).unwrap());
        assert!(!ByteMerkleTree::verify(SHA2, tree.root(), &proof).unwrap());
    }

    #[test]
    fn merkle_verification_should_reject_tampered_payload() {
        let ctx = sample_context();
        let values = vec![ctx.one(), ctx.from_u64(2), ctx.from_u64(3), ctx.from_u64(4)];
        let tree = build_oracle_tree(BLAKE3, &ctx, &values).unwrap();
        let mut proof = tree.open(&[1]).unwrap();
        proof.leaf_payloads[0][0] ^= 1;

        assert!(!ByteMerkleTree::verify(BLAKE3, tree.root(), &proof).unwrap());
    }
}
