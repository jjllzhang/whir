use crate::{
    algebra::galois_ring::{GrContext, GrElem, GrError, Result},
    hash::Hash,
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ByteMerkleTree {
    leaf_count: usize,
    layers: Vec<Vec<Hash>>,
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
    pub fn commit(leaf_payloads: Vec<Vec<u8>>) -> Result<Self> {
        if leaf_payloads.is_empty() {
            return Err(GrError::InvalidDomain(
                "Merkle tree requires at least one leaf",
            ));
        }

        let leaf_count = leaf_payloads.len();
        let padded_leaf_count = leaf_count.next_power_of_two();
        let mut leaf_hashes = leaf_payloads
            .iter()
            .map(|payload| hash_leaf(payload))
            .collect::<Vec<_>>();
        leaf_hashes.resize(padded_leaf_count, Hash::default());

        let mut layers = vec![leaf_hashes];
        while layers.last().expect("layers is nonempty").len() > 1 {
            let previous = layers.last().expect("layers is nonempty");
            let mut next = Vec::with_capacity(previous.len() / 2);
            for pair in previous.chunks_exact(2) {
                next.push(hash_node(&pair[0], &pair[1]));
            }
            layers.push(next);
        }

        Ok(Self {
            leaf_count,
            layers,
            leaf_payloads,
        })
    }

    pub fn root(&self) -> Hash {
        self.layers
            .last()
            .and_then(|layer| layer.first())
            .copied()
            .unwrap_or_default()
    }

    pub const fn leaf_count(&self) -> usize {
        self.leaf_count
    }

    pub fn open(&self, indices: &[usize]) -> Result<MerkleProof> {
        if indices.iter().any(|&index| index >= self.leaf_count) {
            return Err(GrError::IndexOutOfRange {
                index: indices.iter().copied().max().unwrap_or(0) as u64,
                size: self.leaf_count as u64,
            });
        }

        let depth = self.layers.len().saturating_sub(1);
        let mut sibling_hashes = Vec::with_capacity(indices.len() * depth);
        for &index in indices {
            let mut current_index = index;
            for layer in self.layers.iter().take(depth) {
                sibling_hashes.push(layer[current_index ^ 1]);
                current_index >>= 1;
            }
        }

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

    pub fn verify(root: Hash, proof: &MerkleProof) -> Result<bool> {
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

        let mut siblings = proof.sibling_hashes.iter();
        for (&index, payload) in proof.queried_indices.iter().zip(&proof.leaf_payloads) {
            let mut current_index = index;
            let mut current = hash_leaf(payload);
            for sibling in siblings.by_ref().take(depth) {
                current = if current_index & 1 == 0 {
                    hash_node(&current, sibling)
                } else {
                    hash_node(sibling, &current)
                };
                current_index >>= 1;
            }
            if current != root {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

pub fn build_oracle_leaves(ctx: &GrContext, oracle_evals: &[GrElem]) -> Vec<Vec<u8>> {
    oracle_evals
        .iter()
        .map(|value| ctx.serialize(value))
        .collect()
}

pub fn build_oracle_tree(ctx: &GrContext, oracle_evals: &[GrElem]) -> Result<ByteMerkleTree> {
    ByteMerkleTree::commit(build_oracle_leaves(ctx, oracle_evals))
}

fn hash_leaf(payload: &[u8]) -> Hash {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"whir-gr.merkle.leaf.v1");
    hasher.update(&(payload.len() as u64).to_le_bytes());
    hasher.update(payload);
    Hash(hasher.finalize().into())
}

fn hash_node(lhs: &Hash, rhs: &Hash) -> Hash {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"whir-gr.merkle.node.v1");
    hasher.update(&lhs.0);
    hasher.update(&rhs.0);
    Hash(hasher.finalize().into())
}

#[cfg(test)]
mod tests {
    use super::{build_oracle_tree, ByteMerkleTree};
    use crate::algebra::galois_ring::{GrConfig, GrContext};

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

        let lhs = build_oracle_tree(&ctx, &values).unwrap();
        let rhs = build_oracle_tree(&ctx, &values).unwrap();

        assert_eq!(lhs.root(), rhs.root());
    }

    #[test]
    fn merkle_opening_should_verify() {
        let ctx = sample_context();
        let values = vec![ctx.one(), ctx.from_u64(2), ctx.from_u64(3), ctx.from_u64(4)];
        let tree = build_oracle_tree(&ctx, &values).unwrap();

        let proof = tree.open(&[1, 3]).unwrap();

        assert!(ByteMerkleTree::verify(tree.root(), &proof).unwrap());
    }

    #[test]
    fn merkle_verification_should_reject_tampered_payload() {
        let ctx = sample_context();
        let values = vec![ctx.one(), ctx.from_u64(2), ctx.from_u64(3), ctx.from_u64(4)];
        let tree = build_oracle_tree(&ctx, &values).unwrap();
        let mut proof = tree.open(&[1]).unwrap();
        proof.leaf_payloads[0][0] ^= 1;

        assert!(!ByteMerkleTree::verify(tree.root(), &proof).unwrap());
    }
}
