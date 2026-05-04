use blake3::Hasher;

use crate::algebra::galois_ring::{is_teichmuller_element, GrContext, GrElem, GrError, Result};

#[derive(Clone, Debug)]
pub struct Transcript {
    hasher: Hasher,
}

impl Transcript {
    pub fn new(domain: &[u8]) -> Self {
        let mut transcript = Self {
            hasher: Hasher::new(),
        };
        transcript.absorb_labeled_bytes(b"domain", domain);
        transcript
    }

    pub fn absorb_labeled_bytes(&mut self, label: &[u8], bytes: &[u8]) {
        absorb_labeled_bytes(&mut self.hasher, label, bytes);
    }

    pub fn absorb_ring_element(&mut self, ctx: &GrContext, label: &[u8], value: &GrElem) {
        self.absorb_labeled_bytes(label, &ctx.serialize(value));
    }

    pub fn absorb_ring_vector(&mut self, ctx: &GrContext, label: &[u8], values: &[GrElem]) {
        let mut bytes = Vec::with_capacity(8 + values.len() * ctx.elem_bytes());
        bytes.extend_from_slice(&(values.len() as u64).to_le_bytes());
        for value in values {
            bytes.extend_from_slice(&ctx.serialize(value));
        }
        self.absorb_labeled_bytes(label, &bytes);
    }

    pub fn challenge_bytes(&self, label: &[u8], out: &mut [u8]) {
        let mut hasher = self.hasher.clone();
        absorb_labeled_bytes(&mut hasher, b"challenge", label);
        hasher.finalize_xof().fill(out);
    }

    pub fn challenge_u64(&self, label: &[u8]) -> u64 {
        let mut bytes = [0; 8];
        self.challenge_bytes(label, &mut bytes);
        u64::from_le_bytes(bytes)
    }

    pub fn challenge_index(&self, label: &[u8], modulus: u64) -> Result<u64> {
        if modulus == 0 {
            return Err(GrError::InvalidDomain("challenge modulus must be nonzero"));
        }
        Ok(self.challenge_u64(label) % modulus)
    }

    pub fn derive_unique_positions(
        &self,
        label_prefix: &[u8],
        modulus: u64,
        count: u64,
    ) -> Result<Vec<u64>> {
        if modulus == 0 {
            return Err(GrError::InvalidDomain("position modulus must be nonzero"));
        }

        let target = count.min(modulus) as usize;
        let mut positions = Vec::with_capacity(target);
        let mut attempt = 0u64;
        while positions.len() < target {
            let mut label = Vec::with_capacity(label_prefix.len() + 8);
            label.extend_from_slice(label_prefix);
            label.extend_from_slice(&attempt.to_le_bytes());
            let candidate = self.challenge_index(&label, modulus)?;
            if !positions.contains(&candidate) {
                positions.push(candidate);
            }
            attempt = attempt.saturating_add(1);
        }
        Ok(positions)
    }

    pub fn challenge_teichmuller(&self, ctx: &GrContext, label: &[u8]) -> Result<GrElem> {
        for attempt in 0..4096u64 {
            let mut attempt_label = Vec::with_capacity(label.len() + 8);
            attempt_label.extend_from_slice(label);
            attempt_label.extend_from_slice(&attempt.to_le_bytes());
            let mut bytes = vec![0; ctx.elem_bytes()];
            self.challenge_bytes(&attempt_label, &mut bytes);
            let base = ctx.deserialize(&bytes)?;
            if !ctx.is_unit(&base) {
                continue;
            }
            let candidate = ctx.pow(&base, 1u128 << (ctx.config().k_exp - 1));
            if candidate != ctx.zero() && is_teichmuller_element(ctx, &candidate) {
                return Ok(candidate);
            }
        }

        Err(GrError::InvalidDomain(
            "failed to sample Teichmuller transcript challenge",
        ))
    }
}

fn absorb_labeled_bytes(hasher: &mut Hasher, label: &[u8], bytes: &[u8]) {
    hasher.update(&(label.len() as u64).to_le_bytes());
    hasher.update(label);
    hasher.update(&(bytes.len() as u64).to_le_bytes());
    hasher.update(bytes);
}

#[cfg(test)]
mod tests {
    use crate::{
        algebra::galois_ring::{is_teichmuller_element, GrConfig, GrContext},
        protocols::whir_gr::transcript::Transcript,
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
    fn transcript_challenges_should_be_deterministic() {
        let mut lhs = Transcript::new(b"test");
        let mut rhs = Transcript::new(b"test");
        lhs.absorb_labeled_bytes(b"root", b"abc");
        rhs.absorb_labeled_bytes(b"root", b"abc");

        assert_eq!(lhs.challenge_u64(b"alpha"), rhs.challenge_u64(b"alpha"));
    }

    #[test]
    fn transcript_labels_should_domain_separate_challenges() {
        let transcript = Transcript::new(b"test");

        assert_ne!(
            transcript.challenge_u64(b"alpha"),
            transcript.challenge_u64(b"beta")
        );
    }

    #[test]
    fn unique_positions_should_not_repeat() {
        let transcript = Transcript::new(b"positions");

        let positions = transcript.derive_unique_positions(b"query", 16, 8).unwrap();

        let mut sorted = positions.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), positions.len());
    }

    #[test]
    fn teichmuller_challenge_should_land_in_teichmuller_set() {
        let ctx = sample_context();
        let transcript = Transcript::new(b"teich");

        let challenge = transcript.challenge_teichmuller(&ctx, b"gamma").unwrap();

        assert!(is_teichmuller_element(&ctx, &challenge));
    }
}
