use serde::Serialize;

use crate::{
    algebra::galois_ring::{is_teichmuller_element, Domain, GrContext, GrElem, GrError, Result},
    transcript::{
        codecs::Empty, DomainSeparator, Encoding, NargDeserialize, ProverState, VerificationError,
        VerificationResult, VerifierMessage,
    },
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TranscriptFrame {
    label: Vec<u8>,
    payload: Vec<u8>,
}

#[derive(Serialize)]
struct TranscriptDomain {
    protocol: &'static str,
    domain: Vec<u8>,
}

pub struct Transcript {
    state: ProverState,
}

impl TranscriptFrame {
    pub fn new(label: &[u8], payload: &[u8]) -> Self {
        Self {
            label: label.to_vec(),
            payload: payload.to_vec(),
        }
    }

    pub fn label(&self) -> &[u8] {
        &self.label
    }

    pub fn payload(&self) -> &[u8] {
        &self.payload
    }
}

impl Encoding<[u8]> for TranscriptFrame {
    fn encode(&self) -> impl AsRef<[u8]> {
        let mut encoded = Vec::with_capacity(16 + self.label.len() + self.payload.len());
        encoded.extend_from_slice(&(self.label.len() as u64).to_le_bytes());
        encoded.extend_from_slice(&self.label);
        encoded.extend_from_slice(&(self.payload.len() as u64).to_le_bytes());
        encoded.extend_from_slice(&self.payload);
        encoded
    }
}

impl NargDeserialize for TranscriptFrame {
    fn deserialize_from_narg(buf: &mut &[u8]) -> VerificationResult<Self> {
        let label = read_len_prefixed(buf)?;
        let payload = read_len_prefixed(buf)?;
        Ok(Self { label, payload })
    }
}

impl Transcript {
    pub fn new(domain: &[u8]) -> Self {
        let instance = Empty;
        let config = TranscriptDomain {
            protocol: "whir-gr.transcript.v1",
            domain: domain.to_vec(),
        };
        let domain_separator = DomainSeparator::protocol(&config).instance(&instance);
        Self {
            state: ProverState::new_std(&domain_separator),
        }
    }

    pub fn absorb_labeled_bytes(&mut self, label: &[u8], bytes: &[u8]) {
        self.state
            .public_message(&TranscriptFrame::new(label, bytes));
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

    pub fn challenge_bytes(&mut self, label: &[u8], out: &mut [u8]) {
        self.absorb_labeled_bytes(b"challenge", label);
        let mut filled = 0usize;
        let mut block_index = 0u64;
        while filled < out.len() {
            self.absorb_labeled_bytes(b"challenge.block", &block_index.to_le_bytes());
            let block: [u8; 32] = self.state.verifier_message();
            let take = (out.len() - filled).min(block.len());
            out[filled..filled + take].copy_from_slice(&block[..take]);
            filled += take;
            block_index = block_index.saturating_add(1);
        }
    }

    pub fn challenge_u64(&mut self, label: &[u8]) -> u64 {
        let mut bytes = [0; 8];
        self.challenge_bytes(label, &mut bytes);
        u64::from_le_bytes(bytes)
    }

    pub fn challenge_index(&mut self, label: &[u8], modulus: u64) -> Result<u64> {
        if modulus == 0 {
            return Err(GrError::InvalidDomain("challenge modulus must be nonzero"));
        }
        Ok(self.challenge_u64(label) % modulus)
    }

    pub fn derive_unique_positions(
        &mut self,
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

    pub fn challenge_teichmuller(&mut self, ctx: &GrContext, label: &[u8]) -> Result<GrElem> {
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

    pub fn challenge_teichmuller_outside_domain(
        &mut self,
        ctx: &GrContext,
        label: &[u8],
        domain: &Domain,
    ) -> Result<GrElem> {
        if domain.context().config() != ctx.config() {
            return Err(GrError::DifferentRings);
        }

        for attempt in 0..4096u64 {
            let mut attempt_label = Vec::with_capacity(label.len() + 8);
            attempt_label.extend_from_slice(label);
            attempt_label.extend_from_slice(&attempt.to_le_bytes());
            let candidate = self.challenge_teichmuller(ctx, &attempt_label)?;
            if !domain.contains(&candidate) {
                return Ok(candidate);
            }
        }

        Err(GrError::InvalidDomain(
            "failed to sample out-of-domain Teichmuller transcript challenge",
        ))
    }
}

fn read_len_prefixed(buf: &mut &[u8]) -> VerificationResult<Vec<u8>> {
    let len = read_u64(buf)?;
    let len = usize::try_from(len).map_err(|_| VerificationError)?;
    if buf.len() < len {
        return Err(VerificationError);
    }
    let (head, tail) = buf.split_at(len);
    *buf = tail;
    Ok(head.to_vec())
}

fn read_u64(buf: &mut &[u8]) -> VerificationResult<u64> {
    if buf.len() < 8 {
        return Err(VerificationError);
    }
    let (head, tail) = buf.split_at(8);
    *buf = tail;
    let bytes: [u8; 8] = head.try_into().map_err(|_| VerificationError)?;
    Ok(u64::from_le_bytes(bytes))
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::{
        algebra::galois_ring::{is_teichmuller_element, Domain, GrConfig, GrContext},
        protocols::whir_gr::transcript::{Transcript, TranscriptFrame},
        transcript::{NargDeserialize, NargSerialize},
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
        let mut lhs = Transcript::new(b"test");
        let mut rhs = Transcript::new(b"test");

        assert_ne!(lhs.challenge_u64(b"alpha"), rhs.challenge_u64(b"beta"));
    }

    #[test]
    fn unique_positions_should_not_repeat() {
        let mut transcript = Transcript::new(b"positions");

        let positions = transcript.derive_unique_positions(b"query", 16, 8).unwrap();

        let mut sorted = positions.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), positions.len());
    }

    #[test]
    fn teichmuller_challenge_should_land_in_teichmuller_set() {
        let ctx = sample_context();
        let mut transcript = Transcript::new(b"teich");

        let challenge = transcript.challenge_teichmuller(&ctx, b"gamma").unwrap();

        assert!(is_teichmuller_element(&ctx, &challenge));
    }

    #[test]
    fn ood_teichmuller_challenge_should_avoid_domain() {
        let ctx = Arc::new(sample_context());
        let domain = Domain::teichmuller_subgroup(Arc::clone(&ctx), 9).unwrap();
        let mut transcript = Transcript::new(b"ood");

        let challenge = transcript
            .challenge_teichmuller_outside_domain(&ctx, b"eta", &domain)
            .unwrap();

        assert!(is_teichmuller_element(&ctx, &challenge));
        assert!(!domain.contains(&challenge));
    }

    #[test]
    fn transcript_frame_should_roundtrip_through_narg() {
        let frame = TranscriptFrame::new(b"label", b"payload");
        let mut bytes = Vec::new();
        frame.serialize_into_narg(&mut bytes);

        let mut read = bytes.as_slice();
        let decoded = TranscriptFrame::deserialize_from_narg(&mut read).unwrap();

        assert_eq!(decoded.label(), b"label");
        assert_eq!(decoded.payload(), b"payload");
        assert!(read.is_empty());
    }
}
