# WHIR_GR Rust Port Execution Plan

Status: draft implementation tracker
Branch: `feature/whir-gr`
Primary target: prototype Rust implementation of the C++ `whir_gr_ud` unique-decoding WHIR-over-Galois-ring PCS.

## 1. Scope

This plan targets a pure Rust prototype of the WHIR_GR protocol currently implemented in C++ under
`$HOME/STIR&WHIRoverGR`. The behavior reference is the C++ `whir_gr_ud` path, not the existing
finite-field WHIR implementation in this repository.

The first implementation target is:

- Base ring: `GR(2^s, r)` with `p = 2`.
- Protocol family: unique-decoding WHIR PCS over Galois rings, matching the paper/prototype scope
  represented by C++ `WhirProver`, `WhirVerifier`, and `select_whir_unique_decoding_parameters`.
- Implementation style: self-contained Rust modules, with C++ used only as reference material and
  parity oracle.
- Initial arithmetic representation: coefficients modulo `2^s` stored in `u64`, so the first
  supported range is `s <= 64`.

Non-goals for the prototype:

- No C++ or NTL FFI in the Rust runtime.
- No claim of production-ready security audit.
- No finite-field WHIR rewrite as the primary deliverable.
- No ZK WHIR extension.
- No Johnson/list-decoding WHIR variant unless explicitly added after the unique-decoding prototype.
- No generic computer algebra system for arbitrary Galois rings in phase 1.

## 2. Current Repository Baseline

The current Rust repository already provides useful infrastructure, but the WHIR_GR algebra and
protocol core must be added rather than obtained by a field-type substitution.

Reusable pieces:

- CLI and benchmark organization can host a new WHIR_GR benchmark command or binary.
- Existing hash, Merkle, and transcript patterns are useful as engineering references.
- Existing test and crate layout can host focused algebra/protocol tests.

Pieces that do not directly carry over:

- The current WHIR configuration assumes finite-field WHIR with power-of-two domains and binary
  folding. WHIR_GR needs ternary domains and repeated ternary folding.
- Existing algebra is built around `ark_ff::Field`. `GR(2^s, r)` is not a field because many
  nonzero elements are non-units, so it needs separate ring traits and APIs.
- Existing sumcheck paths are binary/quadratic. WHIR_GR needs the C++ multiquadratic and ternary
  constraint machinery.

## 3. C++ Reference Map

Use these C++ files as the behavior source while porting:

| C++ source | Rust target | Purpose |
| --- | --- | --- |
| `include/algebra/gr_context.hpp`, `src/algebra/gr_context.cpp` | `src/algebra/galois_ring/{mod.rs,context.rs,element.rs}` | Ring context, element representation, add/sub/mul/neg/pow, units, inverses, serialization |
| `include/algebra/teichmuller.hpp`, `src/algebra/teichmuller.cpp` | `src/algebra/galois_ring/teichmuller.rs` | Teichmuller lifts, generator selection, challenge sampling support |
| `include/domain.hpp`, `src/domain.cpp` | `src/algebra/galois_ring/domain.rs` | Multiplicative and affine domains over the ring |
| `include/crypto/*`, `src/crypto/*` | `src/protocols/whir_gr/transcript.rs`, Merkle helpers | Hash transcript, Merkle roots, challenge derivation |
| `include/whir/common.hpp`, `src/whir/common.cpp` | `src/protocols/whir_gr/common.rs` | Public parameters, proof structs, proof size accounting |
| `include/whir/multiquadratic.hpp`, `src/whir/multiquadratic.cpp` | `src/protocols/whir_gr/multiquadratic.rs` | `Pow_m`, multiquadratic evaluation, prefix restriction |
| `include/whir/constraint.hpp`, `src/whir/constraint.cpp` | `src/protocols/whir_gr/constraint.rs` | Sumcheck constraints and verifier identities |
| `include/whir/folding.hpp`, `src/whir/folding.cpp` | `src/protocols/whir_gr/folding.rs` | Repeated ternary folding and virtual fold query calculation |
| `include/whir/soundness.hpp`, `src/whir/soundness.cpp` | `src/protocols/whir_gr/soundness.rs` | Unique-decoding parameter selector and security accounting |
| `include/whir/prover.hpp`, `src/whir/prover.cpp` | `src/protocols/whir_gr/prover.rs` | Commit and open algorithms |
| `include/whir/verifier.hpp`, `src/whir/verifier.cpp` | `src/protocols/whir_gr/verifier.rs` | Verification algorithm |
| `bench/bench_time.cpp`, `bench/presets/whir.json` | `src/bin/whir_gr_benchmark.rs` or existing benchmark CLI | Benchmark and parity surface |
| `tests/test_whir_*.cpp`, `tests/test_domain.cpp`, `tests/test_crypto.cpp` | `tests/whir_gr_*.rs` or module tests | Acceptance and regression tests |

## 4. Rust Module Layout

Add the following modules:

```text
src/algebra/galois_ring/
  mod.rs
  context.rs
  element.rs
  poly_f2.rs
  teichmuller.rs
  domain.rs

src/protocols/whir_gr/
  mod.rs
  common.rs
  config.rs
  serialization.rs
  transcript.rs
  multiquadratic.rs
  constraint.rs
  folding.rs
  soundness.rs
  prover.rs
  verifier.rs

src/bin/whir_gr_benchmark.rs
```

Public API shape:

```rust
pub struct GrContext {
    pub k_exp: u32,
    pub degree: usize,
    // modulus 2^k, defining polynomial, Teichmuller data, serializer widths
}

pub struct GrElem {
    // coefficients in Z/(2^k)Z modulo the context polynomial
}

pub struct WhirGrPublicParameters { /* ring, dimensions, domains, folding schedule */ }
pub struct WhirGrCommitment { /* root and metadata */ }
pub struct WhirGrProof { /* all Merkle openings and ring payloads */ }

pub struct WhirGrProver<'a> { /* borrowed public parameters */ }
pub struct WhirGrVerifier<'a> { /* borrowed public parameters */ }
```

Prefer context-owned operations such as `ctx.add(&a, &b)` and `ctx.mul(&a, &b)` over embedding an
`Arc<GrContext>` in every element. This avoids high-frequency reference-count churn in hot loops and
makes serialization width explicit.

Use `Result<T, WhirGrError>` for fallible production APIs. `unwrap` and `expect` are acceptable only
inside tests or clearly unreachable internal invariants.

## 5. Phase Plan

### P0. Baseline and Tracking

Deliverables:

- Keep implementation work on `feature/whir-gr`.
- Commit this plan before starting code changes if the user wants a checkpoint.
- Record every phase completion in this document.

Acceptance:

- `git status --short --branch` shows `feature/whir-gr`.
- This document exists at `docs/whir_gr_rust_port_plan.md`.

### P1. Galois Ring Core

Deliverables:

- Implement `GrContext`, `GrElem`, and `poly_f2` support.
- Support `GR(2^s, r)` for `s <= 64`.
- Implement zero, one, add, sub, neg, mul, square, pow, equality, unit test, inverse for units,
  random element sampling, and canonical byte serialization.
- Add a small registry or deterministic generator for irreducible binary polynomials needed by
  the reference benchmark parameters.

Acceptance:

- Rust algebra tests match C++ behavior for small rings such as `GR(2^4, 2)`, `GR(2^8, 3)`, and
  one target benchmark ring.
- Ring axioms pass randomized tests for add/mul associativity, distributivity, identities, and
  serialization round trips.
- Unit inversion tests cover units and non-units.
- Command: `cargo test --lib galois_ring`.

Implementation notes:

- Do not implement `ark_ff::Field` for ring elements.
- Keep reduction modulo the defining polynomial structured, not string-based.
- For `s == 64`, avoid `1u64 << 64`; use explicit masking helpers.

### P2. Teichmuller and Domains

Deliverables:

- Port Teichmuller lift/generator logic.
- Implement deterministic challenge sampling into the Teichmuller set.
- Implement multiplicative and affine domain construction.
- Port C++ domain divisibility checks used by the unique-decoding selector.

Acceptance:

- Teichmuller generator order and domain sizes match C++ for selected `(s, r)`.
- `omega`, `omega^2`, and the ternary grid `{1, omega, omega^2}` match C++ serialized values.
- Challenge sampling is deterministic for fixed transcript bytes.
- Command: `cargo test --lib teichmuller domain`.

### P3. WHIR_GR Serialization, Transcript, and Merkle Layer

Deliverables:

- Define canonical byte encoding for `GrElem`, vectors, domains, Merkle leaves, commitments, and
  proof payloads.
- Add WHIR_GR transcript helpers for absorbing public parameters, roots, points, evaluations, and
  challenges.
- Add byte-oriented Merkle helpers or a WHIR_GR adapter around the existing Merkle code.

Acceptance:

- Merkle roots are deterministic for fixed ring payloads.
- Tampering with a leaf, path, root, or serialized proof component makes verification fail.
- Proof byte-size accounting is based on actual serialized bytes.
- Command: `cargo test --lib whir_gr_serialization whir_gr_merkle`.

### P4. Multilinear to Multiquadratic Layer

Deliverables:

- Port `Pow_m`, multiquadratic evaluation, and prefix restriction.
- Implement conversion from multilinear coefficients/evaluations to multiquadratic form as used by
  the C++ prover.
- Add deterministic polynomial/test-vector helpers.

Acceptance:

- Rust `Pow_m` and multiquadratic evaluations match C++ fixtures for small dimensions.
- Prefix restriction agrees with direct evaluation.
- Command: `cargo test --lib whir_gr_multiquadratic`.

### P5. Constraints and Ternary Sumcheck

Deliverables:

- Port WHIR_GR constraint objects and verifier identity checks.
- Implement the ternary grid based checks used by C++ verification.
- Keep challenge derivation byte-for-byte deterministic inside the Rust implementation.

Acceptance:

- Honest constraints pass for small dimensions.
- Tampered claimed evaluation, challenge, or intermediate polynomial fails.
- Command: `cargo test --lib whir_gr_constraint`.

### P6. Repeated Ternary Folding

Deliverables:

- Port repeated ternary fold table construction.
- Port virtual fold query index calculation.
- Port leaf-payload evaluation for folded values.

Acceptance:

- Rust output matches C++ `test_whir_folding.cpp` style cases for `b = 1`, `b = 2`, and `b = 3`.
- Virtual query positions are stable under fixed seeds.
- Command: `cargo test --lib whir_gr_folding`.

### P7. Unique-Decoding Parameter Selector

Deliverables:

- Port `WhirUniqueDecodingInputs`.
- Port input validation, candidate analysis, domain divisibility checks, and
  `select_whir_unique_decoding_parameters`.
- Expose the selected public parameters for the prover, verifier, and benchmark CLI.

Acceptance:

- Rust selector agrees with C++ for the reference preset rows.
- Invalid inputs reject with structured errors rather than panics.
- Command: `cargo test --lib whir_gr_soundness`.

### P8. Prover and Verifier

Deliverables:

- Implement `WhirGrProver::commit`.
- Implement `WhirGrProver::open`.
- Implement `WhirGrVerifier::verify`.
- Implement proof structs with explicit serialization and size accounting.

Acceptance:

- Honest round trip verifies for multiple small deterministic cases.
- Negative tests reject modified roots, evaluations, Merkle paths, folding payloads, and final values.
- Command: `cargo test --lib whir_gr_roundtrip`.

### P9. Benchmark and C++ Parity Surface

Deliverables:

- Add `src/bin/whir_gr_benchmark.rs` or an equivalent subcommand.
- Support the important C++ benchmark knobs: `p`, `k_exp`, `r`, `n`, security target, rate, seed,
  and repetition count.
- Emit CSV columns comparable to the C++ `bench_time.cpp` `whir_gr_ud` path:
  `protocol`, `p`, `k_exp`, `r`, `n`, `rate`, `lambda`, `effective_security_bits`,
  `commit_ms`, `open_ms`, `verify_ms`, `serialized_bytes_actual`.

Acceptance:

- A small Rust benchmark row can be compared against a C++ row with the same parameters.
- The benchmark clearly labels unsupported settings instead of silently changing parameters.
- Command: `cargo run --release --bin whir_gr_benchmark -- --help`.

### P10. Release Candidate Sweep

Deliverables:

- Run formatting, tests, lint, and at least one small parity benchmark.
- Update this document with completed statuses, known deviations, and next-phase work.

Acceptance:

```bash
cargo fmt --check
cargo test --all-targets
cargo clippy --all-targets --all-features --locked -- -D warnings
```

The prototype is acceptable when:

- The verifier accepts honest Rust proofs and rejects targeted tampering.
- The parameter selector agrees with C++ for selected reference inputs.
- Ring/domain/folding/multiquadratic tests have C++ parity coverage.
- Benchmark output can be joined with C++ `whir_gr_ud` rows by matching parameters.

## 6. Cross-Implementation Parity Strategy

Use C++ as a reference oracle in three layers:

1. Algebra fixtures: serialized ring elements, products, inverses, Teichmuller powers, domains.
2. Protocol fixtures: multiquadratic evaluations, fold payloads, selected public parameters.
3. End-to-end fixtures: proof acceptance, failure after tampering, proof size, benchmark rows.

Recommended fixture policy:

- Keep fixtures small and deterministic.
- Store only semantic expected values, not large opaque proof dumps unless needed for debugging.
- Prefer regeneratable fixture exporters in the C++ repo over manually copied constants.
- Label fixture files as C++ reference fixtures.

## 7. Open Decisions Before Heavy Coding

These are the only decisions that can materially change the implementation shape:

1. Pure Rust vs FFI: this plan assumes pure Rust. Choosing C++/NTL FFI would change P1, packaging,
   tests, and portability.
2. Initial precision range: this plan assumes `s <= 64`. Supporting `s > 64` requires a bigint or
   limb-vector representation before P1 is stable.
3. Benchmark surface: this plan recommends a new `whir_gr_benchmark` binary first, then integration
   into the existing benchmark CLI after the protocol is stable.
4. Performance order: this plan recommends correctness-first Horner/direct evaluation where simpler,
   then porting C++ `fft3` and cache optimizations only after parity tests are in place.
5. Transcript compatibility: this plan requires deterministic Rust transcript behavior, but not
   byte-for-byte compatibility with C++ unless explicitly requested. Byte-for-byte C++ compatibility
   would require stricter serialization and transcript fixture work in P3.

## 8. Main Risks

- Algebra risk: incorrect Galois-ring construction or Teichmuller lifting will invalidate every
  higher layer.
- Protocol risk: finite-field WHIR abstractions in this repository can look reusable but encode
  assumptions that are false over `GR(2^s, r)`.
- Transcript risk: proof acceptance depends on stable serialization and challenge ordering.
- Performance risk: a direct correctness implementation may be much slower than C++ until `fft3`,
  batch inversion, and cache-friendly evaluation paths are ported.
- Scope risk: the unique-decoding PCS prototype should not be described as full WHIR or production
  WHIR-over-GR.

## 9. Tracking Checklist

| Status | Phase | Deliverable | Acceptance |
| --- | --- | --- | --- |
| Done | P0 | Branch and plan document | `feature/whir-gr`, this file exists |
| Done | P1 | Galois ring core | `cargo test --lib galois_ring` |
| Done | P2 | Teichmuller and domains | `cargo test --lib galois_ring` |
| Done | P3 | Serialization, transcript, Merkle | `cargo test --lib whir_gr` |
| Done | P4 | Multiquadratic layer | `cargo test --lib multiquadratic` |
| Todo | P5 | Constraints and ternary sumcheck | `cargo test --lib whir_gr_constraint` |
| Todo | P6 | Repeated ternary folding | `cargo test --lib whir_gr_folding` |
| Todo | P7 | Unique-decoding selector | `cargo test --lib whir_gr_soundness` |
| Todo | P8 | Prover and verifier | `cargo test --lib whir_gr_roundtrip` |
| Todo | P9 | Benchmark parity surface | `cargo run --release --bin whir_gr_benchmark -- --help` |
| Todo | P10 | Release candidate sweep | fmt, tests, clippy, parity row |

## 10. Immediate Next Step

Start P5 by porting the WHIR_GR constraint and ternary sumcheck helper layer. Keep prover/verifier
code out of scope until P5 and P6 establish deterministic constraint and folding helpers.

## 11. Phase Review Log

### P1. Galois Ring Core

Status: complete in this branch; review gates passed.

Implemented:

- `src/algebra/galois_ring/context.rs`: `GrConfig`, `GrContext`, structured `GrError`, metadata,
  element construction, add/sub/neg/mul/square/pow, unit detection, Newton-lifted inverse,
  batch inverse, deterministic random sampling, and canonical little-endian serialization.
- `src/algebra/galois_ring/element.rs`: owned coefficient-vector `GrElem`.
- `src/algebra/galois_ring/poly_f2.rs`: binary polynomial arithmetic, irreducibility testing,
  deterministic irreducible polynomial selection, and GF(2) inverse modulo the defining
  polynomial.
- `src/algebra/mod.rs`: public `galois_ring` module export.

Review evidence:

- `cargo fmt --check`: passed.
- `cargo clippy --lib --all-features --locked -- -D warnings`: passed.
- `cargo test --lib galois_ring`: passed, 12 tests.
- `cargo test --lib`: passed, 124 passed and 25 ignored.
- `git diff --check`: passed.
- C++ reference smoke: `$HOME/STIR&WHIRoverGR/build/test_gr_basic` passed all tests.

Known boundary:

- P1 supports `p = 2` and `s <= 64`.
- The Rust ring uses deterministic irreducible binary polynomial selection. It is suitable for the
  Rust prototype and P1 algebra tests. Byte-for-byte C++ defining-polynomial compatibility remains a
  P2/P3 decision if strict C++ transcript/proof compatibility is required.

### P2. Teichmuller and Domains

Status: complete in this branch; review gates passed.

Implemented:

- `src/algebra/galois_ring/teichmuller.rs`: deterministic Teichmuller projection, subgroup-size
  support checks for divisibility by `2^r - 1`, subgroup generator search, exact-order checks,
  Teichmuller membership, subgroup enumeration, and index-based Teichmuller element access.
- `src/algebra/galois_ring/domain.rs`: `Domain` over `Arc<GrContext>` with subgroup/coset
  constructors, element access, full enumeration, membership check, Teichmuller-subset check,
  `scale`, `scale_offset`, `pow_map`, and disjointness checks.
- `src/algebra/galois_ring/context.rs`: limb-vector exponentiation helper for large
  `2^r - 1`-style exponents and structured errors for domain/subgroup failures.
- `src/algebra/galois_ring/mod.rs`: public exports for Teichmuller and domain APIs.

Review evidence:

- `cargo fmt --check`: passed.
- `cargo clippy --lib --all-features --locked -- -D warnings`: passed.
- `cargo test --lib galois_ring`: passed, 24 tests.
- `cargo test --lib`: passed, 136 passed and 25 ignored.
- `git diff --check`: passed.
- C++ reference smoke: `$HOME/STIR&WHIRoverGR/build/test_domain` passed all tests.

Confirmed boundary:

- The user confirmed Rust does not need byte-for-byte compatibility with C++ defining polynomials,
  Teichmuller generator bytes, domain root bytes, or serialized `omega` bytes.
- P2 therefore keeps Rust deterministic and protocol-self-consistent, with C++ used as behavior and
  test-shape reference rather than as a byte-level fixture oracle.

### P3. Serialization, Transcript, and Merkle

Status: complete in this branch; review gates passed.

Implemented:

- `src/protocols/whir_gr/common.rs`: WHIR_GR public-parameter, commitment, sumcheck polynomial,
  round proof, proof, and opening structs for later protocol phases.
- `src/protocols/whir_gr/serialization.rs`: canonical length-prefixed byte writer plus serializers
  for ring elements, ring vectors, domains, public parameters, sumcheck polynomials, Merkle proofs,
  round proofs, full proofs, and openings.
- `src/protocols/whir_gr/transcript.rs`: BLAKE3-based labeled transcript with deterministic byte
  challenges, index challenges, unique-position derivation, and Teichmuller challenge sampling.
- `src/protocols/whir_gr/merkle.rs`: byte-oriented Merkle tree over serialized ring payloads,
  opening proof generation, verification, tamper rejection, and oracle leaf/tree helpers.
- `src/protocols/mod.rs`: public `whir_gr` protocol module export.

Review evidence:

- `cargo fmt --check`: passed.
- `cargo clippy --lib --all-features --locked -- -D warnings`: passed.
- `cargo test --lib whir_gr`: passed, 10 tests.
- `cargo test --lib`: passed, 146 passed and 25 ignored.
- `git diff --check`: passed.
- C++ reference smoke: `$HOME/STIR&WHIRoverGR/build/test_crypto` passed all tests.

Confirmed boundary:

- P3 provides deterministic Rust-native transcript/Merkle behavior. It intentionally does not
  attempt byte-for-byte C++ transcript compatibility, matching the user's P2 clarification.

### P4. Multilinear and Multiquadratic Layer

Status: complete in this branch; review gates passed.

Implemented:

- `src/protocols/whir_gr/multiquadratic.rs`: checked `pow2`/`pow3`, base-3 little-endian
  encode/decode, `pow_m`, `MultiQuadraticPolynomial`, `MultilinearPolynomial`,
  multiquadratic evaluation, univariate `evaluate_pow`, prefix restriction, multilinear
  evaluation, and multilinear-to-multiquadratic embedding.
- `src/algebra/galois_ring/context.rs`: structured polynomial/overflow error variants used by
  the polynomial layer.
- `src/protocols/whir_gr/mod.rs`: public multiquadratic module export.

Review evidence:

- `cargo fmt --check`: passed.
- `cargo clippy --lib --all-features --locked -- -D warnings`: passed.
- `cargo test --lib multiquadratic`: passed, 8 tests.
- `cargo test --lib`: passed, 154 passed and 25 ignored.
- `git diff --check`: passed.
- C++ reference smoke: `$HOME/STIR&WHIRoverGR/build/test_whir_multiquadratic` passed all tests.
- C++ reference smoke: `$HOME/STIR&WHIRoverGR/build/test_whir_multilinear` passed all tests.

Confirmed boundary:

- P4 uses Rust-native ring/domain choices from P1/P2. It checks structural behavior and algebraic
  identities rather than byte-for-byte C++ coefficient fixtures.
