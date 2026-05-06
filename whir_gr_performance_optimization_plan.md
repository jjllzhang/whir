# WHIR-GR Performance Optimization Plan

Date: 2026-05-06

This plan targets the Rust WHIR-over-GR path under `src/protocols/whir_gr`.
It is based on the benchmark artifact `results/whir_gr_bench_20260506_000223.txt`,
the current Rust implementation, and the already-proven optimization directions
from the C++ WHIR-over-GR prototype in `/home/zjl/STIR&WHIRoverGR`.

The immediate symptom is that `whir_gr_commit` reaches 2.461 h at
`gr216_r162_m9_multilinear` and does not complete the printed `m10` row. The
benchmark only recorded commit timing and allocator statistics; it did not
record open time, verify time, or proof size.

## Scope

Optimize the Rust implementation of the unique-decoding WHIR-GR PCS over
`GR(2^16, r)` for the multilinear benchmark cases currently used by
`benches/whir_gr.rs`.

In scope:

- `WhirGrProver::commit_multilinear`
- `WhirGrProver::open`
- `WhirGrVerifier::verify`
- WHIR-GR benchmark/profiling surfaces
- Galois-ring arithmetic and serialization helpers used by WHIR-GR
- Merkle tree construction and proof-size reporting for WHIR-GR

Out of scope:

- Changing protocol semantics, transcript labels, query policy, or proof format
- Changing selected benchmark parameters to make numbers look better
- Implementing finite-field WHIR, ZK WHIR, or Johnson/list-decoding WHIR
- Removing correctness checks just to improve timings

## Current Evidence

The current commit benchmark behaves as follows:

| case | n | mean commit time | total alloc | peak live alloc |
|---|---:|---:|---:|---:|
| `m4` | 189 | 143.9 ms | 38 MB | 150.4 KB |
| `m5` | 513 | 989.7 ms | 264.1 MB | 422.2 KB |
| `m6` | 1539 | 8.794 s | 2.249 GB | 1.51 MB |
| `m7` | 4617 | 1.417 min | 19.83 GB | 3.744 MB |
| `m8` | 13203 | 16.43 min | 168.8 GB | 10.63 MB |
| `m9` | 39609 | 2.461 h | 1.515 TB | 32.37 MB |
| `m10` | 124173 | incomplete | incomplete | incomplete |

Interpretation:

- The bottleneck is not persistent memory usage. Peak live allocation remains
  small compared with total allocated bytes.
- The large `total alloc` values show repeated short-lived heap allocations.
- Runtime grows mainly because the current commit path performs dense Horner
  evaluation over a large ternary coefficient space and each Galois-ring
  operation allocates.

Relevant current code paths:

- `benches/whir_gr.rs`: ignored Divan benchmarks for commit/open/verify.
- `src/protocols/whir_gr/prover.rs`: `commit`, `commit_multilinear`,
  `encode_oracle`, and `open`.
- `src/protocols/whir_gr/multiquadratic.rs`: multilinear-to-multiquadratic
  embedding and dense Horner evaluation.
- `src/algebra/galois_ring/context.rs`: `add`, `mul`, `square`, `serialize`,
  `deserialize`.
- `src/protocols/whir_gr/merkle.rs`: oracle leaf serialization and Merkle tree
  construction.

## Root Causes

1. `commit_multilinear` first embeds a `2^m` multilinear polynomial into a
   sparse subset of a `3^m` multiquadratic coefficient vector.

2. `encode_oracle` then treats that representation as dense and evaluates it at
   every domain point with Horner evaluation.

3. For `m9`, this means roughly:

   ```text
   n = 39609
   dense coefficient count ~= (3^9 + 1) / 2 = 9842
   Horner steps ~= 39609 * 9842 = 3.9e8
   ```

   Each step performs a Galois-ring multiply and add.

4. `GrElem` stores coefficients in `Vec<u64>`. `GrContext::add` allocates a new
   vector for every addition, and `GrContext::mul` allocates a temporary vector
   of length about `2r - 1`. For `r = 162`, this is expensive.

5. `Domain::element(index)` computes `root^index` for each point. Sequential
   encoding should instead advance by repeated multiplication from a chunk
   start.

6. Merkle construction serializes oracle values into `Vec<Vec<u8>>`, then hashes
   them. This adds copying and allocation after the expensive encode phase.

7. Divan allocation profiling records every allocation. That is useful for
   diagnosis, but it can distort wall-clock timing once a case performs hundreds
   of millions of allocations.

## Execution Rules

- Optimize one item at a time.
- For every item, record before/after numbers on small cases before running
  larger cases.
- Do not use `m8+` as inner-loop validation.
- Treat proof bytes and transcript output as compatibility contracts unless the
  change explicitly documents a proof-format change.
- Keep public API changes narrow and aligned with current module boundaries.
- Prefer borrowed slices and reusable buffers over cloning or new owned vectors
  in hot loops.
- Add regression tests before replacing a reference implementation path.

## Phase 0: Measurement Harness

Goal: make performance work measurable before changing hot code.

### Tasks

1. Add a focused profiling binary, for example `src/bin/whir_gr_profile.rs`.

   Required CLI options:

   ```text
   --case m4|m5|m6|m7|m8|m9|m10
   --phase commit|open|verify|roundtrip
   --reps N
   --allocator-stats true|false
   --format csv|text|json
   ```

2. Reuse the benchmark cases from `benches/whir_gr.rs` without duplicating
   unchecked constants. A small shared helper module is acceptable, but keep it
   private to benchmark/profile surfaces if it is not a library API.

3. Report at least these columns:

   ```text
   case,k_exp,r,n,variable_count,max_layer_width,lambda_target,rho0,
   phase,reps,commit_ms,open_ms,verify_ms,encode_oracle_ms,merkle_ms,
   to_multiquadratic_ms,open_fold_ms,verify_algebra_ms,serialized_opening_bytes,
   accepted
   ```

   If a bucket is not yet available, print an empty field rather than a fake
   zero.

4. Add proof-size reporting by serializing `WhirGrOpening` with
   `serialize_opening` and printing `serialized.len()`.

5. Keep the existing Divan bench, but document the safer commands:

   ```bash
   cargo run --release --bin whir_gr_profile -- --case m4 --phase commit --reps 3 --format csv
   cargo run --release --bin whir_gr_profile -- --case m4 --phase open --reps 3 --format csv
   cargo run --release --bin whir_gr_profile -- --case m4 --phase verify --reps 10 --format csv
   ```

   For allocation profiling on small cases only:

   ```bash
   cargo bench --bench whir_gr -- --ignored whir_gr_commit_small --max-time 1
   ```

   Do not use the full ignored Divan bench as the default inner loop:

   ```bash
   cargo bench --bench whir_gr -- --ignored whir_gr_commit
   ```

   Divan's function-name filter does not filter individual `WhirGrBenchCase`
   arguments, so this command enters `m7+`, `m9`, and `m10`. Use it only when a
   full allocation-profile sweep is intended.

### Validation

Run:

```bash
cargo fmt --check
cargo test --lib whir_gr
cargo run --release --bin whir_gr_profile -- --case m4 --phase roundtrip --reps 1 --format csv
```

Acceptance:

- `m4 roundtrip` verifies successfully.
- The profile output includes non-empty `commit_ms`, `open_ms`, `verify_ms`,
  and `serialized_opening_bytes`.
- The profile command can run one selected case without entering later cases.

### Phase 0 Review

Status: complete.

Implemented files:

- `src/protocols/whir_gr/bench_support.rs`
- `src/bin/whir_gr_profile.rs`
- `benches/whir_gr.rs`
- `src/protocols/whir_gr/mod.rs`

Protocol-drift review:

- No changes were made to `WhirGrProver::commit`, `WhirGrProver::open`,
  `WhirGrVerifier::verify`, transcript labels, Merkle hashing, serialization,
  soundness selection, query selection, or proof structures.
- The original Divan bench case constants and input-construction helpers were
  moved into `bench_support` and reused by both `benches/whir_gr.rs` and the
  new profile binary.
- The profile binary calls the existing `commit_multilinear`, `open`, `verify`,
  and `serialize_opening` APIs. It does not implement an alternate prover,
  verifier, serializer, transcript, or parameter selector.
- Internal bucket columns such as `encode_oracle_ms` and `merkle_ms` are printed
  as empty fields until instrumentation exists, matching the Phase 0 rule to
  avoid fake zeros. Phase 2 added commit-side instrumentation for
  `encode_oracle_ms`, `merkle_ms`, and `to_multiquadratic_ms`.
- `--allocator-stats` is accepted for command stability, but detailed allocator
  tallies remain delegated to the existing Divan `AllocProfiler` bench.

Validation evidence:

```bash
cargo fmt --check
cargo test --lib whir_gr
cargo run --release --bin whir_gr_profile -- --case m4 --phase roundtrip --reps 1 --format csv
cargo clippy --all-targets --all-features --locked -- -D warnings
git diff --check
```

Observed `m4` profile output included non-empty `commit_ms`, `open_ms`,
`verify_ms`, `serialized_opening_bytes=162788`, and `accepted=true`. A later
Phase 3 check confirmed that Divan argument filtering is not case-selective, so
the safe inner-loop commands are the per-case `whir_gr_profile` commands above.

## Phase 1: Multilinear-Specialized Commit Encoding

Goal: avoid dense `3^m` Horner evaluation for multilinear inputs.

### Design

Add a multilinear-specific encoding path:

```rust
fn encode_multilinear_oracle(
    params: &WhirGrPublicParameters,
    domain: &Domain,
    polynomial: &MultilinearPolynomial,
) -> Result<Vec<GrElem>>
```

For each domain element `x`, evaluate the multilinear polynomial at:

```text
(x, x^3, x^(3^2), ..., x^(3^(m-1)))
```

This should match the current embedding semantics, because the current
multilinear-to-multiquadratic embedding maps a binary monomial index into the
corresponding ternary monomial index.

Use folded multilinear evaluation, not coefficient-by-coefficient monomial
reconstruction:

```text
scratch = coefficients padded to 2^m
for y_i in (x, x^3, ..., x^(3^(m-1))):
    scratch[j] = scratch[2*j] + y_i * scratch[2*j + 1]
return scratch[0]
```

Expected effect:

- For `m9`, dense Horner uses about 9842 coefficient steps per domain point.
- Folded multilinear evaluation uses 511 multiply-add steps per domain point.
- This reduces the dominant arithmetic count by about 19x before lower-level
  arithmetic improvements.

### Tasks

1. Add a reference test that compares the new multilinear encoder with the
   current `polynomial.to_multi_quadratic(...); encode_oracle(...)` path for
   small cases:

   - `m1`, `m2`, `m3`, `m4`
   - at least two deterministic coefficient seeds
   - at least one `GR(2^16, r)` context compatible with the WHIR-GR path

2. Implement the new encoder in `src/protocols/whir_gr/prover.rs` or a small
   helper module under `src/protocols/whir_gr`.

3. Change `commit_multilinear` to use the specialized encoder directly instead
   of converting to a dense `MultiQuadraticPolynomial` just to encode the first
   oracle.

4. Preserve `WhirGrCommitmentState` semantics. If `open` still needs a
   `MultiQuadraticPolynomial`, construct and store it after the initial oracle
   has been encoded, or delay construction until `open`. Measure both choices
   before deciding.

5. Keep the generic `commit(&MultiQuadraticPolynomial)` path unchanged unless a
   later phase optimizes it independently.

### Validation

Run:

```bash
cargo fmt --check
cargo test --lib whir_gr
cargo run --release --bin whir_gr_profile -- --case m4 --phase commit --reps 3 --format csv
cargo run --release --bin whir_gr_profile -- --case m5 --phase commit --reps 3 --format csv
cargo run --release --bin whir_gr_profile -- --case m6 --phase commit --reps 1 --format csv
```

Acceptance:

- All old WHIR-GR tests pass.
- The new encoder matches the old dense encoder on small cases.
- Commitment roots for the same small inputs are unchanged.
- `m6 commit` improves materially before moving to `m7+`.

### Phase 1 Review

Status: complete.

Implemented files:

- `src/protocols/whir_gr/prover.rs`
- `whir_gr_performance_optimization_plan.md`

Protocol-drift review:

- The generic `commit(&MultiQuadraticPolynomial)` path is unchanged.
- `commit_multilinear` now encodes the initial oracle with
  `encode_multilinear_oracle`, but it still stores the same embedded
  `MultiQuadraticPolynomial` in `WhirGrCommitmentState`; `open` therefore sees
  the same polynomial representation as before.
- The new test compares the optimized multilinear oracle against the old dense
  `to_multi_quadratic(...); encode_oracle(...)` reference for `m1`, `m2`,
  `m3`, and the benchmark `m4` case, using two deterministic seeds.
- The same test compares `commit_multilinear` roots against the dense
  `commit(&embedded)` roots. This is the key transcript compatibility check,
  because the commitment root is the first committed transcript message.
- No transcript labels, Merkle hashing, serialization format, soundness
  selection, query policy, opening proof structure, or verifier logic changed.

Validation evidence:

```bash
cargo fmt --check
cargo test --lib whir_gr
cargo run --release --bin whir_gr_profile -- --case m4 --phase commit --reps 3 --format csv
cargo run --release --bin whir_gr_profile -- --case m5 --phase commit --reps 3 --format csv
cargo run --release --bin whir_gr_profile -- --case m6 --phase commit --reps 1 --format csv
cargo clippy --all-targets --all-features --locked -- -D warnings
git diff --check
```

Observed commit timings:

| case | before | after Phase 1 | speedup |
|---|---:|---:|---:|
| `m4` | 143.9 ms | 95.2 ms | 1.51x |
| `m5` | 989.7 ms | 394.0 ms | 2.51x |
| `m6` | 8.794 s | 1.978 s | 4.45x |

The speedup increases with `m`, which is consistent with replacing sparse dense
Horner over the ternary embedding by folded multilinear evaluation.

## Phase 2: Sequential Domain Iteration

Goal: remove repeated exponentiation from oracle encoding.

### Tasks

1. Add a domain iterator or helper:

   ```rust
   impl Domain {
       pub fn iter_elements(&self) -> impl Iterator<Item = GrElem> + '_
   }
   ```

   It should start at `offset` and repeatedly multiply by `root`.

2. For chunked/parallel encoding, add a helper that computes one starting point
   with `offset * root^begin` and then advances sequentially inside the chunk.

3. Update `encode_oracle` and `encode_multilinear_oracle` to use sequential
   advancement instead of `domain.element(index)` inside the innermost loop.

4. Keep `Domain::element(index)` for random access and tests.

### Validation

Run:

```bash
cargo test --lib galois_ring::domain
cargo test --lib whir_gr
cargo run --release --bin whir_gr_profile -- --case m5 --phase commit --reps 3 --format csv
```

Acceptance:

- Iterated elements match `Domain::element(i)` for representative domains.
- Commitment roots remain unchanged.
- `commit_ms` decreases measurably. Phase 2 also populates `encode_oracle_ms`
  so later phases can compare the encoder bucket directly.

### Phase 2 Review

Status: complete.

Implemented files:

- `src/algebra/galois_ring/domain.rs`
- `src/protocols/whir_gr/prover.rs`
- `src/bin/whir_gr_profile.rs`
- `whir_gr_performance_optimization_plan.md`

Protocol-drift review:

- `Domain::element(index)` remains unchanged and is still used for sparse
  Merkle-query and verifier paths where random access is semantically needed.
- `Domain::iter_elements()` and `iter_elements_from(begin)` only enumerate the
  same coset sequence by repeated multiplication from `offset * root^begin`.
  They do not change domain construction, subgroup checks, offsets, roots, or
  sizes.
- `encode_oracle` and `encode_multilinear_oracle` now scan whole domains with
  `iter_elements()` instead of repeatedly calling `element(index)`. The oracle
  order is unchanged because the new domain test compares full iterator output
  against `element(i)` for both a subgroup and a coset.
- The existing multilinear encoder test still compares optimized commitment
  roots against dense `commit(&embedded)` roots. This covers the transcript
  compatibility point affected by oracle ordering.
- The new `commit_multilinear_profiled` method only wraps the same
  `commit_multilinear` implementation with timing buckets. It does not alter
  transcript labels, Merkle leaf bytes, proof format, soundness selection,
  query selection, `open`, or `verify`.

Validation evidence:

```bash
cargo fmt --check
cargo test --lib galois_ring::domain
cargo test --lib whir_gr
cargo run --release --bin whir_gr_profile -- --case m5 --phase commit --reps 3 --format csv
cargo run --release --bin whir_gr_profile -- --case m6 --phase commit --reps 1 --format csv
```

Observed commit timings:

| case | after Phase 1 | after Phase 2 | speedup |
|---|---:|---:|---:|
| `m5` | 394.0 ms | 311.3 ms | 1.27x |
| `m6` | 1.978 s | 1.679 s | 1.18x |

Observed Phase 2 commit bucket timings:

| case | `encode_oracle_ms` | `merkle_ms` | `to_multiquadratic_ms` |
|---|---:|---:|---:|
| `m5` | 300.7 ms | 0.71 ms | 0.18 ms |
| `m6` | 1.652 s | 16.4 ms | 0.59 ms |

The commit path is still dominated by oracle encoding, which makes Phase 3
in-place arithmetic the next highest-leverage step.

## Phase 3: In-Place Galois-Ring Arithmetic

Goal: remove heap allocation from hot arithmetic loops.

### Design

Introduce scratch-based arithmetic APIs in `GrContext`:

```rust
pub fn add_into(&self, out: &mut GrElem, lhs: &GrElem, rhs: &GrElem);
pub fn sub_into(&self, out: &mut GrElem, lhs: &GrElem, rhs: &GrElem);
pub fn mul_into(&self, out: &mut GrElem, lhs: &GrElem, rhs: &GrElem, scratch: &mut [u64]);
pub fn square_into(&self, out: &mut GrElem, value: &GrElem, scratch: &mut [u64]);
pub fn serialize_into(&self, out: &mut [u8], value: &GrElem);
```

The current allocating APIs can remain as compatibility wrappers. Hot paths
should use the `*_into` versions.

### Tasks

1. Add an internal scratch-size helper:

   ```rust
   pub const fn mul_scratch_len(&self) -> usize
   ```

   The required length is `2 * r - 1`.

2. Add tests comparing `add_into`, `sub_into`, `mul_into`, `square_into`, and
   `serialize_into` against the current allocating APIs.

3. Convert the following first:

   - `MultiQuadraticPolynomial::evaluate_pow`
   - the new multilinear oracle encoder
   - `Domain` sequential element advancement

4. Convert secondary hot paths in the phase that owns their protocol surface
   after commit improves:

   - `fold_eval`
   - `ternary_fold_table`
   - `evaluate_repeated_ternary_fold_from_values`
   - verifier equality-polynomial evaluation
   - Merkle leaf serialization

   The first three belong with Phase 5 open-path work; verifier equality
   belongs with Phase 6; Merkle leaf serialization can be pulled earlier if
   `merkle_ms` becomes material.

5. Avoid changing public representation in this phase. Keep `GrElem { coefficients:
   Vec<u64> }` until the scratch APIs are proven.

### Validation

Run:

```bash
cargo fmt --check
cargo test --lib galois_ring
cargo test --lib whir_gr
cargo run --release --bin whir_gr_profile -- --case m5 --phase commit --reps 3 --format csv
cargo run --release --bin whir_gr_profile -- --case m6 --phase commit --reps 1 --format csv
```

Acceptance:

- Arithmetic tests pass against old APIs.
- Commitment roots remain unchanged.
- `total alloc` in Divan drops substantially on `m5` and `m6`.
- Timing improves without depending on allocator profiler artifacts.

### Phase 3 Review

Status: complete for the commit hot loop. Secondary open/verifier conversions
remain assigned to Phase 5 and Phase 6, where their protocol surfaces are
reviewed directly.

Implemented files:

- `src/algebra/galois_ring/element.rs`
- `src/algebra/galois_ring/context.rs`
- `src/algebra/galois_ring/domain.rs`
- `src/protocols/whir_gr/multiquadratic.rs`
- `src/protocols/whir_gr/prover.rs`
- `src/protocols/whir_gr/bench_support.rs`
- `benches/whir_gr.rs`
- `whir_gr_performance_optimization_plan.md`

Protocol-drift review:

- Added `add_into`, `sub_into`, `mul_into`, `square_into`,
  `mul_base_scalar_into`, and `serialize_into` as allocation-reducing
  equivalents of existing arithmetic/serialization operations. The original
  allocating APIs remain available.
- Added tests comparing the in-place arithmetic APIs against the old allocating
  APIs, including the new base-scalar multiplication path.
- `Domain::iter_elements()` still matches `Domain::element(i)`, and now uses
  `mul_into` internally for sequential advancement.
- `MultiQuadraticPolynomial::evaluate_pow` uses scratch arithmetic but keeps
  the same Horner order.
- The multilinear oracle encoder uses scratch arithmetic and a base-scalar fast
  path only when a coefficient is provably in the base ring. Otherwise it falls
  back to full `mul_into`.
- Existing WHIR-GR tests still compare optimized multilinear oracle output and
  commitment roots against the dense reference path.
- No transcript labels, Merkle leaf format, Merkle hashing, proof structure,
  soundness selection, query selection, opening transcript, or verifier logic
  changed.

Validation evidence:

```bash
cargo fmt --check
cargo test --lib galois_ring
cargo test --lib whir_gr
cargo run --release --bin whir_gr_profile -- --case m5 --phase commit --reps 3 --format csv
cargo run --release --bin whir_gr_profile -- --case m6 --phase commit --reps 1 --format csv
cargo bench --bench whir_gr -- --ignored whir_gr_commit_small --max-time 1
```

Observed commit timings:

| case | after Phase 2 | after Phase 3 | speedup |
|---|---:|---:|---:|
| `m5` | 311.3 ms | 202.1 ms | 1.54x |
| `m6` | 1.679 s | 999.1 ms | 1.68x |

Observed Phase 3 commit bucket timings:

| case | `encode_oracle_ms` | `merkle_ms` | `to_multiquadratic_ms` |
|---|---:|---:|---:|
| `m5` | 191.4 ms | 0.65 ms | 0.18 ms |
| `m6` | 972.1 ms | 16.5 ms | 0.58 ms |

Observed `whir_gr_commit_small` Divan allocation totals:

| case | original total alloc | after Phase 3 total alloc | reduction |
|---|---:|---:|---:|
| `m4` | 38 MB | 8.046 MB | 4.7x |
| `m5` | 264.1 MB | 29.91 MB | 8.8x |
| `m6` | 2.249 GB | 152.5 MB | 14.7x |

The remaining commit time is still dominated by `encode_oracle_ms`, but the
dominant cost has shifted toward actual ring multiplication rather than heap
allocation churn.

## Phase 4: Parallel Encoding and Merkle Construction

Goal: use multiple cores after allocation pressure has been reduced.

### Tasks

1. Add a `rayon`-backed chunked encoder behind the existing `parallel` feature.

2. Use enough work per chunk to avoid scheduler overhead. Start with:

   ```text
   target_chunks = rayon::current_num_threads() * 4
   chunk_size = ceil(domain_size / target_chunks)
   ```

3. Each worker should allocate its own reusable scratch buffers.

4. Add deterministic tests comparing single-thread and parallel oracle outputs.

5. Parallelize Merkle leaf hashing and parent hashing only after measuring
   `merkle_ms`. Do not prioritize this if `encode_oracle_ms` still dominates.

6. Avoid global mutable caches in this phase unless profiling shows plan
   construction is material.

### Validation

Run:

```bash
RAYON_NUM_THREADS=1 cargo run --release --bin whir_gr_profile -- --case m6 --phase commit --reps 1 --format csv
RAYON_NUM_THREADS=8 cargo run --release --bin whir_gr_profile -- --case m6 --phase commit --reps 1 --format csv
cargo test --lib whir_gr
```

Acceptance:

- Oracle outputs and commitment roots are identical for 1 and 8 threads.
- `m6 commit` improves when `RAYON_NUM_THREADS=8`.
- `m4` does not regress badly from parallel overhead; if it does, add a size
  threshold before parallel dispatch.

### Phase 4 Review

Status: complete for oracle encoding. Merkle parallelization is intentionally
deferred because `merkle_ms` is still small compared with `encode_oracle_ms`.

Implemented files:

- `src/protocols/whir_gr/prover.rs`
- `whir_gr_performance_optimization_plan.md`

Protocol-drift review:

- The parallel path is behind the existing `parallel` feature and is used only
  when the initial domain size is at least `1024` and Rayon has more than one
  worker thread.
- Each Rayon chunk starts from `Domain::iter_elements_from(begin)` and uses its
  own scratch buffers. Chunks are collected in deterministic start-order before
  building the oracle vector.
- The sequential encoder remains available as
  `encode_multilinear_oracle_sequential`, and the new test compares direct
  parallel output against direct sequential output on the `m6` benchmark case.
- The commitment path still builds the same Merkle tree from the same ordered
  oracle values. No transcript labels, Merkle leaf bytes, Merkle hashing,
  proof format, soundness selection, query policy, `open`, or `verify` changed.
- `m4` remains below the parallel threshold, so small-case latency does not pay
  Rayon scheduling overhead.

Validation evidence:

```bash
cargo test --lib whir_gr
RAYON_NUM_THREADS=1 cargo run --release --bin whir_gr_profile -- --case m6 --phase commit --reps 1 --format csv
RAYON_NUM_THREADS=8 cargo run --release --bin whir_gr_profile -- --case m6 --phase commit --reps 1 --format csv
RAYON_NUM_THREADS=1 cargo run --release --bin whir_gr_profile -- --case m4 --phase commit --reps 3 --format csv
RAYON_NUM_THREADS=8 cargo run --release --bin whir_gr_profile -- --case m4 --phase commit --reps 3 --format csv
```

Observed commit timings:

| case | threads | commit time | `encode_oracle_ms` |
|---|---:|---:|---:|
| `m6` | 1 | 988.2 ms | 975.6 ms |
| `m6` | 8 | 157.3 ms | 144.3 ms |
| `m4` | 1 | 53.8 ms | 41.3 ms |
| `m4` | 8 | 51.5 ms | 41.3 ms |

The `m6` 8-thread run gives a 6.3x commit speedup over the 1-thread run. The
`m4` rows are effectively unchanged because the threshold keeps it on the
sequential path.

## Phase 5: Open Path Optimization

Goal: reduce full proof-generation time after commit is no longer the only
blocker.

### Tasks

1. Measure `open_ms`, `open_fold_ms`, and per-round oracle rebuild time from the
   Phase 0 runner.

2. Reduce cloning in `WhirGrProver::open`:

   - Avoid cloning trees unless a later round really needs an owned tree.
   - Avoid cloning full oracle vectors when the current round can borrow state.
   - Keep proof output owned.

3. Replace per-fiber temporary `Vec` allocations in `fold_eval` and
   `ternary_fold_table` with fixed-size arrays or scratch buffers for ternary
   fibers.

4. Reuse C++-proven direction: avoid full-table work where the proof only needs
   sparse virtual fold query evidence. Keep full-table folding only where it is
   necessary to derive the next committed oracle.

5. Preserve transcript messages and Merkle payload bytes exactly.

### Validation

Run:

```bash
cargo test --lib whir_gr
cargo run --release --bin whir_gr_profile -- --case m4 --phase open --reps 3 --format csv
cargo run --release --bin whir_gr_profile -- --case m5 --phase open --reps 1 --format csv
```

Acceptance:

- Honest openings still verify.
- `serialized_opening_bytes` is unchanged for fixed inputs.
- `open_ms` and allocation count decrease on `m4` and `m5`.

### Phase 5 Review

Status: complete for the measured `m4`/`m5` open bottleneck.

Implemented files:

- `src/protocols/whir_gr/prover.rs`
- `src/protocols/whir_gr/constraint.rs`
- `src/protocols/whir_gr/folding.rs`
- `src/bin/whir_gr_profile.rs`

Measurement changes:

- Added `WhirGrProver::open_profiled` and internal open timing buckets:
  `open_clone_ms`, `open_init_ms`, `open_sumcheck_ms`, `open_restrict_ms`,
  `open_fold_ms`, `open_merkle_open_ms`, `open_constraint_ms`, and
  `open_final_ms`.
- The ordinary `WhirGrProver::open` API still calls the same implementation
  without collecting timings.
- The profile CSV schema was extended with the new open buckets. Existing
  unmeasured fields still print as empty fields rather than fake zeros.

Hotspot evidence before Phase 5 optimization:

| case | open_ms | open_sumcheck_ms | open_fold_ms | encode_oracle_ms | merkle_ms | proof bytes |
|---|---:|---:|---:|---:|---:|---:|
| `m4 open, reps=3` | 1273.196 | 1111.423 | 113.291 | 3.532 | 0.118 | 162788 |
| `m5 open, reps=1` | 5623.550 | 5151.787 | 358.494 | 20.954 | 0.334 | 478260 |

The profile showed that clone time, oracle rebuild, and Merkle rebuild were not
the main causes. The dominant cost was repeated constraint evaluation inside
the honest prover's sumcheck loop, followed by sparse virtual fold query
evaluation.

Optimizations applied:

- Added a local `SumcheckConstraintPlan` for honest sumcheck generation. It
  precomputes the constraint contribution for each interpolation point and
  ternary suffix assignment, then reuses those values inside the `evaluate_f`
  loop.
- Added a focused regression test that compares the precomputed plan against
  direct `WhirConstraint::evaluate_a` evaluation on the same points.
- Added a ternary-specialized `fold_eval` path that avoids the generic
  allocation-heavy interpolation loop for the protocol's three-point fibers.
- Left clone/tree/oracle ownership structure unchanged after measurement showed
  those buckets were small relative to sumcheck and fold cost.

Observed Phase 5 timings:

| case | before Phase 5 | after sumcheck plan | after fold fast path | total speedup |
|---|---:|---:|---:|---:|
| `m4 open, reps=3` | 1273.196 ms | 423.534 ms | 409.055 ms | 3.11x |
| `m5 open, reps=1` | 5623.550 ms | 1743.148 ms | 1703.188 ms | 3.30x |

Final open bucket timings:

| case | open_ms | open_sumcheck_ms | open_fold_ms | open_constraint_ms | encode_oracle_ms | merkle_ms | proof bytes |
|---|---:|---:|---:|---:|---:|---:|---:|
| `m4 open, reps=3` | 409.055 | 262.251 | 98.213 | 30.684 | 3.554 | 0.116 | 162788 |
| `m5 open, reps=1` | 1703.188 | 1274.968 | 314.323 | 69.367 | 21.441 | 0.281 | 478260 |

Protocol-drift review:

- No transcript labels, challenge derivation, soundness parameters, query
  selection policy, proof structs, serializer, verifier logic, or Merkle hash
  function changed.
- `open_profiled` only instruments the same prover path used by `open`.
- A regression test compares `open_profiled` against plain `open` on the same
  input and requires the entire `WhirGrOpening` payload to match.
- The sumcheck plan replaces repeated calls to `constraint.evaluate_a` with
  precomputed values for identical `(interpolation point, suffix assignment)`
  inputs; the direct-evaluation regression test covers this equivalence.
- The ternary fold fast path is the same Lagrange interpolation formula
  specialized to exactly three fiber points.
- Serialized opening bytes remained unchanged for the fixed benchmark inputs:
  `m4=162788`, `m5=478260`.

Validation evidence:

```bash
cargo fmt
cargo test --lib whir_gr
cargo test --lib galois_ring
cargo run --release --bin whir_gr_profile -- --case m4 --phase open --reps 3 --format csv
cargo run --release --bin whir_gr_profile -- --case m5 --phase open --reps 1 --format csv
cargo run --release --bin whir_gr_profile -- --case m4 --phase roundtrip --reps 3 --format csv
cargo run --release --bin whir_gr_profile -- --case m5 --phase verify --reps 3 --format csv
```

The post-optimization `m4 roundtrip` row reported
`commit_ms=51.483`, `open_ms=407.872`, `verify_ms=182.808`,
`serialized_opening_bytes=162788`, and `accepted=true`. The post-optimization
`m5 verify` row reported `verify_ms=455.806`,
`serialized_opening_bytes=478260`, and `accepted=true`.

## Phase 6: Verify Path Optimization

Goal: speed up verifier algebra after proof generation is measurable.

### Tasks

1. Add or isolate a verifier equality-basis cache, analogous to the retained
   C++ `VerifierEqCache` optimization.

2. Cache ternary-grid Lagrange basis coefficients for the verifier.

3. Use the cache in:

   - constraint restriction
   - final polynomial/equality evaluation
   - repeated verifier checks over the same public parameters

4. Keep cache local to `WhirGrVerifier` or a public-parameter-derived helper.
   Avoid hidden process-global state until repeated-construction overhead is
   measured.

### Validation

Run:

```bash
cargo test --lib whir_gr
cargo run --release --bin whir_gr_profile -- --case m4 --phase verify --reps 10 --format csv
cargo run --release --bin whir_gr_profile -- --case m5 --phase verify --reps 3 --format csv
```

Acceptance:

- Verifier accepts honest proofs and rejects malformed proofs covered by current
  tests.
- `verify_ms` decreases.
- Proof bytes and transcript inputs are unchanged.

## Phase 7: Optional Radix-3 / RS Encoding

Goal: consider a larger algorithmic change after the safer wins land.

The C++ prototype retained a hybrid dispatch:

- Use radix-3 `rs_encode` for single-thread WHIR encoding.
- Keep parallel Horner when multi-threaded Horner is faster.

The Rust repo currently has finite-field NTT infrastructure, but not a
Galois-ring radix-3 encoder for `GR(2^s, r)`. Implementing one may become
worthwhile after the multilinear-specialized path and in-place arithmetic are
done.

### Tasks

1. Implement a separate `fft3`/`inverse_fft3` style module for Galois-ring
   domains only after profiling shows dense multiquadratic encoding remains a
   hotspot.

2. Validate against the current `encode_oracle` reference on small powers of 3.

3. Use hybrid dispatch:

   - single-thread: choose radix-3 only if benchmark-proven faster
   - multi-thread: keep chunked Horner/multilinear path if faster

4. Do not replace the multilinear-specialized encoder blindly. Dense
   multiquadratic and multilinear inputs may need different best paths.

### Validation

Run:

```bash
cargo test --lib whir_gr
RAYON_NUM_THREADS=1 cargo run --release --bin whir_gr_profile -- --case m6 --phase commit --reps 1 --format csv
RAYON_NUM_THREADS=8 cargo run --release --bin whir_gr_profile -- --case m6 --phase commit --reps 1 --format csv
```

Acceptance:

- Radix-3 output matches the reference encoder.
- Hybrid dispatch does not regress the best previously measured path.

## Benchmark Discipline

Use this progression after each phase:

1. Unit tests:

   ```bash
   cargo fmt --check
   cargo test --lib whir_gr
   cargo test --lib galois_ring
   ```

2. Small performance check:

   ```bash
   cargo run --release --bin whir_gr_profile -- --case m4 --phase roundtrip --reps 3 --format csv
   cargo run --release --bin whir_gr_profile -- --case m5 --phase commit --reps 3 --format csv
   ```

3. Medium performance check:

   ```bash
   cargo run --release --bin whir_gr_profile -- --case m6 --phase roundtrip --reps 1 --format csv
   ```

4. Large check only after small and medium cases improve:

   ```bash
   timeout 1800s cargo run --release --bin whir_gr_profile -- --case m7 --phase commit --reps 1 --format csv
   timeout 7200s cargo run --release --bin whir_gr_profile -- --case m8 --phase commit --reps 1 --format csv
   ```

Do not run `m9` or `m10` as routine validation. Use them only for final
milestone evidence.

## Acceptance Criteria

The optimization work is successful when all of the following hold:

- `cargo fmt --check` passes.
- `cargo test --lib whir_gr` passes.
- `cargo test --lib galois_ring` passes.
- `m4` and `m5` roundtrip profile rows include commit, open, verify, proof size,
  and `accepted=true`.
- Commitment roots remain stable across optimized and reference paths for small
  deterministic cases.
- Serialized opening bytes remain stable for unchanged proof format.
- `m6 commit` improves by at least 2x after Phase 1 and Phase 3 combined.
- `total alloc` drops substantially after Phase 3.
- Parallel execution is deterministic and improves `m6+` without causing large
  `m4` regressions.

## Risks and Controls

| risk | control |
|---|---|
| Optimized multilinear encoding changes commitment root | Keep old dense encoder as reference in tests for small cases |
| In-place arithmetic introduces aliasing bugs | Document aliasing contract; add tests where `out` aliases neither input first; only support aliasing after explicit tests |
| Parallel encoding becomes nondeterministic | Compare full oracle vectors and roots between 1-thread and 8-thread runs |
| Allocation profiler distorts timing | Separate timing runner from allocation-focused Divan runs |
| Large cases consume a full workday | Gate `m8+` behind explicit timeout commands |
| Proof-size changes silently | Print and test `serialized_opening_bytes` |

## Suggested Implementation Order

1. Phase 0 measurement harness and proof-size reporting.
2. Phase 1 multilinear-specialized commit encoder.
3. Phase 2 sequential domain iteration.
4. Phase 3 in-place Galois-ring arithmetic for commit hot loops.
5. Phase 4 parallel encoding.
6. Phase 5 open path clone/fold allocation cleanup.
7. Phase 6 verifier equality cache.
8. Phase 7 radix-3 encoding only if fresh profiles still justify it.

This order is intentional: reduce algorithmic work first, then remove
allocation pressure, then parallelize. Parallelizing the current allocation-heavy
dense Horner path risks scaling allocator overhead instead of useful work.
