# WHIR-GR Oracle Encoding Optimization Plan

Date: 2026-05-07

Target repository: `/home/zjl/whir`

Current baseline commit:

```text
558f7c0 perf: parallelize whir-gr hot paths
```

This plan covers one narrow optimization line: replacing pointwise oracle
evaluation in WHIR-GR with a protocol-equivalent structured encoding path.

The optimization must be rejected if it does not deliver at least `1.5x`
end-to-end speedup on the accepted m10 benchmark gate.

## Decision Summary

The next large speedup opportunity is real, but only if it optimizes both:

- initial commit oracle encoding
- prover `open` oracle rebuild

Optimizing only `open` or only `commit` is unlikely to pass the `1.5x`
end-to-end gate.

The concrete target is a Reed-Solomon-style domain evaluation routine for
Galois-ring Teichmuller cosets:

```text
coefficients of P(X) -> [P(offset * root^i)]_i
```

This must produce exactly the same oracle vector as the current Horner/direct
evaluation path.

The optimization is allowed to change only how oracle values are computed. It
must not change protocol parameters, transcript contents, Merkle payloads,
proof structure, query selection, soundness parameters, or serialization.

## Execution Result

Status: accepted.

The implementation completed this plan with a conservative large-domain
dispatch threshold. The structured encoder is enabled only when:

```text
domain.size() >= 10000
coefficients.len() * domain.size() >= 16384
```

This keeps m6/m7 on the previous folded/Horner paths while allowing the large
m10 commit and open oracle rebuilds to use structured encoding.

Final m10 acceptance result with `RAYON_NUM_THREADS=8`:

| phase | commit | open | verify | total |
|---|---:|---:|---:|---:|
| baseline roundtrip | `125.721 s` | `190.528 s` | `3.644 s` | `319.893 s` |
| final roundtrip | `46.918 s` | `84.670 s` | `3.638 s` | `135.225 s` |

Speedup:

```text
319.893 / 135.225 = 2.366x
```

The final result is below the required `213.262 s` threshold, so the
optimization is accepted under the rejection rule.

Final protocol gates:

- `m10 accepted=true`
- `m10 serialized_opening_bytes == 13732680`
- `m10 RAYON_NUM_THREADS=32 roundtrip == 83.270 s`, accepted with unchanged bytes
- `m8 serialized_opening_bytes == 10522264`
- `m7 serialized_opening_bytes == 4908200`
- no changes to `transcript.rs`, `serialization.rs`, or `merkle.rs`
- `cargo test --lib whir_gr` passed
- `cargo test --lib galois_ring` passed
- `cargo clippy --all-targets --all-features --locked -- -D warnings` passed
- `cargo fmt --check` and `git diff --check` passed

Final benchmark artifacts:

```text
results/oracle_encode_trials/m10_roundtrip_threads8_threshold10k.csv
results/oracle_encode_trials/m10_roundtrip_threads32_threshold10k.csv
results/oracle_encode_trials/m10_open_threads8_threshold10k.csv
results/oracle_encode_trials/m10_commit_threads8_threshold10k.csv
results/oracle_encode_trials/m8_roundtrip_threads8_threshold10k.csv
results/oracle_encode_trials/m7_roundtrip_threads8_threshold10k.csv
results/oracle_encode_trials/m6_commit_threads8_threshold10k.csv
```

An earlier `domain.size() >= 1024` dispatch trial was rejected because it made
m6/m7 use structured encoding too early. The retained implementation uses the
`10000` threshold above.

## Current Baseline

Current m10 data with `RAYON_NUM_THREADS=8`:

| phase | current time | important buckets |
|---|---:|---|
| roundtrip | `319.893 s` | commit `125.721 s`, open `190.528 s`, verify `3.644 s` |
| open | `191.165 s` | encode oracle `129.649 s`, sumcheck `57.519 s`, fold `2.830 s` |
| open sumcheck | `57.519 s` | constraint plan `47.959 s`, poly eval `7.359 s` |
| verify | `3.623 s` | fold `2.800 s`, constraint `0.638 s` |

Proof bytes:

```text
m10 serialized_opening_bytes = 13732680
accepted = true
```

The current `open` oracle rebuild path is:

```text
WhirGrProver::open
  prove_round
    encode_oracle
      encode_oracle_parallel
        encode_oracle_chunk
          for point in domain.iter_elements_from(...)
            polynomial.evaluate_pow(ctx, &point)
```

The current multi-quadratic point evaluator is Horner over univariate power
coefficients:

```rust
for coefficient in self.coefficients.iter().rev() {
    acc = acc * x + coefficient;
}
```

For m10, the first `open` rebuild evaluates a degree-bounded polynomial with up
to `3^7 = 2187` coefficients over a domain of size `41391`. That is roughly:

```text
41391 * 2187 ~= 90.5M Galois-ring Horner steps
```

This explains why `open encode_oracle_ms` dominates.

The current initial commit encoder is separate:

```text
commit_multilinear
  encode_multilinear_oracle
    encode_multilinear_oracle_parallel
      encode_multilinear_oracle_chunk
        pow_m_into
        evaluate_multilinear_folded
```

For m10, this is roughly:

```text
124173 * 2^10 ~= 127M folded multilinear steps
```

So commit also needs to be covered if the goal is `>= 1.5x` end-to-end.

## Required Speedup Gate

The accepted baseline is:

```text
m10 roundtrip, RAYON_NUM_THREADS=8: 319.893 s
```

The optimization is accepted only if:

```text
new_m10_roundtrip_time <= 319.893 / 1.5 = 213.262 s
```

In other words:

```text
required m10 roundtrip speedup >= 1.5x
```

This gate is end-to-end. A local bucket win is not enough.

Additional required gates:

- `accepted=true`
- `serialized_opening_bytes == 13732680` for m10
- proof bytes stable across `open`, `verify`, and `roundtrip`
- no transcript label changes
- no Merkle payload or serialization changes
- no soundness parameter or query policy changes
- no more than `10%` regression on m6/m7 roundtrip

If any protocol gate fails, reject immediately regardless of speed.

If protocol gates pass but m10 end-to-end speedup is `< 1.5x`, reject the
optimization and do not commit it.

## Why This Can Be Reasonable

Current oracle encoding is pointwise evaluation:

```text
for each domain point x_i:
    compute P(x_i)
```

For a multiplicative Teichmuller coset:

```text
x_i = offset * root^i
```

the full oracle is a structured Reed-Solomon evaluation. A transform-style
encoder can reuse the domain structure across all points instead of recomputing
the full Horner chain per point.

This is protocol-equivalent because the committed oracle is still:

```text
[P(offset * root^i)]_i
```

Only the algorithm used to compute the same vector changes.

There is also precedent in the C++ WHIR-over-GR path: the retained C++ speedup
used `rs_encode(domain, polynomial.to_univariate_pow_polynomial(ctx))` for
single-thread WHIR encode, with tests comparing it against direct
`evaluate_pow`. That precedent supports the idea, but the Rust implementation
must still be independently validated against current Rust oracle semantics.

## What Is Not Enough

### Incremental Point Generation

Updating points by:

```text
x_{i+1} = x_i * root
```

is already effectively what `DomainElements` does. This is not the main
bottleneck.

Incrementally maintaining powers may save some point-preparation work, but it
does not remove the dominant `domain_size * coefficient_count` evaluation
cost. It is unlikely to pass the `1.5x` gate by itself.

### Generic Cross-Round Oracle Caching

Do not assume that a folded table from the current oracle is the next committed
oracle. In the current protocol flow:

```rust
let next_polynomial = sumcheck_polynomial;
let next_domain = state.domain.pow_map(3)?;
let next_oracle = encode_oracle(params, &next_domain, &next_polynomial)?;
```

The next oracle commits to the restricted next polynomial over `pow_map(3)`.
This is not automatically identical to the sparse shift-fold values or to a
full `pow_map(3^width)` fold table.

Caching is allowed only if there is a proof-backed equality to the direct
`encode_oracle` vector and tests check exact equality. Treat generic caching as
deferred, not as the first implementation.

### Proof-Format or Parameter Changes

Do not optimize by changing:

- `WhirGrPublicParameters`
- `WhirGrProof`
- `WhirGrRoundProof`
- transcript labels
- challenge sampling
- shift/final query positions
- Merkle tree layout
- serialized opening format
- soundness selector output

Those would be protocol changes, not this optimization.

## Expected Impact Model

Current m10 8-thread total:

```text
roundtrip ~= 319.9 s
commit oracle bucket ~= 125.5 s
open oracle bucket ~= 129.6 s
non-oracle total ~= 64.8 s
```

If both commit and open oracle encoding improve by factor `S`, then roughly:

```text
new roundtrip ~= 64.8 + (125.5 + 129.6) / S
```

| combined oracle speedup | estimated roundtrip | end-to-end speedup | accept? |
|---:|---:|---:|---|
| `1.5x` | `~234.9 s` | `~1.36x` | no |
| `2.0x` | `~192.4 s` | `~1.66x` | yes |
| `3.0x` | `~149.8 s` | `~2.14x` | yes |
| `4.0x` | `~128.6 s` | `~2.49x` | yes |

Therefore the practical internal target is:

```text
combined commit+open oracle encoding speedup >= 2x
```

The external acceptance target remains:

```text
m10 roundtrip speedup >= 1.5x
```

If only `open encode_oracle_ms` improves by `2x`, the new roundtrip is roughly:

```text
319.9 - 129.6 / 2 ~= 255.1 s
```

That is only about `1.25x`, so it must be rejected under this plan.

## Implementation Strategy

Work on a separate branch:

```bash
git switch -c perf/whir-gr-rs-encode
```

Do not commit intermediate experiments unless the final acceptance gate passes.

### Phase 0: Baseline Capture

Goal: freeze the baseline used by the rejection gate.

Commands:

```bash
mkdir -p results/oracle_encode_baseline

RAYON_NUM_THREADS=8 cargo run --release --bin whir_gr_profile -- --case m7 --phase roundtrip --reps 1 --format csv | tee results/oracle_encode_baseline/m7_roundtrip_threads8.csv
RAYON_NUM_THREADS=8 timeout 3600s cargo run --release --bin whir_gr_profile -- --case m10 --phase roundtrip --reps 1 --format csv | tee results/oracle_encode_baseline/m10_roundtrip_threads8.csv
RAYON_NUM_THREADS=8 timeout 3600s cargo run --release --bin whir_gr_profile -- --case m10 --phase open --reps 1 --format csv | tee results/oracle_encode_baseline/m10_open_threads8.csv
RAYON_NUM_THREADS=8 timeout 900s cargo run --release --bin whir_gr_profile -- --case m10 --phase verify --reps 1 --format csv | tee results/oracle_encode_baseline/m10_verify_threads8.csv
```

Acceptance:

- Baseline `m10 roundtrip` is close to `319.893 s`.
- Baseline m10 proof bytes are `13732680`.
- `accepted=true`.

If the baseline differs by more than `15%`, recompute the acceptance threshold
from the fresh baseline and write the exact threshold into the experiment notes.

### Phase 1: Reference Oracle Equality Harness

Goal: add tests and internal helpers before changing behavior.

Tasks:

1. Keep the current pointwise encoder as the reference:

   ```rust
   encode_oracle_horner_reference(params, domain, polynomial)
   encode_multilinear_oracle_reference(params, domain, polynomial)
   ```

2. Add focused equality tests:

   - `multi_quadratic_oracle_encoding_should_match_horner_reference`
   - `multilinear_oracle_encoding_should_match_existing_reference`
   - domain sizes covering:
     - pure `3^k`
     - current bench-like `3^a * 511`
     - small cosets with non-one offset

3. Add exact Merkle-root equality tests after building trees from both oracle
   vectors.

Acceptance:

```bash
cargo test --lib whir_gr::prover
cargo test --lib whir_gr::multiquadratic
cargo clippy --all-targets --all-features --locked -- -D warnings
```

No behavior change is allowed in this phase.

### Phase 2: Domain Factorization and Transform API

Goal: add a transform API that can decide whether a domain is supported.

Proposed internal API:

```rust
pub(crate) fn rs_encode_teichmuller_coset(
    ctx: &GrContext,
    domain: &Domain,
    coefficients: &[GrElem],
) -> Result<Option<Vec<GrElem>>>
```

Return semantics:

- `Ok(Some(oracle))`: transform path supports the domain and succeeded.
- `Ok(None)`: unsupported domain shape; caller must use Horner fallback.
- `Err(_)`: invalid polynomial/domain/ring state.

The function must not panic in production code.

Tasks:

1. Factor `domain.size()` into small radices.
2. Support the benchmark sizes first:

   ```text
   124173 = 3^5 * 511
   41391  = 3^4 * 511
   13797  = 3^3 * 511
   4599   = 3^2 * 511
   1533   = 3   * 511
   511    = 7 * 73
   ```

3. Implement a conservative mixed-radix plan:

   ```text
   [3, 3, 3, 3, 3, 7, 73]
   ```

   or the corresponding suffix for smaller domains.

4. Add a memory guard. The transform may allocate one or two full oracle-sized
   buffers, but must not accidentally materialize many full copies.

5. Keep output order exactly equal to `domain.iter_elements()`.

Acceptance:

- plan construction tests for all m4..m10 domain sizes
- unsupported domains fall back cleanly
- no protocol code calls this path yet

### Phase 3: Sequential RS Encode Prototype

Goal: implement correctness-first transform encoding.

Implementation notes:

- Pad coefficients with zero to `domain.size()`.
- Evaluate on `offset * root^i`, not only on the subgroup.
- Handle coset offset either by pre-scaling coefficients:

  ```text
  c_j <- c_j * offset^j
  ```

  then evaluating on powers of `root`, or by an equivalent transform formula.

- Use mixed-radix Cooley-Tukey stages.
- For radix `3`, write a specialized kernel.
- For radix `7` and `73`, a naive local DFT kernel is acceptable initially if
  it is still much cheaper than full pointwise Horner.
- Use scratch-backed `GrContext::mul_into` in inner kernels.
- Do not introduce `unsafe`.

Correctness tests:

```bash
cargo test --lib whir_gr::prover::tests::rs_encode
cargo test --lib whir_gr::multiquadratic
cargo test --lib galois_ring
```

Required equality checks:

- `rs_encode_teichmuller_coset(...) == encode_oracle_horner_reference(...)`
- equality for dense multi-quadratic polynomials
- equality for sparse multilinear-embedded polynomials
- equality over nontrivial coset offset
- equality for each current bench case m4..m10 at reduced or sampled sizes

No integration into commit/open yet.

### Phase 4: Integrate Open Oracle Rebuild Behind Fallback

Goal: use transform encoding for `encode_oracle` only when supported.

Integration logic:

```rust
if rs_encode is supported and estimated worthwhile {
    use rs_encode
} else {
    use current parallel Horner
}
```

Do not delete the current Horner path.

Threshold rules:

- Enable only for `domain.size() >= 1024`.
- Enable only when `coefficients.len() * domain.size()` is large enough to
  dominate transform overhead.
- Keep a debug/test-only equality assertion option for small cases.

Protocol checks:

- `m7` proof bytes unchanged: `4908200`
- `m8` proof bytes unchanged: `10522264`
- `m10` proof bytes unchanged: `13732680`
- `accepted=true`

Benchmark commands:

```bash
RAYON_NUM_THREADS=8 cargo run --release --bin whir_gr_profile -- --case m7 --phase roundtrip --reps 1 --format csv
RAYON_NUM_THREADS=8 timeout 900s cargo run --release --bin whir_gr_profile -- --case m8 --phase open --reps 1 --format csv
RAYON_NUM_THREADS=8 timeout 3600s cargo run --release --bin whir_gr_profile -- --case m10 --phase open --reps 1 --format csv
```

Expected outcome:

- `m10 open encode_oracle_ms` should fall by at least `2x`.
- End-to-end m10 roundtrip may still fail the `1.5x` gate if commit is not
  optimized.

If this phase improves open but the total projected roundtrip remains below
`1.5x`, continue to Phase 5. Do not accept the optimization yet.

### Phase 5: Integrate Commit Oracle Encoding

Goal: make `commit_multilinear` benefit from the same transform.

Approach:

1. Convert the multilinear polynomial to its multi-quadratic embedding.
2. Encode the embedded univariate power coefficients using the same
   `rs_encode_teichmuller_coset`.
3. Compare against the existing `encode_multilinear_oracle_parallel` reference.
4. Use the transform only when it wins under the benchmark thresholds.

Important caveat:

Current multilinear commit encoder evaluates the multilinear polynomial without
materializing all `3^m` coefficients in the hot loop. A transform path must be
benchmarked because sparse embedded coefficients and large transform buffers may
change the tradeoff.

Required tests:

```bash
cargo test --lib whir_gr::prover::tests::parallel_multilinear_oracle_encoding_should_match_sequential
cargo test --lib whir_gr::prover::tests::multilinear_oracle_encoding_should_match_dense_embedding
```

Add one direct test:

```text
rs_multilinear_oracle_encoding_should_match_current_encoder
```

Benchmark commands:

```bash
RAYON_NUM_THREADS=8 cargo run --release --bin whir_gr_profile -- --case m8 --phase commit --reps 1 --format csv
RAYON_NUM_THREADS=8 timeout 1200s cargo run --release --bin whir_gr_profile -- --case m10 --phase commit --reps 1 --format csv
RAYON_NUM_THREADS=8 timeout 3600s cargo run --release --bin whir_gr_profile -- --case m10 --phase roundtrip --reps 1 --format csv
```

Expected outcome:

- combined `commit encode_oracle_ms + open encode_oracle_ms` speedup should be
  at least `2x`
- m10 roundtrip should fall below `213.262 s`

If m10 roundtrip remains above `213.262 s`, reject the optimization.

### Phase 6: Parallelize Transform Stages

Goal: recover scaling when the sequential transform is correct but not enough.

Only start this phase if:

- Phase 3 correctness is complete
- Phase 4 or 5 shows the transform is correct but not yet fast enough

Parallel-safe boundaries:

- Stage-local independent butterflies may run in Rayon.
- Output order must be deterministic.
- Each worker owns scratch buffers.
- No global mutable cache.

Benchmark matrix:

```bash
RAYON_NUM_THREADS=1 cargo run --release --bin whir_gr_profile -- --case m8 --phase open --reps 1 --format csv
RAYON_NUM_THREADS=8 cargo run --release --bin whir_gr_profile -- --case m8 --phase open --reps 1 --format csv
RAYON_NUM_THREADS=32 cargo run --release --bin whir_gr_profile -- --case m8 --phase open --reps 1 --format csv

RAYON_NUM_THREADS=8 timeout 3600s cargo run --release --bin whir_gr_profile -- --case m10 --phase roundtrip --reps 1 --format csv
RAYON_NUM_THREADS=32 timeout 3600s cargo run --release --bin whir_gr_profile -- --case m10 --phase roundtrip --reps 1 --format csv
```

Acceptance:

- 8-thread m10 roundtrip still passes `>=1.5x`.
- 32-thread run is informational, not the primary gate.

### Phase 7: Final Protocol-Drift Review

Before accepting:

```bash
git diff -- src/protocols/whir_gr/serialization.rs
git diff -- src/protocols/whir_gr/transcript.rs
git diff -- src/protocols/whir_gr/merkle.rs
```

These diffs must be empty unless there is a separate explicit protocol review.

Check transcript-sensitive strings:

```bash
rg -n "whir\\.|challenge|derive_unique_positions|serialize_|Merkle|opening" src/protocols/whir_gr
```

Review requirements:

- `opening_transcript(...)` unchanged
- sumcheck polynomial absorption order unchanged
- alpha/gamma challenge labels unchanged
- `whir.g_root` absorption unchanged
- shift query derivation unchanged
- final query derivation unchanged
- Merkle tree built over byte-identical leaves
- opening serialization byte-identical

Proof-byte requirements:

```text
m7  serialized_opening_bytes == 4908200
m8  serialized_opening_bytes == 10522264
m10 serialized_opening_bytes == 13732680
```

## Final Acceptance Commands

Run all:

```bash
cargo fmt --check
cargo test --lib whir_gr
cargo test --lib galois_ring
cargo clippy --all-targets --all-features --locked -- -D warnings
git diff --check

RAYON_NUM_THREADS=8 cargo run --release --bin whir_gr_profile -- --case m7 --phase roundtrip --reps 1 --format csv
RAYON_NUM_THREADS=8 timeout 1800s cargo run --release --bin whir_gr_profile -- --case m8 --phase roundtrip --reps 1 --format csv
RAYON_NUM_THREADS=8 timeout 3600s cargo run --release --bin whir_gr_profile -- --case m10 --phase roundtrip --reps 1 --format csv
RAYON_NUM_THREADS=8 timeout 3600s cargo run --release --bin whir_gr_profile -- --case m10 --phase open --reps 1 --format csv
RAYON_NUM_THREADS=8 timeout 1200s cargo run --release --bin whir_gr_profile -- --case m10 --phase commit --reps 1 --format csv
```

Accept only if:

```text
m10 roundtrip <= 213.262 s
m10 speedup >= 1.5x
m10 serialized_opening_bytes == 13732680
m10 accepted == true
```

Reject if:

```text
m10 roundtrip > 213.262 s
```

even if individual oracle buckets improved.

## Rejection Procedure

If the optimization fails protocol gates:

1. Stop immediately.
2. Do not tune thresholds to hide the failure.
3. Revert the transform integration and any protocol-adjacent changes.
4. Keep only standalone tests if they are useful and do not affect production
   behavior.

If the optimization is correct but speedup is `<1.5x`:

1. Reject the optimization.
2. Do not commit the performance path.
3. Record the measured results in `results/` only.
4. Return to baseline branch:

   ```bash
   git switch feature/whir-gr
   ```

5. Delete the experiment branch only after saving the measured evidence:

   ```bash
   git branch -D perf/whir-gr-rs-encode
   ```

If the optimization is already committed locally before the final gate and then
fails, revert it with a normal revert commit rather than rewriting unrelated
history:

```bash
git revert <optimization-commit>
```

## Commit Rule

Only one final commit is allowed for this optimization line, and only after all
acceptance gates pass.

Suggested commit message:

```text
perf: use structured oracle encoding for whir-gr
```

The commit must include:

- transform implementation
- fallback path
- oracle equality tests
- end-to-end tests
- benchmark evidence in the commit message body or an accompanying tracked
  benchmark note if the repo policy allows it

Do not commit raw `results/` unless explicitly requested.

## Expected Final Outcome

Reasonable target:

```text
m10 roundtrip: 319.9 s -> 130-210 s
speedup: 1.5x-2.5x
```

Conservative pass line:

```text
<= 213.262 s
```

Best case if both commit and open oracle encoding get strong transform wins:

```text
~130-150 s
```

If the final result is worse than `213.262 s`, this plan requires rejecting the
optimization, even if the implementation is correct.
