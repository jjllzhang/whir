# Anonymous WHIR Review Artifact

This repository is an anonymous review artifact containing a modified WHIR
implementation, benchmark code, result files, and supporting scripts used for
the reported measurements. It is intended for artifact review and measurement
reproduction only; the paper remains self-contained and does not rely on this
repository.

This artifact is derived from the public WHIR implementation. The upstream
copyright and license notices are retained unchanged in `LICENSE-APACHE` and
`LICENSE-MIT`; see `THIRD_PARTY_LICENSES.md` for details.

**WARNING:** This is an academic prototype and has not received careful code
review. This implementation is not ready for production use.

## Contents

- `src/`: protocol and algebra implementation.
- `benches/`: benchmark entry points.
- `results/`: result files used for reported measurements.
- `examples/`: small executable examples.
- `docs/`: supplementary implementation notes.

## Build

```bash
cargo build --release
```

Run tests:

```bash
cargo test
```

Run the finite-field WHIR command-line benchmark/help entry point:

```bash
cargo run --release -- --help
```

Run the Galois-ring WHIR benchmark:

```bash
cargo bench --bench whir_gr --features "parallel asm" -- --ignored
```

For a fixed Rayon thread count, set `RAYON_NUM_THREADS` explicitly:

```bash
RAYON_NUM_THREADS=8 cargo bench --bench whir_gr --features "parallel asm" -- --ignored
```

## Finite-Field WHIR CLI

```text
Usage: main [OPTIONS]

Options:
  -t, --type <PROTOCOL_TYPE>             [default: PCS]
  -l, --security-level <SECURITY_LEVEL>  [default: 100]
  -p, --pow-bits <POW_BITS>
  -d, --num-variables <NUM_VARIABLES>    [default: 20]
  -e, --evaluations <NUM_EVALUATIONS>    [default: 1]
  -r, --rate <RATE>                      [default: 1]
      --reps <VERIFIER_REPETITIONS>      [default: 1000]
  -k, --fold <FOLDING_FACTOR>            [default: 4]
      --sec <SOUNDNESS_TYPE>             [default: ConjectureList]
      --fold_type <FOLD_OPTIMISATION>    [default: ProverHelps]
  -f, --field <FIELD>                    [default: Goldilocks2]
      --hash <MERKLE_TREE>               [default: Blake3]
  -h, --help                             Print help
  -V, --version                          Print version
```

Options:

- `-t` can be either `PCS` or `LDT`.
- `-l` sets the overall security level.
- `-p` sets the number of query-phase proof-of-work bits.
- `-d` sets the number of variables.
- `-e` sets the number of evaluations to prove in PCS mode.
- `-r` sets the inverse-rate logarithm.
- `-k` sets the number of variables folded at each iteration.
- `--sec` selects `UniqueDecoding`, `ProvableList`, or `ConjectureList`.
- `--fold_type` selects `Naive` or `ProverHelps`.
- `-f` selects `Goldilocks2`, `Goldilocks3`, `Field192`, or `Field256`.
- `--hash` selects `SHA3` or `Blake3`.
