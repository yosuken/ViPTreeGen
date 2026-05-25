# Changelog

All notable changes to ViPTreeGen are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and ViPTreeGen
uses [semantic versioning](https://semver.org/). The DuckDB `run.duckdb` schema is versioned
independently via the `schema_version` row in `run_metadata` and bumps follow the same
breaking-vs-additive rule (a `schema_version` major bump means previously-generated
`run.duckdb` files are not accepted by the new release, e.g. as `--ref-duckdb` input).

## [Unreleased] — v2.0.0

Major performance + capability overhaul. Backward-compatible CLI by default, but the
default search engine, schema, and on-disk output structure all change.

### Breaking
- **DuckDB schema bumped to 2.0** (`schema_version` in `run_metadata`). The `sequences`
  table gains a `seq VARCHAR NOT NULL` column holding the actual nucleotide sequence,
  so a `run.duckdb` is now self-sufficient as a reference set without the original
  FASTA alongside. `run.duckdb` produced by v1.x is *not* accepted by `--ref-duckdb`.
- **Intermediate tables removed.** `blast_hits`, `blast_hits_filtered`, and
  `sbed_entries` are gone — the Rust binary streams that work in memory and writes
  only the final `summary_pre`. Anything that queried those tables must be ported.
- **Default `--mode` changed from `tblastx` to `mmseqs-tblastx`**. Result matrices are
  no longer bit-identical to v1.x by default. Pass `--mode tblastx` to reproduce
  pre-v2.0 numbers byte-for-byte.

### Added
- **Tool provenance in `run.duckdb`** — a new `tools` table records the resolved
  path and reported version of every external tool the run actually used
  (engine binaries, `duckdb`, `parallel`, `R`, `ruby`, and the bundled
  `viptreegen-summary-pre`), selected mode-aware. `run_metadata` additionally
  stores the full `command_line` and every effective option (idt, aalen, ncpus,
  notree, resume, mmseqs split/chunk sizes, 2D query path), so a
  finished `run.duckdb` fully captures how it was produced.
  `viptreegen-summary-pre` gained a `--version` flag for this.
- **Gzip-compressed FASTA input** — `<input fasta>` and the `--2D` query file may
  be gzipped. Detection is by the gzip magic bytes (`1f 8b`), so it works
  regardless of extension (`.gz`, `.fasta.gz`, or none). No flag needed.
- **Benchmark subsets under `testdata/benchmark/`** — four gzipped Virus-Host DB
  taxon subsets (Begomovirus, Autographivirales, Orthoherpesviridae, Imitervirales)
  spanning the genome-size regimes from `doc/engine-benchmark.md`. ~5 MB total;
  see `testdata/benchmark/README.md` for provenance and sampling.
- **`--mode {mmseqs-tblastx|tblastx|blastn|last}`** — search engine + algorithm choice.
  - `mmseqs-tblastx` (default): MMseqs2 `--search-type 4`, translated 6-frame protein
    alignment. ~10× faster than legacy tblastx on EVG (1811 viral genomes:
    12.5 min vs 145 min). Proteomic tree.
  - `tblastx`: legacy NCBI BLAST+ tblastx. Slow, byte-exact pre-v2.0 reproduction.
  - `blastn`: NCBI BLAST+ nucleotide-vs-nucleotide. Backend for the DiGAlign web tool —
    output is a *nucleotide tree*, NOT a proteomic tree.
  - `last`: LAST (`lastdb` + `lastal -fBlastTab`), nucleotide-vs-nucleotide; DiGAlign
    backend. BlastTab output is byte-compatible with `blastn -outfmt 6` (no converter).
    Much faster + far lower memory than blastn on **many small** genomes (870 × 2.7 kb:
    20 s / 21 MB vs blastn 154 s / 1.7 GB); blastn is faster on **few large** genomes
    (23 × 1.2 Mb: blastn 123 s vs LAST 355 s). Matrices correlate Pearson 0.86–0.99 with
    blastn; cluster membership ARI 0.81–0.88. No 2D support. Single global `lastal` call.
  - `mmseqs-blastn` (mmseqs `--search-type 3`) and MUMmer4 `nucmer` were evaluated and
    NOT adopted; see `doc/mmseqs-blastn.md`, `doc/mummer4-nucmer.md`, and the
    three-dataset blastn-vs-LAST comparison in `doc/engine-benchmark.md`.
- **`--resume`** — resume a previously-interrupted run from the last completed
  pipeline step. Each Rake task records `step_done:<task.name>` in
  `run_metadata` when it finishes; `--resume` reads those markers and skips
  completed steps. Parameter integrity (mode, matrix, evalue, schema_version,
  ref-duckdb path, …) is verified against the prior run and `--resume` is
  rejected on any mismatch. Useful when a long `mmseqs search` or `tblastx`
  step is killed by OOM / a time limit -- restart with the same arguments
  plus `--resume` and the pipeline continues from where it stopped.
  In `--mode tblastx` / `--mode blastn`, `--resume` is also batch-level: 01-2
  filters out batch lines whose `-out FILE` already exists, so partially-
  completed per-(query, split) batches survive across restarts. To make this
  safe, SearchCmd writes the result to `<FILE>.tmp` and `mv`s to `<FILE>`
  only on a zero exit -- a tblastx/blastn killed mid-write leaves a `.tmp`
  but no `<FILE>`, so the batch is re-run cleanly on resume.
- **`--mmseqs-tblastx-chunk-size INT`** (default 1000) — in `--mode
  mmseqs-tblastx`, split the query fasta into N-sequence chunks and run
  `mmseqs search` per chunk. Each chunk gets its own
  `step_done:01-2.chunk:NNNN` marker so `--resume` can pick up between
  chunks instead of restarting the whole search. The target DB is built
  once and shared across all chunks (no extra index-build overhead). At
  the default chunk size, small datasets (≤1000 seqs) run as a single
  chunk identical to the old behavior, and large datasets (e.g. ICTV
  ~6000 seqs) become a 6-chunk pipeline where a killed chunk loses at
  most a few minutes of work. Pass `0` to disable chunking (one global
  search). 2D mode keeps its existing two-search flow (not chunked).
  Secondary benefit: because each chunk is a separate mmseqs process, a
  smaller chunk size caps per-process CPU time and RAM. This lets a large
  translated search complete under an interactive node's `ulimit -t`
  (e.g. 30 min CPU) without a batch job — verified on 1780 T7-like phage
  genomes (74 MB) where the default 2-chunk run was killed at the 30-min
  CPU limit but `--mmseqs-tblastx-chunk-size 200` (9 chunks) completed.
- **`--ref-duckdb PATH`** — with-reference mode. Reuse a previous v2.0 `run.duckdb`
  as the reference set: ref sequences, `self_scores`, and ref-vs-ref `summary_tsv`
  are ATTACHed and copied; only input-vs-input and input-vs-ref are recomputed. The
  output matrix and tree cover `(ref ∪ input)`. Designed for ViPTree v2 web backend.
- **Bundled Rust binary `viptreegen-summary-pre`** — fuses the legacy 01-3 (cat+rename),
  02-1 (filter), 02-2 (BED conversion), and 02-3 (interval merge) into one parallel
  single-process step. EVG 02-3 wall clock: 114 s → 2.4 s (48× speedup). Built with
  `cargo build --release --manifest-path rust/Cargo.toml` (Bioconda installs it
  automatically).
- **`--mmseqs-split-memory-limit SIZE`** (default `12G`). Caps peak RAM used by
  the mmseqs target index. Pass `0` to let mmseqs use as much RAM as needed.
- **Bundled BLOSUM45 substitution matrix** under `data/blosum45.out` (mmseqs2
  format, freely redistributable), passed to mmseqs via `--sub-mat`. Without this,
  mmseqs-tblastx would default to BLOSUM62 and disagree more with `--mode tblastx`.
- **Ruby CLI entry point** (`./viptreegen`) replaces the previous bash wrapper.
  OptionParser-based argument handling, mode-aware dependency check, structured
  logging via `Open3`.
- **Per-step logs ingested into `run.duckdb`'s `logs` table** (with `source` like
  `tblastx`, `mmseqs:search.node`, `makeblastdb`). Per-node `tblastx.out` files are
  deleted after the Rust binary consumes them — `run.duckdb` is now the single
  authoritative output artifact.
- **Bioconda package**: `conda install -c bioconda viptreegen`.
- **CI smoke tests** (`.github/workflows/test.yml` + `test.sh`) covering all three
  modes, normal/2D/with-reference flows, validation rejects, golden matrix
  comparison for `--mode tblastx`, and tree generation.

### Changed
- **Primary command renamed to lowercase `viptreegen`** (matching the conda package
  name). The previous `ViPTreeGen` is kept as a backward-compatible alias (a guarded
  symlink in the Bioconda package); from a source clone the command is `./viptreegen`.
  The displayed brand name remains "ViPTreeGen".
- **Minimum Ruby bumped to 3.0** (was 2.0). The CLI now rejects Ruby < 3.0.
- **FASTA label normalization** (replace non-`[A-Za-z0-9.\-_]` with `_`, strip
  leading/trailing dots, collapse consecutive separators) is now factored into
  a single `ParseFastaEntries` lambda used by 01-1 / 01-1.2D / 01-1.ref.prep.
- **Rust binary rounding** is half-to-even (banker's rounding, IEEE 754), not
  half-away-from-zero. Affects `%.1f`-formatted `pct_que_len_in_hit` and friends
  at exact-half boundaries (rare in practice).

### Removed
- **`--queue` / `--wtime`** (ICR supercomputer `qsub` batch-submission path) removed.
  Parallelism is now via `--ncpus` only (GNU parallel for tblastx/blastn, internal
  engine threads for mmseqs-tblastx/last). For a scheduler, wrap the whole `ViPTreeGen`
  invocation in your own job script.
- `--mode mmseqs-blastn` was evaluated and dropped. mmseqs nucleotide search
  (`--search-type 3`) at our typical input size was ~9× slower than `blastn`,
  used ~80× more RAM, and was ~16% less sensitive (Pearson r ≈ 0.70 vs blastn).
  See `doc/mmseqs-blastn.md` for the full evaluation.

---

## [1.1.4] — 2026-05-16

### Added
- Bioconda packaging: `meta.yaml` and `environment.yaml` shipped with the repo, with
  `run_exports` set for Bioconda lint compliance.

### Fixed
- Ruby 4.x compatibility (`Open3` / `OptionParser` interactions tightened up).

---

## [1.1.3]

### Fixed
- Misc. small fixes carried from pre-Bioconda development; no schema changes.

---

## [1.1.2] and earlier

See `git log` for the granular history of the v1.x line — `run.duckdb` aggregation
landed in v1.2.0 (unreleased; folded into v2.0), prior to which all intermediate
state lived as files under the output directory.
