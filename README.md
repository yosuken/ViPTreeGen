
# ViPTreeGen - a standalone tool for viral proteomic tree generation

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE)
[![size](https://img.shields.io/github/size/webcaetano/craft/build/phaser-craft.min.js.svg)]()
[![doi](https://img.shields.io/badge/doi-10.1093%2Fbioinformatics%2Fbtx157-blue.svg?style=flat)](https://doi.org/10.1093/bioinformatics/btx157)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/viptreegen)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/viptreegen/badges/version.svg)](https://anaconda.org/bioconda/viptreegen)

ViPTreeGen is a tool for automated generation of viral "proteomic tree" by computing genome-wide sequence similarities based on tBLASTx results.
The original proteomic tree (i.e., "the Phage Proteomic Tree”) was developed by [Rohwer and Edwards, 2002](https://doi.org/10.1128/JB.184.16.4529-4535.2002).
A proteomic tree is a dendrogram that reveals global genomic similarity relationships between tens, hundreds, or thousands of viruses.
It has been shown that viral groups identified in a proteomic tree well correspond to established viral taxonomies.
The proteomic tree approach is effective to investigate genomes of newly sequenced viruses as well as those identified in metagenomes.

ViPTreeGen has been developed as a part of [the ViPTree server project](http://www.genome.jp/viptree).

## requirements
* MMseqs2 (when `--mode mmseqs-tblastx`, the default)
* BLAST+ (when `--mode tblastx` or `--mode blastn`)
* Ruby (ver >=2.0; tested with 2.x/3.x/4.x)
* R (ver >=3.0)
* DuckDB CLI (ver >=1.0) -- aggregates run state into `${outdir}/run.duckdb`
* `viptreegen-summary-pre` -- bundled Rust binary (new in v2.0.0). Installed automatically by Bioconda; for source builds run `cargo build --release --manifest-path rust/Cargo.toml`

## install

### via Bioconda
```
conda install -c bioconda viptreegen
```

### from source
```
git clone https://github.com/yosuken/ViPTreeGen.git
cd ViPTreeGen
cargo build --release --manifest-path rust/Cargo.toml
./ViPTreeGen --help
```
The bundled Rust binary is located at `rust/target/release/viptreegen-summary-pre`. The `ViPTreeGen` wrapper finds it via PATH, the `VIPTREEGEN_SUMMARY_PRE` env var, or this dev path.

## usage 
```
### ViPTreeGen ver 2.0.0 (2026-05-17) ###

[description]
ViPTreeGen - tool for viral proteomic tree generation from viral genomic sequences.
ViPTreeGen has been developed as the ViPTree server project (http://www.genome.jp/viptree).

ViPTreeGen first computes genome-wide sequence similarity distance based on tBLASTx results,
then construct (bio)nj tree based on the distance (1 - similarity) matrix.

Complete genomes are recommended as input sequence for accurate distance/tree computation,
though genome fragments are also acceptable.

If you compute many sequences (e.g. n > 100) or large sequences (e.g. NCLDV genomes), it may take a long time.
In those cases, use '--ncpus' or '--queue' for parallel computating.

[usage]
$ ViPTreeGen [options] <input fasta> <output dir>

- <input fasta> should be in nucleotide FASTA format and include at least 3 sequences.
- If sequence name (before the first space in the header line) includes a character other than
  alphabets, numbers, dot(.), hyphen(-), or underscore(_), it will be replaced with underscore.
- <output dir> should not exist.

[dependencies]
    - tblastx                 -- included in the BLAST+ program;
                                 https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
    - ruby (ver >=2.0)
    - duckdb (>=1.0)          -- CLI; aggregates run state into ${outdir}/run.duckdb
    - viptreegen-summary-pre  -- bundled Rust binary; new in v2.0.0

  (for tree generation)
    - R (ver >=3.0)
    - R package 'ape'      -- try install.packages('ape') in R terminal
    - R package 'phangorn' -- try install.packages('phangorn') in R terminal

[options]
  (general)
    -h, --help
    -v, --version

  (tblastx)
    --cutlen       [>=10000]   (default: 100000)    -- length to split input sequences for faster tblastx computation
    --dbsize       [int]       (default: 200000000)
    --matrix       [str]       (default: BLOSUM45)
    --evalue       [num]       (default: 1e-2)
    --min-idt      [int]       (default: 30)
    --min-aalen    [int]       (default: 30)

  (tree)
    --notree                   (default: off)       -- generate only similarity/distance matrix
    --method       [nj|bionj]  (default: bionj)     -- proteomic tree generation method

  (2D mode)
    --2D           [query fasta] (default: off)     -- do not generate tree but similarity matrix of: 'query sequences' against 'input sequences'. 2D mode is designed to find the most related 'input sequence' for each 'query sequence'

  (resume)
    --resume                                        -- continue from the last completed pipeline step. <output dir> must
                                                       already exist and contain run.duckdb from the prior run; <input fasta>
                                                       is required syntactically but ignored at runtime (the input is read
                                                       from cat/all/all.fasta in the existing dir). Parameter integrity is
                                                       enforced: --mode/--evalue/--matrix/--ref-duckdb/etc. must match the
                                                       prior run, otherwise --resume is rejected.

  (with-reference mode)
    --ref-duckdb   [path/run.duckdb]                -- a previous ViPTreeGen run.duckdb that supplies the reference set.
                                                       ref sequences, self_scores, and ref-vs-ref summary_tsv are imported
                                                       from it; input is searched against (input ∪ ref) so that input×ref
                                                       pairs are computed but ref×ref pairs are NOT recomputed. Output is
                                                       the combined sim/dist matrix and tree over (ref ∪ input). Requires
                                                       ref schema_version >= 2.0. Mutually exclusive with --2D.

  (use GNU parallel)
    --ncpus        [int]                            -- number of jobs in parallel

  (for icr user)                                    -- for computation in the ICR supercomputer system
    --queue        [JP1]                            -- queue for computation

[output files]
  (normal mode)
    result/all.sim.matrix              -- similarity (SG score) matrix
    result/all.dist.matrix             -- distance (1-SG score) matrix
    result/all.[bio]nj.[a|de]sc.newick -- Newick files of the viral proteomic tree, midpoint rooted and ladderized

  (2D mode)
    result/2D.sim.matrix               -- similarity (SG score) matrix (row-wise: query fasta, column-wise: input fasta)
    result/top10.sim.list              -- top10 SG scores for each query sequence. format: 1st column - query ID; from 2nd to 11st columns - ID:SG_score

  (aggregated state)
    run.duckdb                         -- single DuckDB file containing summary_pre, self_scores, summary_tsv, sequences,
                                          run_metadata, logs. Inspect with `duckdb run.duckdb`.
                                          Per-node tblastx.out files are deleted after the Rust binary consumes them.
```

## what changed in v2.0.0

Three changes — the first two improve performance, the third adds a new analysis mode. All are backward-compatible by default (existing invocations behave unchanged):

### 1. Rust binary for post-tBLASTx pipeline

The legacy four pipeline steps `01-3.cat_and_rename_split_tblastx` + `02-1.tblastx_filter` + `02-2.make_sbed_of_blast` + `02-3.make_summary_pre` are fused into a single Rust binary `viptreegen-summary-pre` that:

- reads per-node `tblastx.out` files directly (split-aware)
- applies the q.start/q.end shift, the pident/alen filter, the BED conversion, and the sweep-line interval merge in one streaming pass
- runs all nodes in parallel via rayon (threads = `--ncpus`)
- emits a single TSV that is COPYed into `run.duckdb`'s `summary_pre` table

Empirical benchmarks (1,811 viral genomes, ~8.7M filtered HSPs, 12 threads):

| step | v1.2.0 (Ruby + SQL) | v2.0.0 (Rust binary) | speedup |
|---|---|---|---|
| 02-3 wall clock | 114 s | 2.4 s | **48x** |
| peak RAM | (per-node) | 280 MB | bounded |

The DuckDB schema lost `blast_hits`, `blast_hits_filtered`, and `sbed_entries` (intermediate hit tables that the legacy SQL pipeline materialized). `summary_pre` and downstream tables are unchanged.

### 2. `--mode {mmseqs-tblastx|tblastx}` — search engine choice (default: mmseqs-tblastx)

```
ViPTreeGen [other options] <input fasta> <output dir>            # uses mmseqs-tblastx, the new default
ViPTreeGen --mode tblastx [other options] <input fasta> <output dir>   # legacy BLAST+ tblastx
```

Both modes compute the same kind of result: translated 6-frame all-vs-all alignment HSPs (i.e., tblastx-style). They differ only in implementation:

- **`mmseqs-tblastx` (default)** invokes one global `mmseqs search --search-type 4` against `all.fasta`, then `mmseqs convertalis` to BLAST-tab, then splits the result by query into per-node `tblastx.out` files for the existing Rust binary to consume. Recommended for all but the smallest datasets.
- **`tblastx`** runs the legacy NCBI BLAST+ `tblastx` binary per (query, split) batch via GNU parallel. Slower; needed for byte-exact reproduction of pre-v2.0 outputs.

| | mmseqs-tblastx (default) | tblastx |
|---|---|---|
| Algorithm | MMseqs2 prefilter + Smith-Waterman, translated 6-frame | NCBI BLAST+ heuristic, translated 6-frame |
| Matrix default | BLOSUM45 (bundled under `data/blosum45.out`, mmseqs2 format) | BLOSUM45 (NCBI BLAST+ default for tblastx is BLOSUM62; ViPTreeGen sets BLOSUM45 via `--matrix`) |
| `--dbsize` | ignored (mmseqs uses true target DB size) | applied |
| Parallelism | one mmseqs process with `--threads N` | per-(query,split) via GNU parallel × `--ncpus` |
| EVG 1811 seq | **12.5 min (5 threads)** — 11.6x speedup, measured | 145 min (12 threads) |
| Peak RAM (EVG, defaults) | ~4 GB (mmseqs target index in memory) | < 1 GB |
| Peak RAM cap | `--mmseqs-split-memory-limit` (default `12G`); mmseqs auto-splits the target index if it would exceed the cap | (n/a) |
| Result fidelity | reference for v2.0+ | NOT bit-identical to mmseqs-tblastx; matches pre-v2.0 byte-for-byte. SG matrix Pearson r typically > 0.99 between the two engines |

Result fidelity caveat: the two engines differ algorithmically (different matrices, heuristics, score normalization), so SG scores do not match bit-for-bit. The proteomic-tree topology is preserved (only a few branch swaps in deep clades). Choose `--mode tblastx` if you need to compare against pre-v2.0 runs; otherwise stay on the default `mmseqs-tblastx`.

### 2b. `--mode blastn` — nucleotide tree (DiGAlign backend)

In addition to the two proteomic-tree modes above, ViPTreeGen v2.0 supports nucleotide-vs-nucleotide search for use as the backend of the **DiGAlign** web tool. The output structure is identical (`result/all.{sim,dist}.matrix` + Newick tree), but the tree is a **nucleotide-identity tree, NOT a proteomic tree** — biological interpretation differs.

```
ViPTreeGen --mode blastn [options] <input.fasta> <output dir>
```

Notes:

- Feeds the same SG-style scoring pipeline (Rust binary → DuckDB → matrix) as the proteomic-tree modes, but at the nucleotide level.
- `--matrix` is **ignored** (blastn uses its `-reward`/`-penalty` defaults).
- `--min-aalen` filter is in **nucleotide units** in this mode (despite the option name); the default `30` is very lenient for nucleotide alignments. Consider raising it (e.g., `--min-aalen 90`) for stricter HSP filtering.
- Use the proteomic-tree modes (`tblastx` / `mmseqs-tblastx`) for **ViPTree** analyses; use this mode for **DiGAlign** analyses.
- An `mmseqs-blastn` mode was evaluated (mmseqs `--search-type 3`) but removed — at our typical input size it was ~9× slower than `blastn`, used ~80× more RAM, and was ~16% less sensitive (Pearson r ≈ 0.70 vs `blastn`). See [`doc/mmseqs-blastn.md`](doc/mmseqs-blastn.md) for the full evaluation.

A `search_mode` row is written to `run.duckdb`'s `run_metadata` table so each output dir self-documents which engine produced it:
```
duckdb path/to/output/run.duckdb "SELECT key, value FROM run_metadata WHERE key='search_mode'"
```

### 3. `--ref-duckdb` — with-reference mode

```
ViPTreeGen --ref-duckdb path/to/ref/run.duckdb [other options] <input.fasta> <output dir>
```

For repeated analyses against a fixed reference panel (e.g., ICTV viruses, an institutional virus catalog, or a previous internal dataset), `--ref-duckdb` reuses a previously computed ViPTreeGen `run.duckdb` so that **only input-vs-input and input-vs-ref pairs are computed**; ref-vs-ref similarities, self-scores, and the ref FASTA itself are taken from the reference's `run.duckdb`. The output is a full combined sim/dist matrix and proteomic tree over `(ref ∪ input)`.

Workflow:
```
# (one-time) build a reference run.duckdb from a large set:
ViPTreeGen large_reference.fasta /data/refset_v1

# (per request) classify a small input against that reference:
ViPTreeGen --ref-duckdb /data/refset_v1/run.duckdb new_input.fasta /tmp/result
```

What gets imported from `--ref-duckdb` (via DuckDB `ATTACH ... (READ_ONLY)`):
- `sequences` (including the nucleotide `seq` column → `all.fasta` is rebuilt without needing the original ref FASTA file)
- `self_scores`
- `summary_tsv` (ref-vs-ref pairs, already bidirectional)

Notes:
- Requires the reference to have been generated with **ViPTreeGen v2.0+** (`schema_version >= 2.0`); older `run.duckdb` files do not include the FASTA. The CLI prints a clear error and exits if the schema is too old.
- Input seq IDs (after normalization) must not collide with ref IDs; the CLI lists conflicts and exits.
- `--ref-duckdb` is mutually exclusive with `--2D`.
- Works with both `--mode mmseqs-tblastx` (default) and `--mode tblastx`.
- SG scores for `input × ref` pairs are not bit-identical to a baseline all-vs-all run because only one HSP direction (input → ref) is computed in this mode, while a baseline run computes both directions. In testing the residual difference is small (Pearson r ≈ 0.9993, max |Δ| ≈ 0.005 on small testdata); ref-vs-ref and input-vs-input cells match exactly.

A `ref_mode=true` and `ref_duckdb_path=...` and `ref_count=N` are recorded in the output `run.duckdb`'s `run_metadata` for provenance.

## inspecting aggregated data
```
duckdb path/to/output/run.duckdb '.tables'
duckdb path/to/output/run.duckdb 'SELECT COUNT(*) FROM summary_pre'
duckdb path/to/output/run.duckdb "SELECT que, sub, mean_idt_of_que FROM summary_pre WHERE que='SEQ_ID'"
duckdb path/to/output/run.duckdb 'SELECT node, target, sg FROM summary_tsv LIMIT 10'
```

## citation
If you use results (data / figures) genereted by ViPTree in your research, please cite:
```
ViPTree: the viral proteomic tree server. Bioinformatics 33:2379–2380 (2017), doi:10.1093/bioinformatics/btx157
Yosuke Nishimura, Takashi Yoshida, Megumi Kuronishi, Hideya Uehara, Hiroyuki Ogata, and Susumu Goto
```
[Author manuscript](http://www.genome.jp/viptree/img/AM_Nishimura_Bioinformatics_2017.pdf) is freely available. 
