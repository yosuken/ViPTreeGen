
# ViPTreeGen - a standalone tool for viral proteomic tree generation

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE)
[![size](https://img.shields.io/github/size/webcaetano/craft/build/phaser-craft.min.js.svg)]()
[![doi](https://img.shields.io/badge/doi-10.1093%2Fbioinformatics%2Fbtx157-blue.svg?style=flat)](https://doi.org/10.1093/bioinformatics/btx157)

ViPTreeGen is a tool for automated generation of viral "proteomic tree" by computing genome-wide sequence similarities based on tBLASTx results.
The original proteomic tree (i.e., "the Phage Proteomic Tree”) was developed by [Rohwer and Edwards, 2002](https://doi.org/10.1128/JB.184.16.4529-4535.2002).
A proteomic tree is a dendrogram that reveals global genomic similarity relationships between tens, hundreds, or thousands of viruses.
It has been shown that viral groups identified in a proteomic tree well correspond to established viral taxonomies.
The proteomic tree approach is effective to investigate genomes of newly sequenced viruses as well as those identified in metagenomes.

ViPTreeGen has been developed as a part of [the ViPTree server project](http://www.genome.jp/viptree).

## requirements
* BLAST+ (only when `--mode tblastx`, the default)
* MMseqs2 (only when `--mode mmseqs`)
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

Two performance changes, both backward-compatible by default:

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

### 2. `--mode {tblastx|mmseqs}` — alternative search engine

```
ViPTreeGen --mode mmseqs [other options] <input fasta> <output dir>
```

`--mode mmseqs` invokes a single global `mmseqs search --search-type 4` (translated 6-frame, BLOSUM62 default) against `all.fasta`, then `mmseqs convertalis` to BLAST-tab, then splits the result by query into per-node `tblastx.out` files for the existing Rust binary to consume. This is the recommended path for large datasets:

| | tblastx (default) | mmseqs |
|---|---|---|
| Algorithm | NCBI BLAST+ heuristic | MMseqs2 prefilter + Smith-Waterman, translated 6-frame |
| Matrix default | BLOSUM45 | BLOSUM62 (mmseqs2 built-in; BLOSUM45 not bundled) |
| `-dbsize` | applied | ignored (mmseqs uses true target DB size) |
| Parallelism | per-(query,split) via GNU parallel × `--ncpus` | one mmseqs process with `--threads N` |
| EVG 1811 seq | 145 min (12 threads) | **12.5 min (5 threads)** — 11.6x speedup, measured |
| Peak RAM (EVG, defaults) | < 1 GB | ~4 GB (mmseqs target index in memory) |
| Peak RAM cap | (n/a) | `--mmseqs-split-memory-limit` (default `12G`); mmseqs auto-splits the target index if it would exceed the cap |
| Result fidelity vs tblastx | (n/a) | NOT bit-identical; biological signal (SG matrix Pearson) typically > 0.99 |

Result fidelity caveat: mmseqs translated mode and tblastx differ algorithmically. The proteomic-tree topology is preserved (only a few branch swaps in deep clades), but exact SG scores will differ. Choose `--mode mmseqs` when you need throughput on >500 sequences; stick with `--mode tblastx` for maximum compatibility with prior runs.

A `search_mode` row is written to `run.duckdb`'s `run_metadata` table so each output dir self-documents which engine produced it:
```
duckdb path/to/output/run.duckdb "SELECT key, value FROM run_metadata WHERE key='search_mode'"
```

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
