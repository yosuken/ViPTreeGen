
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

## install

### via Bioconda (recommended)
```
conda install -c bioconda viptreegen
```
or
```
mamba install -c bioconda viptreegen
```

### from source
```
git clone https://github.com/yosuken/ViPTreeGen.git
cd ViPTreeGen
./ViPTreeGen --help
```

## requirements
* MMseqs2 (when `--mode mmseqs-tblastx`, the default)
* BLAST+ (when `--mode tblastx` or `--mode blastn`)
* LAST (when `--mode last`)
* Ruby (ver >=3.0; tested with 3.x/4.x)
* R (ver >=3.0)
* DuckDB CLI (ver >=1.0) -- aggregates run state into `${outdir}/run.duckdb`
* `viptreegen-summary-pre` -- bundled Rust binary (new in v2.0.0). Installed automatically by Bioconda; for source builds run `cargo build --release --manifest-path rust/Cargo.toml`

## install

### via Bioconda
```
conda install -c bioconda viptreegen
```
The command is installed as `ViPTreeGen`; a lowercase `viptreegen` alias (matching the package name) is also provided, so either works.

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
In those cases, use '--ncpus' for parallel computation.

[usage]
$ ViPTreeGen [options] <input fasta> <output dir>

- <input fasta> should be in nucleotide FASTA format and include at least 3 sequences.
  Gzip-compressed input is accepted (auto-detected; `.gz`, `.fasta.gz`, or plain).
- If sequence name (before the first space in the header line) includes a character other than
  alphabets, numbers, dot(.), hyphen(-), or underscore(_), it will be replaced with underscore.
- <output dir> should not exist.

[dependencies]
    - tblastx                 -- included in the BLAST+ program;
                                 https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
    - ruby (ver >=3.0)
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
                                                       prior run, otherwise --resume is rejected. For --mode mmseqs-tblastx,
                                                       resume is also chunk-level (per --mmseqs-tblastx-chunk-size). For
                                                       --mode tblastx/blastn it is batch-level (per per-query/split file).

  (with-reference mode)
    --ref-duckdb   [path/run.duckdb]                -- a previous ViPTreeGen run.duckdb that supplies the reference set.
                                                       ref sequences, self_scores, and ref-vs-ref summary_tsv are imported
                                                       from it; input is searched against (input ∪ ref) so that input×ref
                                                       pairs are computed but ref×ref pairs are NOT recomputed. Output is
                                                       the combined sim/dist matrix and tree over (ref ∪ input). Requires
                                                       ref schema_version >= 2.0. Mutually exclusive with --2D.

  (parallelization)
    --ncpus        [int]                            -- number of CPUs to use for the all-vs-all computation.
                                                       In --mode tblastx/blastn: number of per-query search jobs run
                                                       in parallel via GNU parallel. In --mode mmseqs-tblastx/last:
                                                       the search engine's internal thread count.

[output files]
  (normal mode)
    result/all.sim.matrix              -- similarity (SG score) matrix
    result/all.dist.matrix             -- distance (1-SG score) matrix
    result/all.[bio]nj.[a|de]sc.newick -- Newick trees, midpoint rooted and ladderized (a proteomic tree in
                                          --mode tblastx / mmseqs-tblastx; a nucleotide tree in --mode blastn / last)

  (2D mode)
    result/2D.sim.matrix               -- similarity (SG score) matrix (row-wise: query fasta, column-wise: input fasta)
    result/top10.sim.list              -- top10 SG scores for each query sequence. format: 1st column - query ID; from 2nd to 11st columns - ID:SG_score

  (aggregated state)
    run.duckdb                         -- single DuckDB file containing summary_pre, self_scores, summary_tsv, sequences,
                                          run_metadata, tools (path/version of tools used), logs. Inspect with `duckdb run.duckdb`.
                                          Per-node search output files are deleted after the Rust binary consumes them.
```

## what changed in v2.0.0

Four headline changes

- **`--mode {mmseqs-tblastx|tblastx}`** — selectable proteomic-tree engine; `mmseqs-tblastx` is the new default (~10× faster than tblastx).
- **Rust binary** — the legacy Ruby steps are fused into one parallel `viptreegen-summary-pre`.
- **`--mode {blastn|last}`** — nucleotide-tree modes for the **DiGAlign** backend (a nucleotide tree, NOT a proteomic tree) (see [`doc/engine-benchmark.md`](doc/engine-benchmark.md)).
- **`--ref-duckdb`** — with-reference mode: reuse a prior run's reference set so only input×input and input×ref are recomputed.

Full details — benchmarks, per-mode comparison tables, and the with-reference workflow — are in [`doc/v2.0.0_update.md`](doc/v2.0.0_update.md) (the complete itemized list is in [`CHANGELOG.md`](CHANGELOG.md)).

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
