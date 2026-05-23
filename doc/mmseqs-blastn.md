# Why `--mode mmseqs-blastn` was removed

**TL;DR**: On the kind of input sizes ViPTreeGen and DiGAlign actually run on (tens to a few hundred viral genomes, ~10-100 kb each), `mmseqs search --search-type 3` (nucleotide vs nucleotide) was **slower, more memory-hungry, and less sensitive** than NCBI `blastn`. We removed `--mode mmseqs-blastn` so the only nucleotide engine is `blastn`. This document records the evidence so we don't re-litigate it.

## Context

When `--mode blastn` was added (DiGAlign backend), `--mode mmseqs-blastn` was added alongside it for symmetry with the proteomic-tree pair `tblastx` / `mmseqs-tblastx`, on the assumption that "mmseqs is faster than BLAST+ for the equivalent search". This is true for translated 6-frame protein alignment (where `mmseqs-tblastx` gives a ~10× speedup over `tblastx`) but turned out NOT to be true for nucleotide search at our scale.

## Test setup

- Dataset: `testdata/1.fasta` concatenated with `testdata/2.fasta`, deduped on seq_id → 15 viral plasmid genomes, ~1.17 Mb total
- Threads: 1 (single-thread comparison)
- Same e-value (1e-2), same `--max-seqs 1000000`
- Modules: `blast+/2.17.0`, `mmseqs2/13-45111`
- Measurement: `/usr/bin/time -f "wall=%e RSS=%M"` on the search step only (DB build and `convertalis` excluded)

## Headline numbers

| | Wall | Peak RAM | Total HSPs |
|---|---:|---:|---:|
| `blastn` | **1.3 s** | **103 MB** | 4,067 |
| `mmseqs --search-type 3` (default) | 11.3 s | 8.4 GB | 5,880 |
| `mmseqs --search-type 3 --alignment-mode 4` | 9.0 s | 8.4 GB | 6,095 |

On this size of input, `blastn` is **~9× faster** and uses **~80× less RAM**. mmseqs has a high fixed setup cost (index build + memory allocation) that does not amortize at this scale; it would amortize if you were searching a single 1 Mb query against a pre-built 100 Gb DB many times in one process, but ViPTreeGen and DiGAlign use one-shot all-vs-all instead.

## Sensitivity (the more important finding)

Comparison on one representative pair (`CP028145.1` vs `CP033140.1`, 75 kb genomes, ~95% identity):

| Mode | # HSPs | Σ bit-score | Σ alen | longest HSP |
|---|---:|---:|---:|---:|
| `blastn` | 17 | **122,196** | **72,042 bp** | **22,015 bp** |
| `mmseqs --search-type 3` (default) | 22 | 97,974 | 57,415 bp | 8,184 bp |
| `mmseqs --search-type 3 --alignment-mode 4` | 23 | 62,826 | 35,511 bp | 6,150 bp |

After interval merging in the Rust binary, this becomes:

| Mode | `que_score` | `que_len_in_hit` | %query covered |
|---|---:|---:|---:|
| `blastn` | 114,609 | 66,630 bp | 90% |
| `mmseqs default` | 95,511 | 55,914 bp | 76% |

mmseqs misses ~16% of the query coverage that blastn finds, **on a clearly homologous pair at ~98% mean identity**. The HSP-length distribution shows why: blastn produces single ~22 kb alignments that stretch through small divergent regions, while mmseqs fragments the same region into several shorter (~6-8 kb) HSPs and stops extending earlier. This is `blastn`'s banded extension with X-drop, which we cannot reproduce with any mmseqs parameter.

## What we tried before giving up

### Reducing `-k` (k-mer size for the prefilter)

mmseqs default is `-k 15` for nucleotide. Smaller k = more sensitive prefilter but slower.

| `-k` | wall | Σ bit-score (the pair) | Σ alen |
|---:|---:|---:|---:|
| 15 (default) | 7.04 s | 97,974 | 57,415 bp |
| 13 | 1.13 s | 98,088 | 57,502 bp |
| 11 | 0.74 s | 98,133 | 57,549 bp |
| 9 | 0.77 s | 98,133 | 57,549 bp |
| 7 | 1.00 s | 98,133 | 57,549 bp |

Plateau at `-k 11`. Lowering `-k` further does not help. The 57 kb mmseqs ceiling stays ~14 kb below blastn's 72 kb. The prefilter is not the limiting factor.

### `-s 7.5 -c 0` (max sensitivity, no coverage filter)

Identical HSP count and Σ bit-score to default. The prefilter setting was not gating these pairs.

### `--alignment-mode 4` (full Smith-Waterman alignment)

Actually got *worse* (Σ alen 57 kb → 36 kb). The alignment-mode flag does not control extension length; if anything it makes mmseqs more conservative about which alignments to keep.

### Other knobs not exhaustively explored

- `--gap-open` / `--gap-extend` (gap penalties) — could in principle help by encouraging longer alignments through divergent regions, but tuning these without a clear scoring-model justification would just produce a non-standard nucleotide score that doesn't match anything else.
- `--num-iterations` — iterative search; would help homology recall but does not address the extension-length issue.

These knobs were not pursued because the headline speed result (mmseqs 9× slower than blastn at our scale) made the cost-benefit moot.

## Correlation across modes

End-to-end comparison of `result/all.sim.matrix` from all four modes on `testdata/1.fasta` (15 sequences):

| pair | Pearson r | Spearman r |
|---|---:|---:|
| `tblastx` vs `mmseqs-tblastx` | 0.994 | 0.958 |
| `tblastx` vs `blastn` | 0.982 | 0.931 |
| `tblastx` vs `mmseqs-blastn` | 0.664 | 0.527 |
| `mmseqs-tblastx` vs `blastn` | 0.988 | 0.982 |
| `mmseqs-tblastx` vs `mmseqs-blastn` | 0.683 | 0.577 |
| `blastn` vs `mmseqs-blastn` | **0.702** | **0.574** |

The pattern is striking: **same-engine comparisons across alphabets (tblastx vs blastn: r=0.98)** are far better than **same-alphabet comparisons across engines (blastn vs mmseqs-blastn: r=0.70)**. This says the BLAST+ tblastx and blastn produce highly consistent scores with each other, and so do the mmseqs translated and (other) tools, but mmseqs's nucleotide-vs-nucleotide search is genuinely producing a different similarity ranking, not just a re-scaled version of blastn's. Even Spearman is low, so trees built from `mmseqs-blastn` would differ from `blastn` trees at the topology level — that's a non-starter for DiGAlign.

## Why mmseqs-tblastx works but mmseqs-blastn does not

mmseqs's design wins on translated protein search because:

- Genomes contain many distinct protein-length k-mers (a 15-mer in protein space has 20^15 alphabet entropy → very few false positives).
- The protein alphabet is small (20), and BLOSUM-scored alignments have a well-tuned X-drop cutoff that mmseqs's translated extension mimics closely.

For nucleotide-vs-nucleotide on viral genomes:

- The alphabet is 4 → k-mers are much less discriminating; the prefilter does less work, and mmseqs ends up doing many alignments anyway.
- blastn has 25+ years of tuning specifically for nucleotide extension through divergent regions (megablast defaults, X-drop, X-drop-final). mmseqs does not replicate this and there is no parameter that would force it to.

## Decision

`--mode mmseqs-blastn` was removed. The four-mode matrix is now three:

| Mode | algorithm | output tree | notes |
|---|---|---|---|
| `mmseqs-tblastx` (default) | translated 6-frame, MMseqs2 | proteomic | recommended for ViPTree analyses |
| `tblastx` | translated 6-frame, BLAST+ | proteomic | byte-exact pre-v2.0 reproduction |
| `blastn` | nucleotide-vs-nucleotide, BLAST+ | nucleotide | DiGAlign backend |

If a future use case needs nucleotide all-vs-all on >10,000 sequences and the mmseqs index-build overhead is finally amortized, `mmseqs-blastn` can be reintroduced — but it would need to ship with non-default parameters that close the sensitivity gap, and we don't currently know what those are.

## Reproducer

```bash
# On this dev host (lustre/aptmp/yosuke/dev/tool/ViPTreeGen):
cat testdata/1.fasta testdata/2.fasta > /tmp/bench29.fasta
awk '/^>/{seen[$1]++; keep=(seen[$1]==1)} keep' /tmp/bench29.fasta > /tmp/bench_uniq.fasta

# blastn
source /usr/appli/freeware/Modules/init/bash
module load blast+/2.17.0 mmseqs2/13-45111

mkdir -p /tmp/bench_blastn
makeblastdb -in /tmp/bench_uniq.fasta -dbtype nucl -out /tmp/bench_blastn/db
/usr/bin/time -f "wall=%es RSS=%MKB" \
  blastn -query /tmp/bench_uniq.fasta -db /tmp/bench_blastn/db \
  -outfmt 6 -evalue 1e-2 -max_target_seqs 1000000 -num_threads 1 \
  -out /tmp/bench_blastn/out.tsv

# mmseqs
mkdir -p /tmp/bench_mm
mmseqs createdb /tmp/bench_uniq.fasta /tmp/bench_mm/db
/usr/bin/time -f "wall=%es RSS=%MKB" \
  mmseqs search /tmp/bench_mm/db /tmp/bench_mm/db /tmp/bench_mm/res /tmp/bench_mm/tmp \
  --search-type 3 -e 1e-2 --max-seqs 1000000 --threads 1
```
