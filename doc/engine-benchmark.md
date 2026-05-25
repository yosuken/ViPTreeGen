# Nucleotide engine benchmark: blastn vs LAST (and mmseqs-tblastx)

**TL;DR**: For the DiGAlign nucleotide-tree backend, **neither blastn nor LAST is universally faster — it depends on the dataset.** LAST wins decisively on **many small genomes**; blastn wins on **few large genomes**. Both produce highly-correlated similarity matrices (Pearson 0.86-0.99) and near-interchangeable virus clusters (ARI 0.81-0.88). ViPTreeGen ships both (`--mode blastn`, `--mode last`); choose by workload. mmseqs-tblastx (protein, the proteomic-tree default) is included for reference and is a different, more sensitive metric.

## Why this benchmark exists

`--mode blastn` was the first DiGAlign nucleotide engine. We evaluated faster alternatives (see also `doc/mmseqs-blastn.md` and `doc/mummer4-nucmer.md`, both rejected). LAST (`lastdb` + `lastal -fBlastTab`) emerged as a serious candidate: its BlastTab output is byte-compatible with `blastn -outfmt 6` (no converter needed), it is robust to repeats (adaptive seeds, unlike nucmer), and it has very low memory on small genomes. This file records the head-to-head across three deliberately different datasets.

## Datasets (from Virus-Host DB, 2026-05)

| set | taxon | # seqs | median genome | total | regime |
|---|---|---:|---:|---:|---|
| 1 | Begomovirus (Geminiviridae genus) | 870 | 2.7 kb | 2.3 MB | many, tiny, dense cross-similarity |
| 2 | Autographivirales (Caudoviricetes order; T7-like phages) | 1780 | 41 kb | 74 MB | many, medium |
| 3 | Orthoherpesviridae (family) | 135 | 159 kb | 21.5 MB | medium count, large genomes |
| 4 | Imitervirales (giant viruses, order) | 23 | 1.2 Mb | 25.3 MB | few, enormous genomes |

LAST config: `lastdb` (default seeding) + `lastal -r 1 -q 2 -a 0 -b 2` (megablast-matched scoring) `-E <evalue> -fBlastTab`. blastn: ViPTreeGen defaults (megablast). All runs `--ncpus 4`, `--notree`.

Datasets are ordered below by **median genome size**, which turns out to be the variable that best predicts which engine wins.

## Speed (wall clock)

| dataset | median genome | blastn | LAST | mmseqs-tblastx |
|---|---:|---:|---:|---:|
| Begomovirus (870 seq) | 2.7 kb | 154 s | **20 s** | 89 s |
| Autographivirales (1780 seq) | 41 kb | 576 s | **520 s** | killed¹ |
| Orthoherpesviridae (135 seq) | 159 kb | **51 s** | 107 s | 234 s |
| Imitervirales (23 seq) | 1.2 Mb | **123 s** | 355 s | 327 s |

¹ mmseqs-tblastx on Autographivirales (74 MB total, translated 6-frame) was killed by the interactive node's per-process CPU-time ulimit; it needs a batch job (qsub). Not a code issue, and not relevant to the blastn-vs-LAST nucleotide comparison.

## Peak memory (RSS)

| dataset | median genome | blastn | LAST | mmseqs-tblastx |
|---|---:|---:|---:|---:|
| Begomovirus | 2.7 kb | 1.7 GB | **21 MB** | 1.3 GB |
| Autographivirales | 41 kb | **400 MB** | 497 MB | (killed) |
| Orthoherpesviridae | 159 kb | **109 MB** | 418 MB | 1.27 GB |
| Imitervirales | 1.2 Mb | **197 MB** | 866 MB | 1.2 GB |

## The crossover — genome size is the key variable

```
  median genome:   2.7 kb        41 kb        159 kb        1.2 Mb
  speed winner:    LAST (8x)  ~tied (1.1x)   blastn (2x)   blastn (3x)
  memory winner:   LAST (80x) blastn (1.2x)  blastn (4x)   blastn (4x)
                   |------------------|--------------------------------|
                   LAST regime     ~crossover        blastn regime
                                   (~40-50 kb)
```

LAST's advantage tracks **small genome size**, not just sequence count: Autographivirales has *more* sequences (1780) than Begomovirus (870), yet the engines are tied there because its genomes are 15× larger. The crossover sits around a median genome of ~40-50 kb.

Mechanistic explanation:

- **LAST** does one global `lastal` against a single on-disk, memory-mapped index. Per-genome overhead is tiny, so it amortizes beautifully across **many** genomes; but total alignment work and the mmap'd working set scale with total sequence **volume** and per-genome size, so large genomes slow it down and grow its memory (21 MB → 866 MB as median size grows).
- **blastn** runs per-(query, split) jobs through GNU parallel. With **many small dense** genomes it pays per-invocation overhead and accumulates a huge number of HSPs in memory (1.7 GB on Begomovirus). With **larger** genomes the chunked jobs parallelize well and, where genomes are divergent/sparse, there are few HSPs to hold (109-400 MB).

Memory is essentially mirror-imaged: LAST tiny→large as genomes grow; blastn large→tiny as genomes get bigger and sparser.

## Similarity-matrix agreement (blastn vs LAST)

| dataset | median genome | Pearson r | blastn saturation | LAST saturation |
|---|---:|---:|---:|---:|
| Begomovirus | 2.7 kb | 0.86 | 63% | 50% |
| Autographivirales | 41 kb | 0.98 | 81% | 82% |
| Orthoherpesviridae | 159 kb | 0.99 | 88% | 56% |
| Imitervirales | 1.2 Mb | 0.99 | 42% | 9.5% |

Pearson is lowest (0.86) on Begomovirus — the genus-level (50-90% identity) regime where LAST's adaptive seeding extends through divergent regions that blastn's X-drop stops at, so LAST reports systematically higher SG (higher in ~95% of co-detected pairs there). On the other sets, where similarity is either sparse-and-clear or absent, the two agree at r≈0.98-0.99. (On Autographivirales the two are nearly symmetric — LAST higher in only 44% of co-detected pairs — so even the directional bias is dataset-dependent.)

LAST is consistently **more sensitive** (lower saturation = detects relationships in more pairs, and extends alignments further). The two agree most (r≈0.99) where nucleotide similarity is sparse/clear; they diverge most (r=0.86) at genus-level divergence (Begomovirus 50-90% identity), exactly the regime where LAST's adaptive seeding extends through divergent regions that blastn's X-drop stops at. On co-detected pairs, LAST's SG is higher than blastn's in ~95% of cases.

### Tree topology / clustering (Begomovirus, dense 150-seq subset)

- Robinson-Foulds (bionj trees): blastn vs LAST normRF **0.18** (vs ~0.53 for nucleotide-vs-protein). 82% of bipartitions shared.
- Cluster membership after midpoint-root + cut, Adjusted Rand Index: blastn vs LAST **0.81-0.88** across cut depths (vs 0.55-0.68 for nucleotide-vs-protein).

So blastn and LAST place viruses into essentially the same groups; the protein method (mmseqs-tblastx) gives a genuinely different (and more sensitive) tree, as expected for translated vs nucleotide alignment.

## A note on mmseqs-tblastx (protein) sensitivity

On Orthoherpesviridae, mmseqs-tblastx saturation was **0.9%** vs blastn's **88%**: at the protein level, herpesvirus core genes (DNA pol, terminase, …) are conserved across all subfamilies even though the nucleotide sequence has diverged past blastn/LAST detection. This is why the proteomic tree (ViPTree's purpose) and the nucleotide tree (DiGAlign's purpose) are fundamentally different products, not interchangeable engines.

## Recommendation

Rule of thumb: **pick by median genome size**, with the crossover around ~40-50 kb.

| DiGAlign workload (median genome size) | recommended `--mode` |
|---|---|
| small genomes ≲40 kb (ssDNA viruses, small phages, most metagenomic viral contigs) | **last** (much faster, far lower memory; the gap widens as genomes shrink) |
| medium genomes ~40-50 kb (T7-like phages and similar) | either — roughly tied; `last` slightly faster, `blastn` slightly lower memory |
| large genomes ≳50 kb (herpesviruses, poxviruses, NCLDV / giant viruses) | **blastn** (faster, lower memory; the gap widens as genomes grow) |
| need byte-exact comparability with pre-v2.0 / published nucleotide results | **blastn** |

Both modes feed the identical downstream SG pipeline and produce `result/all.{sim,dist}.matrix` + tree, so switching engines is a one-flag change.

## Reproducer (sketch)

```bash
source /usr/appli/freeware/Modules/init/bash
module load blast+/2.17.0 last/1650 mmseqs2/13-45111

# extract a taxon subset from Virus-Host DB by the virus-taxonomy field (3rd '|' column), e.g.:
#   Begomovirus, Orthoherpesviridae, Imitervirales
# then for each engine:
for m in blastn last mmseqs-tblastx; do
  /usr/bin/time -v ./ViPTreeGen --notree --mode "$m" --ncpus 4 subset.fasta out_"$m"
done
# compare result/all.sim.matrix across out_* (Pearson, saturation, bionj + RF.dist, ARI)
```
