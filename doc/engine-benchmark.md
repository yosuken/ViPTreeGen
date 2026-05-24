# Nucleotide engine benchmark: blastn vs LAST (and mmseqs-tblastx)

**TL;DR**: For the DiGAlign nucleotide-tree backend, **neither blastn nor LAST is universally faster — it depends on the dataset.** LAST wins decisively on **many small genomes**; blastn wins on **few large genomes**. Both produce highly-correlated similarity matrices (Pearson 0.86-0.99) and near-interchangeable virus clusters (ARI 0.81-0.88). ViPTreeGen ships both (`--mode blastn`, `--mode last`); choose by workload. mmseqs-tblastx (protein, the proteomic-tree default) is included for reference and is a different, more sensitive metric.

## Why this benchmark exists

`--mode blastn` was the first DiGAlign nucleotide engine. We evaluated faster alternatives (see also `doc/mmseqs-blastn.md` and `doc/mummer4-nucmer.md`, both rejected). LAST (`lastdb` + `lastal -fBlastTab`) emerged as a serious candidate: its BlastTab output is byte-compatible with `blastn -outfmt 6` (no converter needed), it is robust to repeats (adaptive seeds, unlike nucmer), and it has very low memory on small genomes. This file records the head-to-head across three deliberately different datasets.

## Datasets (from Virus-Host DB, 2026-05)

| set | taxon | # seqs | genome size | total | regime |
|---|---|---:|---:|---:|---|
| 1 | Begomovirus (Geminiviridae genus) | 870 | ~2.7 kb | 2.3 MB | many, tiny, dense cross-similarity |
| 2 | Orthoherpesviridae (family) | 135 | ~159 kb (med) | 21.5 MB | medium count, large genomes |
| 3 | Imitervirales (giant viruses, order) | 23 | ~1.2 Mb (med) | 25.3 MB | few, enormous genomes |

LAST config: `lastdb` (default seeding) + `lastal -r 1 -q 2 -a 0 -b 2` (megablast-matched scoring) `-E <evalue> -fBlastTab`. blastn: ViPTreeGen defaults (megablast). All runs `--ncpus 4`, `--notree`.

## Speed (wall clock)

| dataset | blastn | LAST | mmseqs-tblastx |
|---|---:|---:|---:|
| Begomovirus (870 × 2.7 kb) | 154 s | **20 s** | 89 s |
| Orthoherpesviridae (135 × 159 kb) | **51 s** | 107 s | 234 s |
| Imitervirales (23 × 1.2 Mb) | **123 s** | 355 s | 327 s |

## Peak memory (RSS)

| dataset | blastn | LAST | mmseqs-tblastx |
|---|---:|---:|---:|
| Begomovirus | 1.7 GB | **21 MB** | 1.3 GB |
| Orthoherpesviridae | **109 MB** | 418 MB | 1.27 GB |
| Imitervirales | **197 MB** | 866 MB | 1.2 GB |

## The crossover

```
        LAST faster / lower-mem  <------------------->  blastn faster / lower-mem
        many small genomes                              few large genomes
        (Begomovirus: LAST 8x faster, 80x less RAM)     (Giant: blastn 3x faster, 4x less RAM)
```

Mechanistic explanation:

- **LAST** does one global `lastal` against a single on-disk, memory-mapped index. The per-genome overhead is tiny, so it amortizes beautifully across **many** genomes; but the total alignment work scales with total sequence volume, so **few but enormous** genomes make it slow, and its mmap'd working set grows (21 MB → 866 MB).
- **blastn** runs per-(query, split) jobs through GNU parallel. With **many small** genomes it pays per-invocation overhead and, on dense cross-similarity, accumulates a huge number of HSPs in memory (1.7 GB on Begomovirus). With **few large** genomes the chunked jobs parallelize well and, where genomes are divergent/sparse, there are few HSPs to hold (109-197 MB).

Memory is essentially mirror-imaged: LAST tiny→large as genomes grow; blastn large→tiny as genome count drops and sparsity rises.

## Similarity-matrix agreement (blastn vs LAST)

| dataset | Pearson r | blastn saturation (dist=1.0) | LAST saturation |
|---|---:|---:|---:|
| Begomovirus | 0.86 | 63% | 50% |
| Orthoherpesviridae | 0.99 | 88% | 56% |
| Imitervirales | 0.99 | 42% | 9.5% |

LAST is consistently **more sensitive** (lower saturation = detects relationships in more pairs, and extends alignments further). The two agree most (r≈0.99) where nucleotide similarity is sparse/clear; they diverge most (r=0.86) at genus-level divergence (Begomovirus 50-90% identity), exactly the regime where LAST's adaptive seeding extends through divergent regions that blastn's X-drop stops at. On co-detected pairs, LAST's SG is higher than blastn's in ~95% of cases.

### Tree topology / clustering (Begomovirus, dense 150-seq subset)

- Robinson-Foulds (bionj trees): blastn vs LAST normRF **0.18** (vs ~0.53 for nucleotide-vs-protein). 82% of bipartitions shared.
- Cluster membership after midpoint-root + cut, Adjusted Rand Index: blastn vs LAST **0.81-0.88** across cut depths (vs 0.55-0.68 for nucleotide-vs-protein).

So blastn and LAST place viruses into essentially the same groups; the protein method (mmseqs-tblastx) gives a genuinely different (and more sensitive) tree, as expected for translated vs nucleotide alignment.

## A note on mmseqs-tblastx (protein) sensitivity

On Orthoherpesviridae, mmseqs-tblastx saturation was **0.9%** vs blastn's **88%**: at the protein level, herpesvirus core genes (DNA pol, terminase, …) are conserved across all subfamilies even though the nucleotide sequence has diverged past blastn/LAST detection. This is why the proteomic tree (ViPTree's purpose) and the nucleotide tree (DiGAlign's purpose) are fundamentally different products, not interchangeable engines.

## Recommendation

| DiGAlign workload | recommended `--mode` |
|---|---|
| many small viral genomes (phages, ssDNA viruses, metagenomic viral contigs; hundreds–thousands) | **last** (much faster, far lower memory) |
| few large viral genomes (herpesviruses, poxviruses, NCLDV / giant viruses) | **blastn** (faster, lower memory) |
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
