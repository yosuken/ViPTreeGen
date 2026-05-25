# Why MUMmer4 (`nucmer`) was not adopted as an alternative `blastn` engine

**TL;DR**: `nucmer` is data-dependent in ways that make it unreliable as a drop-in faster `blastn`. Its default mode rejects k-mer anchors that repeat in the reference, so repetitive regions (and any kind of intra-fasta sequence duplication) silently drop out. The fallback `--maxmatch` mode handles repeats but is dramatically slower on those same inputs. We do not ship a `--mode nucmer` until we have a real-workload benchmark (e.g. ICTV / EVG) that demonstrates a clear win.

## Context

`--mode blastn` (NCBI BLAST+) is the current nucleotide-vs-nucleotide engine for the DiGAlign backend. At small scale it is fast (<1 min for tens of viral genomes) but at large scale (10k+ sequences) it would be desirable to have a faster alternative. MUMmer4's `nucmer` is well-known for whole-genome alignment of bacterial chromosomes and is suffix-array based, so a priori it should outperform BLAST+'s heuristic seeded extension on long sequences with sparse, conserved blocks.

## Test setup

- Module: `MUMmer/4.0.0rc1`, `blast+/2.17.0`
- Three datasets:
  - `ssDNA.prok.8.fasta`: 8 small ssDNA virus genomes, 54 KB total, mostly unrelated.
  - `1.fasta + 2.fasta` (deduped): 15 Vibrio plasmids, 1.17 MB total, with substantial cross-pair similarity.
  - Synthetic: 10× duplicated Vibrio (150 entries, 12 MB total) — intended as a "scaled-up" stress test.
- All commands ran with `--threads 4` where applicable. Timed with `/usr/bin/time -f "wall=%e CPU=%P"`.

## Benchmark numbers

| Dataset | blastn -t4 | nucmer default (`--mum-reference`) | nucmer `--mum` (unique both sides) | nucmer `--maxmatch` |
|---|---:|---:|---:|---:|
| ssDNA 8 seq, 54 KB | 0.56 s, 8 HSPs | **0.02 s**, 8 alignments | n/a | n/a |
| Vibrio 15 seq, 1.17 MB | 0.96 s, 4067 HSPs | 1.79 s, 148 alignments | 1.81 s | 4.23 s |
| Synth 150 seq, 12 MB | 28 s, 406 k HSPs | 2.9 s, **0 alignments** (silent failure) | (same as default for 10× duplicates: 0) | **378 s** at 394% CPU |

On the small unrelated set nucmer wins decisively (28× faster). On a small set with dense cross-pair similarity, blastn is faster. On the synthetic high-repeat set, both nucmer modes fail loudly or silently.

## Root cause of the silent zero on synthetic data

`nucmer --help`:

> By default, nucmer uses anchor matches that are **unique in the reference** but not necessarily unique in the query.

The synthetic dataset has 10 identical copies of every Vibrio sequence (renamed `SEQ_0001..SEQ_0150`). Every k-mer therefore appears at least 10 times in the reference, so **no anchor is unique** and nucmer's default seeding rejects everything. The result is a `.delta` file with 2 lines (header only) and 0 alignments. There is no warning and no nonzero exit -- the run "succeeds" but emits no data, which would corrupt downstream summary_pre population if we trusted it.

The same logic applies to `--mum` (which is stricter: unique in BOTH reference and query).

`--maxmatch` allows non-unique anchors and does find alignments, but every k-mer match generates O(reps²) candidate anchor pairs. With 10× duplicates that is 100 candidate pairs per matching k-mer, leading to combinatorial blowup in the cluster step. Wall time goes from 28 s (blastn baseline) to 378 s (13× slower).

## Why this matters for ViPTree-style workloads

Real ViPTree inputs are not "10× duplicates of the same plasmid". But they DO contain:

- **Repeat elements** within individual viral genomes (terminal repeats, prophage repeats, integrase regions).
- **Highly conserved gene-family motifs** that appear across many genomes in the input (e.g., capsid-protein-coding regions shared by tens of phages in a family).
- For ICTV-style **representative collections**, a fairly clean set with limited repetition.
- For metagenome-derived collections, much more redundancy and contamination.

`nucmer --mum-reference` (default) will systematically drop signal from any region that happens to be conserved across multiple input genomes -- precisely the regions ViPTree most wants to count toward genome similarity. That is the opposite of the desired behavior for a tree-construction backend.

We could force `--maxmatch` to be safe, but then the 13× slowdown observed on duplicates is a worst-case warning of what real repeat-heavy datasets could cost. Without an actual EVG-scale benchmark we cannot calibrate this.

## Implementation cost we would also pay

Even if the speed were favorable, integration is non-trivial:

- `show-coords -T -H -l -c` output: columns are `S1 E1 S2 E2 LEN1 LEN2 %IDY LEN_R LEN_Q COV_R COV_Q REF QUERY`. Order differs from BLAST tab (REF=subject is last, QUERY is after REF), so a converter script is required.
- `nucmer` does not emit bit score or e-value. The Rust binary in 02-3 (`viptreegen-summary-pre`) consumes columns including `bitscore` and `evalue`. We would have to synthesize a bit score (e.g. `alen * idy / 100 * k` for some constant) and zero out the e-value. The SG scores would then differ from blastn in absolute value -- not a deal-breaker (SG is a normalized ratio) but it means we lose bit-exact reproducibility against blastn output and would need a separate golden matrix.
- The `--mum` / `--mum-reference` / `--maxmatch` choice becomes a user-facing knob with significant biological consequences, which is exactly the kind of footgun ViPTreeGen has been trying to remove (see `doc/mmseqs-blastn.md`).

## Decision

`--mode nucmer` is **not implemented**. The three current modes (`mmseqs-tblastx`, `tblastx`, `blastn`) cover the proteomic-tree and nucleotide-tree (DiGAlign) use cases.

If a future use case justifies revisiting MUMmer4:

1. Real-workload benchmark first: ICTV (~6000 reps) or EVG (1811) vs blastn, looking at both wall clock and SG matrix Pearson r.
2. If `nucmer --maxmatch` is fast enough on real data, design around that mode and skip the `--mum-reference` default entirely.
3. Document the bit-score synthesis explicitly and ship a separate golden matrix.

For now, `--mode blastn` (NCBI BLAST+, the established tool for nucleotide HSP-based all-vs-all) is the recommended nucleotide engine.

## Reproducer

```bash
source /usr/appli/freeware/Modules/init/bash
module load MUMmer/4.0.0rc1 blast+/2.17.0

# small, sparse — nucmer wins
fa=testdata/ssDNA.prok.8.fasta
/usr/bin/time -f "blastn wall=%es" \
  blastn -query $fa -db <(makeblastdb -in $fa -dbtype nucl -out /tmp/db; echo /tmp/db) \
         -outfmt 6 -num_threads 4 -out /tmp/b.tsv >/dev/null 2>&1 || \
  ( makeblastdb -in $fa -dbtype nucl -out /tmp/db && \
    /usr/bin/time -f "blastn wall=%es" blastn -query $fa -db /tmp/db -outfmt 6 -num_threads 4 -out /tmp/b.tsv )
mkdir -p /tmp/nuc
/usr/bin/time -f "nucmer wall=%es" \
  nucmer --threads 4 -p /tmp/nuc/o $fa $fa

# dataset with intra-fasta duplicates -- default nucmer returns 0 alignments
python3 -c "
import sys
i = 0
for copy in range(10):
    with open('testdata/1.fasta') as f:
        for line in f:
            if line.startswith('>'):
                i += 1; print(f'>SEQ_{i:04d}')
            else: print(line, end='')
" > /tmp/dup.fasta
mkdir -p /tmp/nuc_dup
nucmer --threads 4 -p /tmp/nuc_dup/o /tmp/dup.fasta /tmp/dup.fasta
show-coords -T -H -l -c /tmp/nuc_dup/o.delta | wc -l   # expect 0
```
