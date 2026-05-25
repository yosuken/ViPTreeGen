# Benchmark subsets

Gzipped FASTA subsets used in the nucleotide-engine benchmark
(`doc/engine-benchmark.md`). Each covers a different genome-size regime so the
`blastn` / `last` / `mmseqs-tblastx` engines can be exercised across the range
where their speed/memory trade-offs cross over.

| file | taxon | # seqs | median genome | regime |
|---|---|---:|---:|---|
| `Begomovirus.40.fasta.gz` | Begomovirus (Geminiviridae genus) | 40 | 2.7 kb | many tiny, dense cross-similarity |
| `Autographivirales.40.fasta.gz` | Autographivirales (T7-like phages) | 40 | 41 kb | medium |
| `Orthoherpesviridae.30.fasta.gz` | Orthoherpesviridae (family) | 30 | 159 kb | large genomes |
| `Imitervirales.10.fasta.gz` | Imitervirales (giant viruses) | 10 | 1.2 Mb | few enormous genomes |

## Provenance

Source: **Virus-Host DB** `virushostdb.genomic.fna` (2026-05 release).
Sequences were selected by the virus-taxonomy field (3rd `|`-delimited column of
the FASTA header) containing the taxon name, then **down-sampled
deterministically**: sort the matching accessions, take evenly-spaced indices
(`round(i*(n-1)/(k-1))`). This keeps the subset reproducible and spread across
the taxon rather than clustered. The full taxa hold 897 / 1780 / 135 / 23
sequences respectively; see `doc/engine-benchmark.md` for the full-set numbers.

The files are committed **gzipped** to keep the repo light (~5 MB total vs ~17 MB
uncompressed); ViPTreeGen reads gzip input transparently (detected by magic
bytes), so they can be passed directly:

```bash
./viptreegen --mode last testdata/benchmark/Begomovirus.40.fasta.gz out_bego
```

FASTA headers are kept intact for provenance; ViPTreeGen uses the first
whitespace-delimited token (the accession) as the sequence label.
