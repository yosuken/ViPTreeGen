
# ViPTreeGen - a standalone tool for viral proteomic tree generation

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](/LICENSE)
[![size](https://img.shields.io/github/size/webcaetano/craft/build/phaser-craft.min.js.svg)]()
[![doi](https://img.shields.io/badge/doi-10.1093%2Fbioinformatics%2Fbtx157-blue.svg?style=flat)](https://doi.org/10.1093/bioinformatics/btx157)

## description
ViPTreeGen has been developed as [the ViPTree server project](http://www.genome.jp/viptree).

## requirements
* BLAST+
* Ruby (ver >=2.0)
* R (ver >=3.0)

## usage 
```
### ViPTreeGen ver 1.0.1 (2018-02-08) ###

[description]
ViPTreeGen - tool for viral proteomic tree generation from viral genomic sequences.
ViPTreeGen has been developed as the ViPTree server project (http://www.genome.jp/viptree).

If you compute many sequences (e.g. >100) or large sequences (e.g. NCLDV genomes), it may take a long time.
In those cases, use '--ncpus' or '--queue' for parallel computating.

[usage]
$ ViPTreeGen <input fasta> <output dir> [options]

[dependencies]
    - tblastx              -- included in the BLAST+ program;
                              https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
    - ruby (ver >=2.0)

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

  (use GNU parallel)
    --ncpus        [int]                            -- number of jobs in parallel

  (for icr user)                                    -- for computation in the ICR supercomputer system
    --queue        [JP1]                            -- queue for computation

[output files]
result/all.sim.matrix              -- similarity (SG score) matrix
result/all.dist.matrix             -- distance (1-SG score) matrix
result/all.[bio]nj.[a|de]sc.newick -- Newick files of the viral proteomic tree, midpoint rooted and ladderized
                                      asc: nodes with fewer children sort before nodes with more children.
                                      desc: nodes with more children sorting before nodes with fewer children.
```

## citation
If you use results (data / figures) genereted by ViPTree in your research, please cite:
```
ViPTree: the viral proteomic tree server. Bioinformatics 33:2379–2380 (2017), doi:10.1093/bioinformatics/btx157
Yosuke Nishimura, Takashi Yoshida, Megumi Kuronishi, Hideya Uehara, Hiroyuki Ogata, and Susumu Goto
```
[Author manuscript](http://www.genome.jp/viptree/img/AM_Nishimura_Bioinformatics_2017.pdf) is freely available. 
