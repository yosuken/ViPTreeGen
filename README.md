
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
* BLAST+
* Ruby (ver >=2.0)
* R (ver >=3.0)

## usage 
```
### ViPTreeGen ver 1.1.0 (2018-02-15) ###

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
                                          asc: (nodes with fewer children) --> (nodes with more children).
                                          desc: (nodes with more children) --> (nodes with fewer children).

  (2D mode)
    result/2D.sim.matrix               -- similarity (SG score) matrix (row-wise: query fasta, column-wise: input fasta)
    result/top10.sim.list              -- top10 SG scores (and relevant sequences) for each query sequence
```

## citation
If you use results (data / figures) genereted by ViPTree in your research, please cite:
```
ViPTree: the viral proteomic tree server. Bioinformatics 33:2379–2380 (2017), doi:10.1093/bioinformatics/btx157
Yosuke Nishimura, Takashi Yoshida, Megumi Kuronishi, Hideya Uehara, Hiroyuki Ogata, and Susumu Goto
```
[Author manuscript](http://www.genome.jp/viptree/img/AM_Nishimura_Bioinformatics_2017.pdf) is freely available. 
