#!/usr/bin/env Rscript --no-save --no-restore
#
#  03-2.matrix_to_nj.R / ViPTreeGen - a command line tool for viral proteomic tree generation
#
#    Copyright: 2017 (C) Yosuke Nishimura (yosuke@kuicr.kyoto-u.ac.jp)
#    License: MIT license
#    Initial version: 2017-02-07
#

## usage
# Rscript <this.R> <fin> <foutpref> <flag> <method>
# fin:      sim/dist matrix
# foutpref: prefix for output Newick tree
# flag:     "sim" or "dist"
# method:   "nj" or "bionj"

## library
library(ape)
library(phangorn)

## args
args <- commandArgs(trailingOnly = TRUE)

## input
fin      = args[1]
foutpref = args[2]
flag     = args[3]
method   = args[4]

## parse matrix
df = read.delim(fin, row.names=1)
m  = as.matrix(df)

## convert similarity to distance by substracting from 1
if (flag == "sim"){
	m = 1 - m
}

## calc bionj tree
if (method == "nj"){
	tr = nj(m)
} else if (method == "bionj"){
	tr = bionj(m)
}

## midpoint rooting
tr = phangorn::midpoint(tr)

## ladderize
tr1 = ape::ladderize(tr)
tr2 = ape::ladderize(tr, right=F)

## make output
fo1 = paste0(foutpref, ".desc.newick")
fo2 = paste0(foutpref, ".asc.newick")
write.tree(tr1, file = fo1)
write.tree(tr2, file = fo2)
