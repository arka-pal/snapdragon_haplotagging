#### compareAlleles: Compute chromosome-wide allele comparisons between Amolle and Amajus
#### author: Arka Pal
#### written: 
#### last update: 

#### usage: Rscript compareAlleles_bash.R baseDIR chrom stitchRun outFile

## Read arguments from the bash
args <- commandArgs(trailingOnly = TRUE)
baseDIR <- as.character(args[1])
chrom <- as.character(args[2])
stitchRun <- as.character(args[3])
outFile <- as.character(args[4])

##
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(tictoc)

source("~/snap_hap/Amajus_alleles-ancestral/_scripts/functions_polarisation.R")

alleleSequence =  c('A','T','C','G','N','del')

alleleDat_chrom = compile_AlleleDat_chromSegments(alleleSequence, baseDIR, chrom, stitchRun)
fwrite(alleleDat_chrom, 
       file = outFile,
       sep = ',', 
       quote = F, 
       row.names = FALSE, 
       col.names = TRUE)