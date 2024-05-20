#!/bin/bash

##### Bash script to extract positions of SNPs to create posfile
##### author: Arka Pal
##### date: 10.08.2023
##### update: 10.08.2023

##### USAGE: bash extract_posfile.sh <options>
##### $1 - chromosome no. 
##### $2 - chromosome start position (bp)
##### $3 - chromosome end position (bp)
##### $4 - buffer (bp)
##### $5 - inVCF
##### $6 - output posfile

##### posfile: text file with columns- <chr> <pos> <REF> <ALT>


## Load modules
module load bcftools/1.18

## Read input
chr=$1
chromStart=$2
chromEnd=$3
buffer=$4
inVCF=$5
outfile=$6

## Initiate variables
if [ $chromStart -lt $buffer ];
then
    chromStart_withBuffer=1;
else
    chromStart_withBuffer=$((chromStart-buffer))
fi
chromEnd_withBuffer=$((chromEnd+buffer))
region=${chr}:${chromStart_withBuffer}-${chromEnd_withBuffer}

## Create posfile
bcftools view -H -r $region $inVCF | cut -f1,2,4,5 > $outfile