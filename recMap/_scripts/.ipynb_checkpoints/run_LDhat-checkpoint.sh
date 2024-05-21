#!/bin/bash

##### Bash script to estimate recMap with LDhat
##### author: Arka Pal
##### written: 29.04.2024
##### update: 16.05.2024

##### USAGE: bash run_LDhat.sh <options>
##### $1 - chrom
##### $2 - flank
##### $3 - num_of_cores

## Load modules
module load bcftools vcftools ldhat

## Read inputs
stitchRun=stitchRun1
flank=$1
chrom=$2
ncores=$3

## Define variables
baseDIR=~/snap_hap/recMap/$flank/$chrom
LDhat=/mnt/nfs/clustersw/Debian/bookworm/ldhat/20240426/bin
stitchVCF=~/snap_hap/variants/stitch/$chrom/Am_all_${stitchRun}_${chrom}.final.vcf.gz
sampleFile=~/snap_hap/sample_info/samples_by_pool/all_${flank}.list
# ncores=10
num_indv=$(cat $sampleFile | wc -l)
theta=0.01

## Preprocessing VCF file
cd $baseDIR
# bcftools view -r $chrom $stitchVCF -S $sampleFile -Oz -o $baseDIR/$flank.$chrom.vcf.gz
# vcftools --gzvcf $baseDIR/$flank.$chrom.vcf.gz \
# 		 --chr $chrom \
# 		 --ldhat-geno \
# 		 --out $baseDIR/$flank.$chrom
# This will output .sites and .locs file

## Create lookup file
#$LDhat/complete -n 25 -rhomax 100 -n_pts 101 -theta 0.01 -split $ncores -element 0 ##didn't run
# pre_lkFile=$LDhat/lk_n100_t0.01
# lkgen -lk $pre_lkFile -nseq $((num_indv*2))
# new_lkFile=new_lk.txt

## Interval mapping
seqFile=$baseDIR/$flank.$chrom.ldhat.sites
locFile=$baseDIR/$flank.$chrom.ldhat.locs
outPrefix=$baseDIR/$flank.$chrom
niter=2000
samp=100
bpen=5
# interval -seq $seqFile -loc $locFile -lk $new_lkFile -its $niter -samp $samp -bpen $bpen -prefix $outPrefix
stat -input $outPrefix.rates.txt -burn 5 -loc $locFile -prefix $outPrefix.interval

## Rhomap
# seqFile=$baseDIR/$chrom/$flank.$chrom.ldhat.sites
# locFile=$baseDIR/$chrom/$flank.$chrom.ldhat.locs
# niter=1000000
# samp=5000
# burn=200000
# rhomap -seq $seqFile -loc $locFile -lk $new_lkFile -its $niter -sampe $samp -burn $burn -prefix $prefix