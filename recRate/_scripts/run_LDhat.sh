#!/bin/bash

##### Bash script to estimate recMap with LDhat
##### author: Arka Pal
##### written: 29.04.2024
##### update: 16.05.2024

##### USAGE: bash run_LDhat.sh <options>
##### $1 - chrom; Chr6
##### $2 - sampPrefix; MF3.pFR.hCov
##### $3 - sampFile; 
##### $4 - number of cores
##### $5 - chrom window

## Load modules
module load bcftools vcftools ldhat
LDhat=/mnt/nfs/clustersw/Debian/bookworm/ldhat/20240426/bin

## Read inputs
stitchRun=stitchRun1
chrom=$1
sampPrefix=$2
sampleFile=$3
# sampleFile=~/snap_hap/sample_info/samples_by_pool/all_${flank}.list
ncores=$4
# ncores=10
window=$5

## Define variables
baseDIR=~/snap_hap/recRate/${chrom}.${sampPrefix}
stitchVCF=~/snap_hap/variants/stitch/$chrom/Am_all_${stitchRun}_${chrom}.final.vcf.gz
buffer=5000

# chromStart
if [ $window -eq 1 ]; 
then 
    # echo "$window is equal to 1"
    start=$((($window-1) * 50000 + 1)); 
else 
    # echo "$window greater than to 1"
    start=$((($window-1) * 50000 + 1 - $buffer));
fi 

# chromEnd
end=$((($window) * 50000 + $buffer)) # chromEnd
num_indv=$(cat $sampleFile | wc -l)
theta=0.01

## Make base directory
if [ ! -d $baseDIR ]; then mkdir $baseDIR; mkdir $baseDIR/windows; fi
cd $baseDIR

## Preprocessing VCF file
bcftools view -r ${chrom}:${start}-${end} -S $sampleFile $stitchVCF -Oz -o $baseDIR/windows/w$window.$start.$end.vcf.gz
# Create .sites and .locs file
vcftools --gzvcf $baseDIR/windows/w$window.$start.$end.vcf.gz \
		 --chr $chrom \
		 --ldhat-geno \
		 --out $baseDIR/windows/w$window.$start.$end


## Create lookup file
#$LDhat/complete -n 25 -rhomax 100 -n_pts 101 -theta 0.01 -split $ncores -element 0 ##didn't run
## Needs to be run once!
# pre_lkFile=$LDhat/lk_n100_t0.01
# lkgen -lk $pre_lkFile -nseq $((num_indv*2))
# mv new_lk.txt lkTable.txt
# new_lkFile=lkTable.txt

## Interval mapping
seqFile=$baseDIR/windows/w$window.$start.$end.ldhat.sites
locFile=$baseDIR/windows/w$window.$start.$end.ldhat.locs
outPrefix=$baseDIR/windows/w$window.$start.$end.
niter=1000000
samp=5000
bpen=5
interval -seq $seqFile -loc $locFile -lk $new_lkFile -its $niter -samp $samp -bpen $bpen -prefix $outPrefix
stat -input $outPrefix.rates.txt -burn 5 -loc $locFile -prefix $outPrefix.interval

## Rhomap
# seqFile=$baseDIR/$chrom/$flank.$chrom.ldhat.sites
# locFile=$baseDIR/$chrom/$flank.$chrom.ldhat.locs
# niter=1000000
# samp=5000
# burn=200000
# rhomap -seq $seqFile -loc $locFile -lk $new_lkFile -its $niter -sampe $samp -burn $burn -prefix $prefix