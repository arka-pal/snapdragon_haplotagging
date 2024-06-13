#!/bin/bash

##### Bash script to estimate recMap with LDhat
##### author: Arka Pal
##### written: 29.04.2024
##### update: 07.06.2024

##### USAGE: bash run_LDhat.sh <options>
##### $1 - chrom; Chr6
##### $2 - sampPrefix; MF3.pFR.hCov
##### $3 - sampleFile; 
##### $4 - bpen
##### $5 - number of cores
##### $6 - chrom window


## Load modules
echo -e '\n\n ----- \nLoad modules\n -----\n'
module load bcftools vcftools ldhat
LDhat=/mnt/nfs/clustersw/Debian/bookworm/ldhat/20240426/bin


## Read inputs
echo -e '\n\n ----- \nRead inputs\n -----\n'
stitchRun=stitchRun1
chrom=$1
sampPrefix=$2
sampleFile=$3
# sampleFile=~/snap_hap/sample_info/samples_by_pool/all_${flank}.list
bpen=$4
ncores=$5
# ncores=10
window=$6


## Define variables
echo -e '\n\n ----- \nDefine variables\n -----\n'
baseDIR=~/snap_hap/recRate/${chrom}.${sampPrefix}
stitchVCF=~/snap_hap/variants/stitch/$chrom/Am_all_${stitchRun}_${chrom}.final.vcf.gz
windowSize=100000
buffer=5000
num_indv=$(cat $sampleFile | wc -l)
theta=0.01

windowStart=$((($window-1) * $windowSize + 1));
if [ $window -eq 1 ]; # start incl. buffer
then bufferStart=$((($window-1) * $windowSize + 1));
else bufferStart=$((($window-1) * $windowSize + 1 - $buffer));
fi

windowEnd=$((($window) * $windowSize));
bufferEnd=$((($window) * $windowSize + $buffer)) # end incl. buffer


## Make base directory
if [ ! -d $baseDIR ]; then mkdir $baseDIR; mkdir $baseDIR/windows; fi
cd $baseDIR


## Preprocessing VCF file
echo -e '\n\n ----- \nMake VCF & input files\n -----\n'
bcftools view -r ${chrom}:${bufferStart}-${bufferEnd} -S $sampleFile $stitchVCF -Oz -o $baseDIR/windows/w$window.$windowStart.$windowEnd.buf$buffer.vcf.gz
# Create .sites and .locs file
vcftools --gzvcf $baseDIR/windows/w$window.$windowStart.$windowEnd.buf$buffer.vcf.gz \
		 --chr $chrom \
		 --ldhat-geno \
		 --out $baseDIR/windows/w$window.$windowStart.$windowEnd.buf$buffer


echo -e '\n\n ----- \nInput sites\n -----\n'
echo -e 'chrom:' $chrom
echo -e 'windowStart:' $windowStart
echo -e 'windowEnd:' $windowEnd
echo -e 'bufferStart:' $bufferStart
echo -e 'bufferEnd:' $bufferEnd
echo -e 'No. of SNPs without buffer:' $(bcftools view -H $baseDIR/windows/w$window.$windowStart.$windowEnd.buf$buffer.vcf.gz | wc -l)
echo -e 'No. of SNPs incl. buffer:' $(bcftools view -H $baseDIR/windows/w$window.$windowStart.$windowEnd.buf$buffer.vcf.gz | wc -l)
echo -e 'No. of lines in siteFiles:' $(wc -l $baseDIR/windows/w$window.$windowStart.$windowEnd.buf$buffer.ldhat.sites)
echo -e 'No. of locs:' $(wc -l $baseDIR/windows/w$window.$windowStart.$windowEnd.buf$buffer.ldhat.locs)
echo -e 'No. of samples:' $num_indv


## Create lookup file
echo -e '\n\n ----- \nCreate lookup file\n -----\n'
#$LDhat/complete -n 25 -rhomax 100 -n_pts 101 -theta 0.01 -split $ncores -element 0 ##didn't run
## Needs to be run once!
lookup_Table=$baseDIR/lkTable.txt
if [ ! -e $lookup_Table ]; 
then
    pre_lkFile=$LDhat/lk_n100_t0.01
    lkgen -lk $pre_lkFile -nseq $((num_indv*2))
    mv new_lk.txt $baseDIR/lkTable.txt
    lookup_Table=$baseDIR/lkTable.txt
else
    echo -e 'Lookup Table exists!\n'
fi


## Interval mapping
echo -e '\n\n ----- \nInterval mapping\n -----\n'
seqFile=$baseDIR/windows/w$window.$windowStart.$windowEnd.buf$buffer.ldhat.sites
locFile=$baseDIR/windows/w$window.$windowStart.$windowEnd.buf$buffer.ldhat.locs
outPrefix=$baseDIR/windows/w$window.$windowStart.$windowEnd.buf$buffer.bpen$bpen.
niter=10000000
samp=5000
# bpen=5
burnIn=500
echo -e 'seqFile:' $seqFile
echo -e 'locFile:' $locFile
echo -e 'Lookup Table:' $lookup_Table
echo -e 'No. of iterations:' $niter
echo -e 'Sampling iterations:' $samp
echo -e 'bpen:' $bpen
echo -e 'burnIn:' $burnIn
interval -seq $seqFile -loc $locFile -lk $lookup_Table -its $niter -samp $samp -bpen $bpen -prefix $outPrefix


## Summarising results
echo -e '\n\n ----- \nSummarise results\n -----\n'
stat -input ${outPrefix}rates.txt -burn $burnIn -loc $locFile -prefix $outPrefix


## Rhomap
# seqFile=$baseDIR/$chrom/$flank.$chrom.ldhat.sites
# locFile=$baseDIR/$chrom/$flank.$chrom.ldhat.locs
# niter=1000000
# samp=5000
# burn=200000
# rhomap -seq $seqFile -loc $locFile -lk $new_lkFile -its $niter -sampe $samp -burn $burn -prefix $prefix


## Cleanup
# rm $baseDIR/windows/w$window.$windowStart.$windowEnd.buf$buffer.vcf.gz