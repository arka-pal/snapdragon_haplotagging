#!/bin/bash

##### Bash Script to postprocess STITCH output
##### Author: Arka Pal
##### Date: 28.09.2023
##### last Update: 09.10.2023

##### USAGE: bash postprocessing_stitch-imputation.sh <options>
##### $1 - chrom 


## Load modules
module load bcftools/1.18 vcftools


## initialize variables
chrom=$1
baseDIR=~/snap_hap/variants/stitch/${chrom}
outVCF=${baseDIR}/Am_all_stitch_${chrom}_SnpGap5_biSNPs_filtered-DP_500-7732_QUAL20_MQ30.vcf.gz

## Make stitch-VCF list
realpath ${baseDIR}/stitch_chromSegments/*/stitch.${chrom}.*.*[^PL].vcf.gz > ${baseDIR}/${chrom}-highConf_stitchVcfs.list
stitchVcfList=${baseDIR}/${chrom}-highConf_stitchVcfs.list

echo -e "\n"
echo chrom: $chrom
echo baseDIR: $baseDIR 
echo stitchVcfList: $stitchVcfList
echo outVCF: $outVCF
echo -e '\n'


## Index all stitch VCFs
echo -e '\nIndex all stitch VCFs'
parallel -j10 "bcftools index --tbi -f {}" ::: ${baseDIR}/stitch_chromSegments/*/stitch.${chrom}.*.*[^PL].vcf.gz


## Merge all stitch VCFs
echo -e '\nMerge all stitch VCFs'
time bcftools concat --allow-overlaps -Oz -o $outVcf -f $stitchVcfList


## Extract all genotypes from highconfidence STITCH vcf based on posterior genotype probabilities.
echo -e '\nExtract all genotypes from STITCH'
time bash ~/snap_hap/_scripts/bash/stitch/add_PG_PL.sh $outVCF


## Make VCF stats file
echo -e '\nbcftools stats'
time bcftools stats --threads 10 -s- --af-bins <(seq 0 0.1 1) --depth 0,25000,25 $outVcf > ${outVcf/.vcf.gz/.stats}
time bcftools stats --threads 10 -s- --af-bins <(seq 0 0.1 1) --depth 0,25000,25 ${outVcf/.vcf.gz/.PL.vcf.gz} > ${outVcf/.vcf.gz/.PL.stats}


## Extract individualised statistics
echo -e '\nCalculating fraction of missing data'
time vcftools --gzvcf $outVCF --missing-indv --out ${chrom}_missing-indv
time vcftools --gzvcf $outVCF --missing-site --out ${chrom}_missing-site

## Calculate heterozygosity per individual
echo -e '\nCalculating heterozygosity'
time vcftools --gzvcf $outVCF --het --out ${chrom}_het-indv