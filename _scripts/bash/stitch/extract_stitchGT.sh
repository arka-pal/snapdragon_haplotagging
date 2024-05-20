##### Bash Script to postprocess STITCH output
##### Author: Arka Pal
##### Date: 23.08.2023
##### last Update: 04.09.2023

##### USAGE: bash extract_stitchGT.sh <options>
##### $1 - chrom 
##### $2 - snpID
##### $3 - KASP coord position
##### $4 - stitchVCF

##### stitchGT: text file with columns- <plantID> <stitch> <stitch_highConf>


## Load modules
module load bcftools/1.16 vcftools/20210912

## initialize variables
chrom=$1
snpID=$2
coord=$3
stitchVCF=$4


## Index VCF file
if [ ! -e $stitchVCF.tbi ];
then
    echo -e '\nIndexing highConf-stitch VCF'
    bcftools tabix -f $stitchVCF
fi

## Get 012 values from stitch-highconfidence VCF
echo -e '\nGet 012 values from stitch-highconfidence VCF'
vcftools --gzvcf $stitchVCF --chr $chrom --from-bp $((coord-2)) --to-bp $((coord+2)) --012 --out ${snpID}-${chrom}-${coord}_highConf

## Extract all genotypes from STITCH vcf based on posterior genotype probabilities.
echo -e '\nExtract all genotypes from STITCH'
time bash ~/snap_hap/_scripts/bash/stitch/add_PG_PL.sh $stitchVCF

## Extract 012 values from stitch-ALL VCF
echo -e '\nGet 012 values from stitch-ALL VCF'
vcftools --gzvcf ${stitchVCF/.vcf.gz/.PL.vcf.gz} --chr $chrom --from-bp $((coord-2)) --to-bp $((coord+2)) --012 --out ${snpID}-${chrom}-${coord}_ALL

## Change longFormat sample names to shortFormat
echo -e '\nChange longFormat sample names to shortFormat'
cat ${snpID}-${chrom}-${coord}_ALL.012.indv | sed 's/_v3.5//g' | sed 's/^.*Am.*_.*_//g' > ${snpID}-${chrom}-${coord}.shortFormat.indv

## Create file with stitch genotypes
echo -e '\nCreate stitch GT file'
paste ${snpID}-${chrom}-${coord}.shortFormat.indv <(cut -f2 ${snpID}-${chrom}-${coord}_ALL.012) <(cut -f2 ${snpID}-${chrom}-${coord}_highConf.012) > stitch_${snpID}-${chrom}-${coord}.gt