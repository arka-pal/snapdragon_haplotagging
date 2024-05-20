#!/bin/bash

##### SLURM script: HapCUT2 molecular phasing
##### author: Arka Pal
##### written: 09.10.2023
##### update: 02.04.2024

##### USAGE: sbatch --array=1-n -J <job-name> job-molphase_allSamples.sbatch <options>
##### $1 - chrom
##### $2 - baseDIR
##### $3 - sampleList
##### $4 - bamList
##### $5 - vcf
##### $6 - threshold


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%a-%A.out
#SBATCH --error=%x_%a-%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------

## not included -----
## SBATCH --exclude=zeta[243-262],bigterra152
## SBATCH --gres=gpu:1
## -----



## Load modules
## -----
module load htslib/1.18 bcftools/1.18 samtools/1.18 datamash/1.8 pysam/0.22
# unset $PYTHONPATH
# export PYTHONPATH=$PYTHONPATH:~/.local/lib/python3.9/site-packages
export PATH=$PATH:~/_softwares/HapCUT2/build
export PATH=$PATH:~/_softwares/HapCUT2/utilities


## Initiate variables
## -----
chrom=$1
baseDIR=$2
sampleList=$3
bamList=$4
vcf=$5
threshold=$6

# chrom=Chr6
# baseDIR=/nfs/scistore18/bartogrp/apal/snap_hap/variants/molphase/${chrom}/indv-molphase
# sampleList=~/snap_hap/variants/molphase/samples_lt5x.list
# bamList=~/snap_hap/sample_info/bam_info/bams_Am_all.txt
# vcf=~/snap_hap/variants/stitch/${chrom}/Am_all_stitch_${chrom}_SnpGap5_biSNPs_filtered-DP_500-7732_QUAL20_MQ30.PL.vcf.gz
# threshold=10
# SLURM_ARRAY_TASK_ID=1

sample=$(cat $sampleList | cut -f4 | sed -n "${SLURM_ARRAY_TASK_ID}p")
bam=$(cat $bamList | grep $sample)

coverage=$(cat $sampleList | cut -f17 | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ ! -d $baseDIR/${chrom}_cov${coverage}_${sample} ]; then mkdir -p $baseDIR/${chrom}_cov${coverage}_${sample}; fi
sampleDIR=$baseDIR/${chrom}_cov${coverage}_${sample}
prefix=${chrom}_cov${coverage}_${sample}
outPhased=${prefix}_thres${threshold}.hapcut2.output
headerFile=~/snap_hap/variants/molphase/molphase.header 

cd $sampleDIR
echo sampleDIR: $sampleDIR
echo chrom: $chrom
echo sample: $sample
# echo coverage: $coverage
echo BAM: $bam
echo VCF: $vcf
echo sample prefix: $prefix
echo hapcut2 output: $outPhased
echo Header File: $headerFile



## Preprocessing 
## -----

# Extract sample VCF
#bcftools view -s $sample $vcf | bcftools view -i "AC=1" -Ov -o ${sampleHetVCF}.het.vcf ## method1
#time bcftools view -s $sample -r Chr6:1-100000 $vcf | awk '/^#/;/Chr/ {OFS="\t"}; !/^#/ && $10~/^0\/1/' > ${sampleHetVCF}.het.vcf ## test for a truncated region
echo -e '\n Extracting sample VCF with heterozygous sites'
srun time bcftools view -s $sample -r $chrom $vcf | awk '/^#/;/Chr/ {OFS="\t"}; !/^#/ && $10~/^0\/1/' > ${prefix}.het.vcf



## HapCUT2
## -----

# Extract unlinked fragments
echo -e '\n Extract unlinked fragments from BAM file'
time extractHAIRS --10X 1 --bam $bam --VCF ${prefix}.het.vcf --region $chrom --out ${prefix}.unlinked.fragments
echo -e "\nNo. of unlinked fragments:" $(wc -l ${prefix}.unlinked.fragments)'\n'

# Convert unlinked to linked fragments
echo -e '\n Convert unlinked to linked fragment'
time python3 ~/_softwares/HapCUT2/utilities/LinkFragments.py  --bam $bam --VCF ${prefix}.het.vcf --fragments ${prefix}.unlinked.fragments --out ${prefix}.linked.fragments -d 50000
echo -e "\nNo. of linked fragments" $(wc -l ${prefix}.linked.fragments)'\n'

# Phase with HapCUT2
echo -e '\n HapCut2 phasing'
time HAPCUT2 --fragments ${prefix}.linked.fragments --VCF ${prefix}.het.vcf --out ${outPhased} --nf 1 --threshold $threshold --error_analysis_mode 1 --call_homozygous 1 --outvcf 1  --v 0
echo -e "\nNo. of total variants" $(bcftools view -H $outPhased.phased.VCF | wc -l)
echo -e "\nNo. of phased variants" $(bcftools view -H $outPhased.phased.VCF | grep -vc "0/1")
echo -e "\nNo. of unphased variants" $(bcftools view -H $outPhased.phased.VCF | grep -c "0/1")



## Postprocessing
## -----

# Convert HapCut2 output to BED file
echo -e '\nConverting HapCut2 output to BED format'
time bash ~/snap_hap/_scripts/bash/molphase/Hapcut2ToBed.sh $outPhased

# Extract n50
echo -e '\nExtract n50'
time bash ~/snap_hap/_scripts/bash/molphase/extract_n50.sh ${outPhased/.output/.bed} 50

# Extract whole sample VCF
echo -e '\nExtract whole sample VCF'
srun time bcftools view -s $sample $vcf | awk '/^#/;/CHROM/ {OFS="\t"}; !/^#/ &&  $10~/^0\/0/ {$10="0|0:"substr($10,5);print $0}; 
									!/^#/ &&  $10~/^1\/1/ {$10="1|1:"substr($10,5);print $0}; 
									!/^#/ && $10~/^0\/1/ {$9=substr($9, 4); 
									$10=substr($10,5);print $0}' > ${prefix}.vcf

# Write annotation file
echo -e '\nWrite annotation file'
srun time bcftools query -f "%CHROM\t%POS[\t%GT\t%PS\t%PQ\t%PD]\n" ${outPhased}.phased.VCF | \
	awk 'BEGIN {OFS = "\t"} {if ($3 == "0/1") $7 = "0"; else $7 = "1"; print}' | \
	bgzip -c > $outPhased.phased.annot.gz
tabix -f -b2 -e2 ${outPhased}.phased.annot.gz

# Combine unphased and phased VCF along with info annotations
echo -e '\nCombine unphased and phased VCF'
srun time bcftools annotate -h $headerFile -a ${outPhased}.phased.annot.gz -c CHROM,POS,FORMAT/GX,FORMAT/PS,FORMAT/PQ,FORMAT/PD,INFO/HAPCUT ${prefix}.vcf | \
	awk '!/<ID=GX/' | \
	sed 's/:GX:/:GT:/' | \
	bcftools view - -Oz -o ${outPhased}.combined.vcf.gz; 
tabix -f ${outPhased}.combined.vcf.gz


## Cleanup
## -----
# rm ${prefix}.het.vcf
# rm ${prefix}.vcf
# rm ${prefix}.*.fragments
# rm ${prefix}.phased.annot*
# rm ${prefix}.phased.VCF