#!/bin/bash

##### 
##### author: Arka Pal
##### written: 20.03.2024
##### update: 20.03.2024

##### USAGE: sbatch -a 1-n -J <job-name> job-postSTITCH-remove_InvariantSites.sbatch.sh <options>
##### $1 - stitchRun; 
#####            stitch Run1: stitch
#####            stitch Run2: stitch_run2
#####            stitch Run3: stitch_run3
##### $2 - stitchRun out prefix;
#####            stitch Run1: stitchRun1
#####            stitch Run2: stitchRun2
#####            stitch Run3: stitchRun3


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x-%a_%A.out
#SBATCH --error=%x-%a_%A.out
#SBATCH --open-mode=append

### #SBATCH --partition=defaultp
#SBATCH --constraint=bookworm
#SBATCH --time=24:00:00
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


## Load modules
module load bcftools/1.18 vcftools

## Read inputs
chrom=Chr$SLURM_ARRAY_TASK_ID
stitchRun=$1
outPrefixRun=$2

## Initiate variables
baseDIR=/nfs/scistore18/bartogrp/apal/snap_hap/variants/$stitchRun
inVCF=$baseDIR/$chrom/Am_all_stitch_${chrom}_SnpGap5_biSNPs_filtered-DP_500-7732_QUAL20_MQ30.PL.tagged.vcf.gz
outVCF=$baseDIR/$chrom/Am_all_${outPrefixRun}_${chrom}.final.vcf.gz

echo -e "\n"
echo chrom: $chrom
echo baseDIR: $baseDIR 
echo inVCF: $inVCF
echo outVCF: $outVCF
echo -e '\n'

## Run
echo -e "Removing invariant sites"
time srun bcftools view -e 'AC==0 | AC==AN' -Oz -o $outVCF $inVCF
echo -e "Indexing final VCF"
time srun bcftools tabix -f $outVCF