#!/bin/bash

##### SLURM script for changing sample names in VCF tools
##### author: Arka Pal
##### written: 14.07.2023
##### update: 14.07.2023

##### Usage: sbatch -J <job-name> job-bcftools_reheaderVCF.sbatch <options>
##### $1: chr; Chr1
##### $2: sample_names; $WORKDIR/sample_info/samples_SnapHap_LastUpdate-2023-07.txt
##### $3: vcfDIR; ~/snap_hap/variants/vcf_bcftools_Am_all/Chr1
##### $4: outVCF=Am_all_bcftools_Chr1.vcf.gz


## Define SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err

#SBATCH --partition=defaultp
#SBATCH --exclude=zeta[243-262],beta[231-235]
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------


## Load modules
module load bcftools/1.16


## Read inputs
chr=$1
sample_names=$2
vcfDIR=$3
vcf=$4


## Define variables
WORKDIR=/nfs/scistore18/bartogrp/apal/snap_hap


## Print variables
echo chr: $chr
echo sample_names Files: $sample_names
echo vcfDIR: $vcfDIR
echo vcf: $vcf


## Run bcftools reheader
echo Changing sample names 
time bcftools reheader -s <(cut -f3,4 ${sample_names}) -o $vcfDIR/tmp_newSamples.vcf.gz  $vcfDIR/$vcf
mv $vcfDIR/tmp_newSamples.vcf.gz $vcfDIR/$vcf

echo Index VCF
time bcftools tabix -f $vcfDIR/$vcf