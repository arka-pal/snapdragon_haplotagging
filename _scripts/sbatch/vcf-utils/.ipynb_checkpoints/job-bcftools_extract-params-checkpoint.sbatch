#!/bin/bash

##### SLURM script to extract VCF stats parameters
##### author: Arka Pal
##### written: 07.05.2023
##### update: 24.05.2023

##### USAGE: sbatch -J <job-name> job-bcftools_extract-params.sbatch <options>
##### $1: inVCF

##### NB: bcftools --exclude=zeta[243-262],beta[231-235]
##### NB: bcftools --include=delta,gamma

## Define SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --exclude=zeta[243-262],beta[231-235]
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
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


## load any module you need here
module load bcftools/1.16

## Read inputs
inVCF=$1

## Print variables
echo inVCF: $inVCF

## Extract params
echo -e "\nExtract params from VCF file"
srun time bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%DP\t%AC\t%AN\t%FS\t%MQ\n' $inVCF > ${inVCF/.vcf.gz/_params.txt}

echo -e "\nDone"