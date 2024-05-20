#!/bin/bash

##### SLURM script for filtering variants from a VCF 
##### author: Arka Pal
##### written: 18.07.2023
##### update: 21.07.2023

##### Filter1: Remove monomorphic REFs
##### Usage: sbatch -J <job-name> job-bcftools_filter1-removeREFs.sbatch <options>
##### $1: inVCF; Am_all_bcftools_Chr1.vcf.gz
##### $2: outVCF; Am_all_bcftools_Chr1_allSNPs+allINDELs.vcf.gz
##### $3: threads


##### NB: bcftools --exclude=zeta[243-262],beta[231-235]
##### NB: bcftools --include=gamma,delta,serbyn,bigterra

## Define SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --exclude=zeta[243-262],beta[231-235]
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
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
inVCF=$1
outVCF=$2
threads=$2


## Print variables
echo inVCF: $inVCF
echo outVCF: $outVCF
echo threads: $threads


## Remove monomorphic REFs
echo -e "\n Remove monomorphic REF sites"
srun time bcftools view --threads $threads -m2 -Oz -o $outVCF $inVCF

echo -e "\n Index VCF"
srun time bcftools tabix -f $outVCF

echo -e "\nDone"