#!/bin/bash

##### SLURM script to calculate Fst with vcftools
##### author: Arka Pal
##### written: 03.04.3034
##### update: 03.04.3034

##### USAGE: sbatch -J <jobname> job-Fst_vcftools.sbatch.sh <options>
##### 


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x-%A.out
#SBATCH --error=%x-%A.out
#SBATCH --open-mode=append
#SBATCH --partition=defaultp
### #SBATCH --constraint=bookworm
### #SBATCH --constraint=bullseye
#SBATCH --time=240:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

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


## Load modules -----
module load bcftools/1.18 vcftools
# export PATH=$PATH:~/genomics_general:~/genomics_general/VCF_processing

## Read variables -----
baseDIR=$1
chrom=$2
inVCF=$3
pop1=$4
pop2=$5
output=$6


## Define variables -----

stitchRun=stitchRun1

## Run commands -----
srun vcftools --gzvcf $inVCF --weir-fst-pop $pop1 --weir-fst-pop $pop2 --out $output