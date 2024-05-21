#!/bin/bash

##### SLURM script to run LDhat and LDhot
##### author: Arka Pal
##### written: 29.04.2024
##### update: 16.05.2024

##### USAGE: sbatch -J <job-name> job-LDhat.sbatch.sh <options>
##### $1: flank
##### $2: chrom

## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x-ldhat-%A.out
#SBATCH --error=%x-ldhat-%A.out
#SBATCH --open-mode=append
#SBATCH --partition=defaultp
### #SBATCH --constraint=bookworm
### #SBATCH --exclude=bigterra152,zeta[243-262]
#SBATCH --time=240:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
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
module load bcftools vcftools ldhat

## Initiate variables
## -----
flank=$1
chrom=$2

## Run LDhat
## -----
bash ~/snap_hap/recMap/_scripts/run_LDhat.sh $flank $chrom