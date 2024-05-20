#!/bin/bash

##### SLURM script to fill AC-AN to VCFs
##### author: Arka Pal
##### written: 18.03.2024
##### update: 18.03.2024

##### USAGE: sbatch -J <job-name> job-postSTITCH-concatVcfs.sbatch <options>
##### $1 - inputVCF


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append

### #SBATCH --partition=defaultp
#SBATCH --partition=debian12
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------


## Load modules

## Read inputs
inputVCF=$1

## Tag VCF with AC, AN
srun bash ~/snap_hap/_scripts/bash/fill-AC-AN_vcfgzFormat.sh $inputVCF ##vcf.gz format
# bash ~/snap_hap/_scripts/bash/fill-AC-AN_bcfFormat.sh $inputVCF ##bcf format