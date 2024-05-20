#!/bin/bash

##### SLURM script to postprocess output from STITCH imputation
##### author: Arka Pal
##### written: 19.09.2023
##### update: 19.09.2023

##### USAGE: sbatch -J <job-name> job-postprocess_STITCH-imputation.sbatch <options>
## $1 - chrom


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x-%a_%A.out
#SBATCH --error=%x-%a_%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --exclude=zeta[243-262]
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G

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
chrom=$1


## Run postprocessing
srun bash ~/snap_hap/_scripts/bash/stitch/postprocessing_stitch-imputation.sh $chrom