#!/bin/bash

##### SLURM script to run LDhat and LDhot
##### author: Arka Pal
##### written: 29.04.2024
##### update: 07.06.2024

##### USAGE: sbatch -J <job-name> job-LDhat.sbatch.sh <options>
##### $1 - chrom; Chr6
##### $2 - sampPrefix; MF3.pFR.hCov
##### $3 - sampFile; 


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x-%a-%A.out
#SBATCH --error=%x-%a-%A.out
#SBATCH --open-mode=append
#SBATCH --partition=defaultp
### #SBATCH --constraint=bookworm
### #SBATCH --exclude=bigterra152,zeta[243-262]
#SBATCH --time=240:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G


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
chrom=$1
sampPrefix=$2
sampleFile=$3
ncores=${SLURM_CPUS_PER_TASK}
window=${SLURM_ARRAY_TASK_ID}


## Run LDhat
## -----
srun bash ~/snap_hap/recMap/_scripts/run_LDhat.sh $chrom $sampPrefix $sampleFile $ncores $window