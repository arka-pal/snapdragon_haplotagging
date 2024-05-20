#!/bin/bash

# Usage: sbatch --array=1-n job-align_only.slurm $1 $2 $3
# $1: <sequencing_batch> 
# $2: <path/to/ref.fa>
# $3: <ref.version>


## Defining options for slurm how to run
#----------------------------------------------------------------

#SBATCH --job-name=$SLURM_JOB_NAME
#SBATCH --output=%x-%a-%A.out
#SBATCH --error=%x-%a-%A.err

#SBATCH --partition=gpu
#SBATCH --exclude=gpu230

#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=20G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

## load any module you need here
module load bioconda samtools sambamba bwa parallel compression-tools/20220329
#----------------------------------------------------------------


## run commands on SLURM's srun
srun time bash /nfs/scistore18/bartogrp/apal/snap_hap/scripts/ema-align_only.sh $SLURM_ARRAY_TASK_ID $1 $2 $3