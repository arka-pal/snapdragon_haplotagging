#!/bin/bash

##### Usage: sbatch --array=1-247:2 job-HapCall.slurm.sh <ref-genome> <input-bam-dir>
##### Note: <input-bam-dir> should be provided without last / -> eg, ~/Proj-Snap_hap/ema-align/bams-Chr6


#SBATCH --job-name=gatk-HapCall
#SBATCH --output=job-HapCall.slurm.out
#SBATCH --error=job-HapCall.slurm.err

#SBATCH --partition=defaultp
#SBATCH --time=24:00:00
#SBATCH -c 1
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=ALL

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

## load any module you need here
module load samtools
module load java
module load gatk

## run commands on SLURM's srun
srun bash gatk-HapCall.sh $1 $2 $SLURM_ARRAY_TASK_ID
