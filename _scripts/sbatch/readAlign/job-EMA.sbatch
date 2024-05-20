#!/bin/bash

##### SLURM script to run read alignment with EMA 
##### author: Arka Pal
##### written: 02.05.2022
##### update: 03.07.2023

##### Usage: sbatch --array=1-n job-align.slurm $1 $2 $3 $4 $5 $6 $7
##### $1: <path/to/sample_list.txt> 
##### $2: <path/to/fastq_directory> 
##### $3: <sequencing_batch> 
##### $4: <sampling_location> 
##### $5: <species>
##### $6: <path/to/ref.fa>
##### $7: <ref.version>


## Define SLURM variables
#----------------------------------------------------------------

#SBATCH --job-name=ema-align
#SBATCH --output=job-ema-%a_%A.out
#SBATCH --error=job-ema-%a_%A.err
#SBATCH --partition=defaultp
#SBATCH --exclude=beta[231-235]
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
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
module load bioconda samtools sambamba bwa parallel compression-tools/20220329

## Run read alignment
srun bash /nfs/scistore18/bartogrp/apal/snap_hap/_scripts/bash/ema.sh $SLURM_ARRAY_TASK_ID $1 $2 $3 $4 $5 $6 $7