#!/bin/bash

##### USAGE: sbatch -J <job-name> job-change_chromName.sbatch <options>
## $1 - filetype (bam, vcf.gz)
## $2 - filepaths (absolute path to all files; put '' around the file path)

##### USAGE Example 2: sbatch -J pos --array=0-7 ~/snap_hap/scripts/sbatch/job-change_chromName.sbatch 'pos' '/nfs/scistore18/bartogrp/apal/snap_hap/impute/impute_all_mpileup/Chr?/*pos'


## Note: make sure to update variable definitions if new parameters are added  
## Note: $SLURM_ARRAY_TASK_ID is 0-indexed. 


## Defining options for slurm how to run
#----------------------------------------------------------------

#SBATCH --output=%x-%a.out
#SBATCH --error=%x-%a.err

#SBATCH --partition=defaultp
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=128G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#----------------------------------------------------------------

## load any module you need here

## Define variables
filetype=$1
filelist=($2)
filepath=${filelist[$SLURM_ARRAY_TASK_ID]}

## Run commands on SLURM's srun
srun bash ~/snap_hap/scripts/change_chromName.sh $filetype $filepath

## BAM files: /nfs/scistore18/bartogrp/apal/snap_hap/bams/v3.5/ema_align/*/*.bam
## all mpileup VCF files: /nfs/scistore18/bartogrp/apal/snap_hap/impute/impute_*_mpileup/*/*vcf.gz