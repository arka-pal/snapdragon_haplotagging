#!/bin/bash

##### Usage: sbatch --array=1-n job-mol_ident.sbatch <ref_gen_ver>
## ref_gen_ver: v3.5, v3.5.SL, SL 

#SBATCH --job-name=mol_ident
#SBATCH --output=job-mol_ident_%A_%a.out
#SBATCH --error=job-mol_ident_%A_%a.err

#SBATCH --partition=gpu
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=ALL

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

## load any module you need here
module load bioconda samtools perl parallel

## Define variable
ref_gen_ver=$1
WORKDIR='/nfs/scistore18/bartogrp/apal/snap_hap'
bamfolder=$WORKDIR/bams/$ref_gen_ver/ema_align
batch=$(sed -n "${SLURM_ARRAY_TASK_ID}p" <(ls $bamfolder))

## run commands on SLURM's srun
srun parallel --compress --bar -j20 \
    "perl /nfs/scistore18/bartogrp/apal/snap_hap/scripts/bed_write.pl" ::: $bamfolder/$batch/*bam