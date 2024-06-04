#!/bin/bash

##### SLURM script to check BAM quality with BAMqc 
##### author: Arka Pal
##### written: 02.05.2022
##### update: 13.07.2023

##### Usage: sbatch -J <job-name> --array=1-n job-bamqc.sbatch <options>
##### $1: batch
##### $2: bamlist (if indv-bamqc) OR config_file (if multi-bamqc)
##### $3: outDIR

##### job-name: indv-qc OR multi-qc 
##### NB: #SBATCH --mem-per-cpu=32G #for individual bamqc
##### NB: #SBATCH --mem-per-cpu=64G #for multi-bamqc


## Define SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%a_%A.out
#SBATCH --error=%x_%a_%A.err

#SBATCH --partition=defaultp
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


## Load module
module load qualimap java


## Read inputs
batch=$1
sample_list=$2 #bamlist OR config_file
outDIR=$3


## Define variables
WORKDIR=/nfs/scistore18/bartogrp/apal/snap_hap
tmp=$outDIR/$(basename $bamfile)
outDIR_bamqc=${tmp/.bam/_stats}
if [ ! -d $outDIR ]; then mkdir -p $outDIR; fi


## Print variables
echo BAM/config File: $sample_list
#echo BAMqc output Folder: $outDIR_bamqc


# ## Run BAM-QC on SLURM
# bamfile=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $sample_list)
# srun time qualimap bamqc -bam $bamfile \
#         --paint-chromosome-limits \
#         -hm 3 -nr 1000 -nw 400 -nt 24 \
#         --collect-overlap-pairs \
#         -outdir $outDIR_bamqc \
#         -outformat PDF:HTML \
#         --skip-dup-mode 2 \
#         --java-mem-size=32G #only set for 60x samples
#         # --output-genome-coverage $outfolder/coverage.txt
        

## Run multi-BAM-QC on SLURM
config_file=${sample_list}
srun time qualimap multi-bamqc --data $config_file \
        -outdir $outDIR \
        -outformat PDF:HTML \
        --java-mem-size=64G