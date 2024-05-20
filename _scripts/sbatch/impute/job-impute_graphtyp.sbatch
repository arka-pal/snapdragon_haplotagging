#!/bin/bash

##### SLURM script to submit jobs for impute SNPs by graphtyper in high-coverage samples 
##### author: Arka Pal
##### written: Feb 2023
##### update: 29.04.2023

##### Usage: sbatch -J ChrX job-graphtyper.sbatch <region>
##### <region>: Chr{1..8}


## Defining options for slurm how to run
#----------------------------------------------------------------
#SBATCH --job-name=%x_graphtyper
#SBATCH --output=%x_graphtyper_%A.out
#SBATCH --error=%x_graphtyper_%A.err

#SBATCH --partition=gpu
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=21
#SBATCH --mem-per-cpu=16G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------


## load any module you need here
module load bioconda


## Define variables
WORKDIR=/nfs/scistore18/bartogrp/apal/snap_hap
outDIR=$WORKDIR/impute/impute_highCov_graphtyp
ref_genome=$WORKDIR/ref_genome/v3.5/Amajus_v3.5.fa


## run commands on SLURM's srun
srun graphtyper genotype $ref_genome \
                    --sams=$outDIR/bamlist/bams_10xM_60xM.txt \
                    --region=Chr1 \
                    --get_sample_names_from_filename
                    --sams_index=$outDIR/bamlist/bams_10xM_60xM_index.txt \
                    #--output=$outDIR \
                    #--log=$outDIR/logfile/Chr1.log \
                    #--verbose \
                    #--threads=21
                    #--avg_cov_by_readlen=$outDIR/bamlist/bams_avg_cov_by_readlen.txt \
                    
#srun bcftools concat --naive -Oz -o $outDIR/highCov_graphtyp_VCFs/highCov_$1.vcf.gz $outDIR/results/$1/*vcf.gz

