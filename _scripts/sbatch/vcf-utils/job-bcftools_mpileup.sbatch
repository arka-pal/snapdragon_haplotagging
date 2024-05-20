#!/bin/bash

##### SLURM script for bcftools pileup and calling variants. 
##### author: Arka Pal
##### written: 10.07.2023
##### update: 10.07.2023

##### Usage: sbatch -J <job-name> job-bamqc.sbatch $1 $2 $3 $4 $5 $6
##### $1: chromosome; Chr1
##### $2: reference genome; $WORKDIR/ref_genome/v3.5/Amajus_v3.5.fa
##### $3: sample_names; $WORKDIR/sample_info/samples_SnapHap_LastUpdate-2023-07.txt
##### $4: bamlist; bamlist=$WORKDIR/sample_info/bam_info/bams_Am_all.txt
##### $5: outDIR=$WORKDIR/variants/vcf_bcftools_Am_all/Chr1
##### $6: outVCF=Am_all_bcftools_Chr1.vcf.gz


## Define SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%a_%A.out
#SBATCH --error=%x_%a_%A.err

#SBATCH --partition=defaultp
#SBATCH --exclude=zeta[243-262],beta[231-235]
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=128G
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
module load bcftools/1.16


## Read inputs
chr=$1
ref_genome=$2
sample_names=$3
bamlist=$4
outDIR=$5
outVCF=$6


## Define variables
WORKDIR=/nfs/scistore18/bartogrp/apal/snap_hap


## Print variables
echo chr: $chr
echo Ref Genome: $ref_genome
echo sample file: $sample_names
echo BAMs: $bamlist
echo output: $outDIR/$outVCF
echo -e "\n"


## Run bcftools mpileup and call variants
echo -e Calling SNPs
srun time bcftools mpileup --annotate AD,ADF,ADR,DP,QS,SP \
        -r $chr \
        -f $ref_genome \
        -Ou -d 500 \
        -S <(cut -f3,4 ${sample_names}) \
        -b $bamlist | \
    bcftools call --annotate GQ,GP \
        -m -Oz \
        -o $outDIR/$outVCF


## Index VCF
echo -e Indexing VCF file
srun time bcftools tabix $outDIR/$outVCF