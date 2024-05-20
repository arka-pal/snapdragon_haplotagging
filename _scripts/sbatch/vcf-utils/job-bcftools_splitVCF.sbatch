#!/bin/bash

##### SLURM script for splitting VCF into smaller chunks
##### author: Arka Pal
##### written: 17.07.2023
##### update: 17.07.2023

##### Usage: sbatch -J <job-name> job-bcftools_splitVCF.sbatch <options>
##### $1: chromosome
##### $2: chrom_segments; ~/snap_hap/ref_genome/chrom_segments/${chr}_segments.txt
##### $3: vcfDIR; $WORKDIR/variants/vcf_bcftools_Am_all/Chr1
##### $4: vcf; Am_all_bcftools_Chr1.vcf.gz


## Define SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%a_%A.out
#SBATCH --error=%x_%a_%A.err

#SBATCH --partition=defaultp
#SBATCH --exclude=zeta[243-262],beta[231-235]
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
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
chrom_segments=$2
vcfDIR=$3
vcf=$4


## Define variables
WORKDIR=/nfs/scistore18/bartogrp/apal/snap_hap
region=$(awk '{print $1"-"$2}' <(sed -n "${SLURM_ARRAY_TASK_ID}p" $chrom_segments))

outVCF=${vcf/$chr/${chr}_${region}}
outDIR=$vcfDIR/vcf_segments
if [ ! -d $outDIR ]; then mkdir -p $outDIR; fi

## Print variables
echo Region: $chr:$region
echo VCF input: $vcfDIR/$vcf
echo VCF output: $outDIR/outVCF 


## Split VCF into 10MB chunks
echo -e "Splitting VCF: $vcfDIR/$vcf\n"
time bcftools view --threads 4 -r $chr:$region -Oz -o $outDIR/$outVCF $vcfDIR/$vcf
time bcftools tabix $outDIR/$outVCF