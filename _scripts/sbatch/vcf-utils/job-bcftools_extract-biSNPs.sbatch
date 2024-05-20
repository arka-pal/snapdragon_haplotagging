#!/bin/bash

##### SLURM script for extracting biallelic SNPs from VCFs
##### author: Arka Pal
##### written: 17.07.2023
##### update: 18.07.2023

##### Usage: sbatch -J <job-name> job-bcftools_extract-biSNPs.sbatch <options>
##### $1: inVCF; Am_all_bcftools_Chr1.vcf.gz

##### NB: bcftools --exclude=zeta[243-262],beta[231-235]
##### NB: bcftools --include=gamma,delta,serbyn,bigterra

## Define SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err

#SBATCH --partition=defaultp
#SBATCH --exclude=zeta[243-262],beta[231-235]
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=10
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
inVCF=$1


## Define variables
outVCF=${inVCF/.vcf.gz/_biSNPs.vcf.gz}


## Print variables
echo inVCF: $inVCF
echo outVCF: $outVCF


## Extract bi-SNPs VCF
echo -e "Extract bi-allelic SNPs from full VCF\n"
time bcftools view --threads 10 -m2 -M2 -v snps -e "AC==0 || AC==AN" -Oz -o $outVCF $inVCF

## Index VCF
echo -e "Index VCF\n"
time bcftools tabix $outVCF

## Run bcftools stats
echo -e "Run bcftools stats\n"
time bcftools stats --threads 10 -s- --af-bins <(seq 0 0.1 1) --depth 0,25000,25 $outVCF > ${outVCF/.vcf.gz/_stats.vchk}


echo "Done"