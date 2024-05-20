#!/bin/bash

##### SLURM script: Merge molphase combined sample VCF into multi-sample VCF
##### author: Arka Pal
##### written: 17.10.2023
##### update: 03.04.2024

##### USAGE: sbatch -J <job-name> job-molphase-mergeVcfs.sbatch <options>
##### $1 - chrom
##### $2 - baseDIR
##### $3 - sampleVcfList


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------

## not included -----
## SBATCH --exclude=zeta[243-262],bigterra152
## SBATCH --gres=gpu:1
## -----



## Load modules
## -----
module load bcftools/1.18


## Initiate variables
## -----
chrom=$1
baseDIR=$2
sampleVcfList=$3

# chrom=Chr6
# baseDIR=/nfs/scistore18/bartogrp/apal/snap_hap/variants/molphase/${chrom}
# realpath ${baseDIR}/indv-molphase/*/*combined.vcf.gz > ${baseDIR}/molphaseVcfs.list
# sampleVcfList=${baseDIR}/molphaseVcfs.list

# outVCF=${baseDIR}/Am_all_molphased_${chrom}_SnpGap5_biSNPs_filtered-DP_500-7732_QUAL20_MQ30.PL.vcf.gz
outVCF=${baseDIR}/Am_all_stitchRun1_${chrom}.molphased.vcf.gz

echo chrom: $chrom
echo sampleVcfList: $sampleVcfList
echo outVCF: $outVCF


## Merge VCFs
echo -e '\nMerging sample VCFs'
srun time bcftools merge --threads ${SLURM_CPUS_PER_TASK} -Oz -o $outVCF -l $sampleVcfList 
echo -e '\nCreating index'
srun time tabix -f $outVCF


## Create GT file
# srun time bcftools query -f '%CHROM\t%POS[\t%GT]\n' $vcf > Am_all_molphased_Chr6.gt