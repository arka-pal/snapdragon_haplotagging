#!/bin/bash

##### SLURM script to merge STITCH vcfs
##### author: Arka Pal
##### written: 02.10.2023
##### update: 21.11.2023

##### USAGE: sbatch -J <job-name> job-merge_vcfs.sbatch <options>
## $1 - STITCH VCF list
## $2 - outPrefix; e.g., Am_all_stitch_${chrom}_SnpGap5_biSNPs_filtered-DP_500-7732_QUAL20_MQ30

## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --exclude=zeta[243-262]
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------


## Load modules
module load bcftools/1.18 vcftools/0.1.16


## Read inputs
vcf=$1
outPrefix=$2


# ## Merge VCFs
# srun bcftools concat -a -Oz -o $outVcf -f $stitchVcfList
# srun tabix $outVcf
# srun bcftools stats --threads 10 -s- --af-bins <(seq 0 0.1 1) --depth 0,25000,25 $outVcf > ${outVcf/.vcf.gz/_stats.vchk}

## Extract allele frequencies and counts
srun vcftools --gzvcf $vcf --out $outPrefix --freq2
srun vcftools --gzvcf $vcf --out $outPrefix --counts2