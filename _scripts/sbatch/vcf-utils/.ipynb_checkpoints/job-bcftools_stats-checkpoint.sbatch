#!/bin/bash

##### SLURM script to run bcftools stats
##### author: Arka Pal
##### written: 17.07.2023
##### update: 21.07.2023

##### Usage: sbatch -J <job-name> job-bcftools_stats.sbatch <options>
##### $1: inVCF
##### $2: threads

##### NB: bcftools --exclude=zeta[243-262],beta[231-235]
##### NB: bcftools --include=delta,gamma

## Define SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --exclude=zeta[243-262],beta[231-235]
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=1
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
module load bcftools/1.16 python3


## Read inputs
inVCF=$1
threads=$2

## Print variables
echo inVCF: $inVCF
echo threads: $threads


## Run bcftools stats
echo -e "\nRun bcftools stats\n"
time bcftools stats --threads 10 -s- --af-bins <(seq 0 0.1 1) --depth 0,25000,25 $inVCF > ${inVCF/.vcf.gz/_stats.vchk}

# ## Make stats plots
# echo -e "Make plots\n"
# plot-vcfstats -p ${inVCF/.vcf.gz/_stats} ${inVCF/.vcf.gz/_stats.vchk}

echo "Done"