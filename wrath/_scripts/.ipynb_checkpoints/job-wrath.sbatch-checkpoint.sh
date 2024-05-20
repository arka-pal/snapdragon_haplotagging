#!/bin/bash

##### SLURM script wrath
##### author: Arka Pal
##### written: 10.04.2024
##### update: 11.04.2024

##### USAGE: sbatch -J $chrom-wrath job-wrath.sbatch.sh <options> 
##### 


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x-%A.out
#SBATCH --error=%x-%A.out
#SBATCH --open-mode=append
#SBATCH --partition=defaultp
### #SBATCH --constraint=bookworm
### #SBATCH --constraint=bullseye
#SBATCH --time=240:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G

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


## Load modules -----
module load samtools bedtools pysam python3 R 
# export PATH=$PATH:~/genomics_general:~/genomics_general/VCF_processing

## Read variables -----
baseDIR=$1
chrom=$2
start=$3
end=$4
windowSize=$5
bamList=$6


## Define variables -----
refGenome=~/snap_hap/ref_genome/v3.5/Amajus_v3.5.fa

## Run commands -----
srun ~/_softwares/wrath/wrath -g $refGenome -c $chrom -s $start -e $end -w $windowSize  -a $bamList -t ${SLURM_CPUS_PER_TASK} -l -v
