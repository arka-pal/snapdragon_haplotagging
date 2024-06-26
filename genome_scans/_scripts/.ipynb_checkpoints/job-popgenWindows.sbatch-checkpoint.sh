#!/bin/bash

##### SLURM script to get popgen stats in windows
##### author: Arka Pal
##### written: 03.04.3034
##### update: 03.04.3034

##### USAGE: sbatch -J <jobname> job-popgenWindows.sbatch.sh <options>
##### $1 - vcf


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
#SBATCH --cpus-per-task=5
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


## Load modules -----
module load bcftools/1.18 python3
# export PATH=$PATH:~/genomics_general:~/genomics_general/VCF_processing

## Read variables -----
baseDIR=$1
chrom=$2
geno=$3
output=$4
window=$5
popFile=$6
p1=$7
p2=$8
windType=$9

## Define variables -----

stitchRun=stitchRun1

## Run commands -----
srun python ~/genomics_general/popgenWindows.py \
    --windType $windType \
    -g $geno \
    -f phased \
	-o $output \
	-w $window \
	-p $p1 -p $p2 \
	--popsFile $popFile \
	--threads ${SLURM_CPUS_PER_TASK} \
    --writeFailedWindow
