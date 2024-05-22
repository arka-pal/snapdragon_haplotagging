#!/bin/bash

##### SLURM script: Infer genealogies with singer
##### author: Arka Pal
##### written: 03.05.2024
##### update: 04.05.2024

##### USAGE: sbatch -J <job-name> job-singer.sbatch.sh <options>
##### $1 - ...

## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x-%A.out
#SBATCH --error=%x-%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --time=240:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------

### not included -----
### SBATCH --exclude=zeta[243-262],bigterra152
### SBATCH --gres=gpu:1
### -----


## Load modules
## -----
module load bcftools/1.18 singer


## Read input
## -----
inVCFprefix=$1
outPrefix=$2
chrom=$3
start=$4
end=$5
samples=$6
Ne=$7
mutRate=$8
ratio=$9
mcmc_iters=${10}
thin=${11}


## Define variables
baseDIR=~/snap_hap/ARGs
stitch=stitchRun1
# if [ ! -d $outPrefix ]; then mkdir -p ; fi  


## singer
## -----
singer_master -Ne $Ne \
              -m $mutRate \
              -vcf $inVCFprefix \
              -ratio $ratio \
              -output $outPrefix \
              -start $start \
              -end $end \
              -n $mcmc_iters \
              -thin $thin
