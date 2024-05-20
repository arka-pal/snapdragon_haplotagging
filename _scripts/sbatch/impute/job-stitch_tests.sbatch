#!/bin/bash

##### SLURM script to test imputation accuracy with different STITCH parameters
##### author: Arka Pal
##### written: 02.05.2023
##### update: 13.09.2023

##### USAGE: sbatch -J <job-name> job-impute_stitch.sbatch <options>
## $1 - chr
## $2 - K
## $3 - posfile
## $4 - downsampletoCov (default is 50)
## $5 - use_bx_tag (TRUE or FALSE)
## $6 - niterations
## $7 - ngen
## $8 - expected recombination rate (default 0.5cM/MB)

## ONLY IN SPECIAL CASE
## $9 - genfile 
## NB: make sure to update variable definitions if new parameters are added  


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=stitch-%x_%A.out
#SBATCH --error=stitch-%x_%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --exclude=zeta[243-262]
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=16G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------


## Load modules
module load stitch/1.6.7


## Read inputs
chr=$1
K=$2
posfile=$3
downsampleToCov=$4
bx_tag=$5
niter=$6
ngen=$7
expRate=$8


## Define variables

#regionStart=$(($SLURM_ARRAY_TASK_ID*10**6+1))
WORKDIR=/nfs/scistore18/bartogrp/apal/snap_hap/impute/stitch_tests

# bamlist=/nfs/scistore18/bartogrp/apal/snap_hap/sample_info/bam_info/old/bams_Am_all_oldlist.txt ##old list
bamlist=~/snap_hap/sample_info/bam_info/bams_Am_all.txt ##Latest new list after 10xNEW samples

# outputdir=$WORKDIR/$chr/run2_2023-04/${posfile##*/}_K${K}_cov${downsampleToCov}_bxTRUE_niter${niter}_ngen${ngen}_r${expRate} ##old
# outputdir=$WORKDIR/$chr/run3_filteredVCFs_2023-08/run3_${posfile##*/}_K${K}_cov${downsampleToCov}_bxTRUE_niter${niter}_ngen${ngen}_r${expRate}_plotFALSE ##latest iter1
outputdir=$WORKDIR/$chr/run3.1_filteredVCFs_2023-08/run3_${posfile##*/}_K${K}_cov${downsampleToCov}_bxTRUE_niter${niter}_ngen${ngen}_r${expRate}_plotFALSE ##latest iter2

regionStart=$(head -1 $posfile | cut -f2)
regionEnd=$(tail -1 $posfile | cut -f2)

echo bamlist: $bamlist
echo posfile: $posfile
echo OUTPUT directory: $outputdir
echo Region: $chr $regionStart $regionEnd 


## Run stitch without genfile
# xvfb-run STITCH.R --chr=$chr \
srun STITCH.R --chr=$chr \
    --nGen=$ngen \
    --K=$K \
    --posfile=$posfile \
    --regionStart=$regionStart \
    --regionEnd=$regionEnd \
    --bamlist=$bamlist \
    --outputdir=$outputdir \
    --downsampleToCov=$downsampleToCov \
    --buffer=5000 \
    --niterations=$niter \
    --nCores=10 \
    --use_bx_tag=TRUE \
    --plotAfterImputation=FALSE \
    --plotHapSumDuringIterations=FALSE \
    --plot_shuffle_haplotype_attempts=FALSE \
    --expRate=0.5 \
    --minRate=0.1 \
    --maxRate=100



## special cases
#   --expRate=10 \
#   --minRate=2 \
#   --maxRate=2000
    
## normal case
#    --expRate=0.5 \
#    --minRate=0.1 \
#    --maxRate=100

## special cases
#    --genfile=$11\
#    --regenerateInput=FALSE \
#    --originalRegionName=GWHBJVT00000006.1.1000 \

## Only generate input files for STITCH runs
# srun STITCH.R --chr=GWHBJVT00000006 \
#     --posfile=pos_chr6.txt \
#     --bamlist=bamlist.txt \
#     --K=2 \
#     --nGen=1000 \
#     --outputdir=./input_files/chr6_$regionStart-$((regionStart+999999)) \
#     --genfile=gen_chr6.txt\
#     --regionStart=$regionStart \
#     --regionEnd=$((regionStart+999999)) \
#     --buffer=100000 \
#     --nCores=10 \
#     --generateInputOnly=TRUE \
#     --use_bx_tag=TRUE