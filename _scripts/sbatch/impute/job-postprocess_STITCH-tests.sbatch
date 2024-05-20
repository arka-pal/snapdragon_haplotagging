#!/bin/bash

##### SLURM script to postprocess output from stitch validation tests
##### author: Arka Pal
##### written: 19.09.2023
##### update: 19.09.2023

##### USAGE: sbatch job-postprocess_STITCHrun.sbatch <options>
## $1 - stitch output folder list
## $2 - KASP coords


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=postSTITCHrun-%a_%A.out
#SBATCH --error=postSTITCHrun-%a_%A.out
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
module load bcftools/1.16 vcftools/20210912


## Read inputs
stitchFolderList=$1
KASPcoords=$2


## Initiate variables
stitchFolder=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $stitchFolderList)
chrom=`echo $stitchFolder | cut -f9 -d/` 
snpID=`basename $stitchFolder | cut -f1 -d- | cut -f2 -d_` 
coord=`cat $KASPcoords | grep -w ${chrom/Chr}_${snpID} | cut -f8`
        
echo -e '\n'$chrom $snpID $coord 
echo $stitchFolder


## Run postprocessing
cd $stitchFolder
bash ~/snap_hap/_scripts/bash/stitch/postprocessing_stitch-tests.sh $chrom $snpID $coord ./stitch.$chrom.*.*[^PL].vcf.gz