#!/bin/bash

##### Extract all genotypes based on posterior probability
##### author: Arka Pal
##### written: 09.10.2023
##### update: 20.03.2024

##### USAGE: sbatch -a 1-n -J <job-name> job-postSTITCH-extractGT.sbatch <options>
##### $1 - chrom
##### $2 - stitchRun; 
#####            stitch Run1: stitch
#####            stitch Run2: stitch_run2
#####            stitch Run3: stitch_run3
##### $3 - stitchVcfList (high confidence VCFs) 

##### Chr1: 72 segments
##### Chr2: 78 segments
##### Chr3: 66 segments
##### Chr4: 55 segments
##### Chr5: 72 segments #71 segemnts since the last segment doesn't have any variant sites in it. 
##### Chr6: 56 segments
##### Chr7: 56 segments
##### Chr8: 58 segments



## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x-%a_%A.out
#SBATCH --error=%x-%a_%A.out
#SBATCH --open-mode=append

### #SBATCH --partition=defaultp
#SBATCH --constraint=bookworm
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------


## Load modules
module load bcftools/1.18 vcftools

## Read inputs
chrom=$1
stitchRun=$2
stitchVcfList=$3

## Initiate variables
baseDIR=~/snap_hap/variants/$stitchRun/$chrom
inVCF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $stitchVcfList)
# $(realpath $baseDIR/stitch_chromSegments/*/stitch.${chrom}.*.*[^PL].vcf.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")

echo -e "\n"
echo chrom: $chrom
echo baseDIR: $baseDIR 
echo stitchVcfList: $stitchVcfList
echo inVCF: $inVCF
echo -e '\n'

## Run postprocessing
echo -e "Extracting GT from highConf stitch VCF"
time srun bash ~/snap_hap/_scripts/bash/stitch/add_PG_PL.sh $inVCF