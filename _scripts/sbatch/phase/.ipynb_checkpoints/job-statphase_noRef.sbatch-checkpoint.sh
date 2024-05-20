#!/bin/bash

##### SLURM script: Statistical Phasing with shapeit5
##### author: Arka Pal
##### written: 03.05.2024
##### update: 03.05.2024

##### USAGE: sbatch -J <job-name> job-statphase_noRef.sbatch <options>
##### $1 - chrom
##### $2 - targetVCF
##### $3 - outputVCF


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x-%A.out
#SBATCH --error=%x-%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --time=240:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G

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
module load shapeit5/5.1.1 bcftools/1.18 vcftools


## Initiate variables
## -----
chrom=$1
targetVCF=$2
outputVCF=$3


## Shapeit5 without ref panel
## -----
phase_common_static --input $targetVCF \
                    --output $outputVCF \
                    --region ${chrom} \
                    --thread 8



# ## Run without reference panel
# targetVCF=~/snap_hap/variants/stitch/${chrom}/Am_all_${stitchRun}_${chrom}.final.vcf.gz
# outputVCF=$baseDIR/$run/$chrom/Am_${samples}_${stitchRun}_${chrom}.noRef.statphased..vcf.gz

# cd $baseDIR/$run/jobs
# sbatch -J statP_${chrom/Chr} ~/snap_hap/_scripts/sbatch/phase/job-statphase_noRef.sbatch.sh $chrom $targetVCF $outputVCF