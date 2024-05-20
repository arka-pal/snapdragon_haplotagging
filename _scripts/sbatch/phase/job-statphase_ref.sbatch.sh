#!/bin/bash

##### SLURM script: Statistical Phasing with ref
##### author: Arka Pal
##### written: 03.05.2024
##### update: 04.05.2024

##### USAGE: sbatch -J <job-name> job-statphase_ref.sbatch <options>
##### $1 - chrom
##### $2 - targetVCF
##### $3 - outputVCF
##### $4 - molVCF

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
#SBATCH --mem-per-cpu=12G

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


## Read input
## -----
chrom=$1
targetVCF=$2
outputVCF=$3
molVCF=$4

## Define variables
baseDIR=~/snap_hap/variants/statphase
run='run-Final'
refPanel=$baseDIR/$run/refPanel/Am_gt5x_refPanel_$chrom.bcf
refScaffold=$baseDIR/$run/refScaffold/Am_gt5x_refScaffold_$chrom.bcf
sampleList=$baseDIR/samples/samples_gt5x.list
cd $baseDIR/$run/$chrom


## Make refPanel
## -----
# echo -e "Subset >5x coverage samples \n"
# time bcftools view -S $sampleList $molVCF -Oz -o ${molVCF/all/gt5x} --write-index --threads ${SLURM_CPUS_PER_TASK}

# echo -e "Making reference panel \n"
# time bcftools filter --set-GTs . -e 'GT~"0/1"' ${molVCF/all/gt5x} -Ob -o $refPanel --write-index --threads ${SLURM_CPUS_PER_TASK}

# echo -e "Making reference scaffold \n"
# time bcftools filter -i 'AN=360' $refPanel -Ob -o $refScaffold --write-index --threads ${SLURM_CPUS_PER_TASK}

echo -e "Comparing no. of  variant sites"
echo -e 'molVCF-gt5x:' $(bcftools view -H ${molVCF/all/gt5x} | wc -l)
echo -e 'refPanel:' $(bcftools view -H $refPanel | wc -l)
echo -e 'refScaffold:' $(bcftools view -H $refScaffold | wc -l)
echo -e "\n"


## Shapeit5
## -----

# With ref Panel
echo -e "Statistical Phasing \n"
phase_common_static --input $targetVCF \
                    --output $outputVCF \
                    --scaffold $refScaffold \
                    --region ${chrom} \
                    --thread ${SLURM_CPUS_PER_TASK}

                    