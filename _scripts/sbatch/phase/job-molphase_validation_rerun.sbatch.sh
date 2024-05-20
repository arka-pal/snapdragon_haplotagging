#!/bin/bash

##### SLURM script: Rerunning HapCUT2 Phasing validation
##### author: Arka Pal
##### written: 03.10.2023
##### update: 03.10.2023

##### Skips all step except the commands with parameter change
##### USAGE: sbatch --array=1-n -J <job-name> job-molphase_validation_rerun.sbatch <options>
##### $1 - baseDIR
##### $2 - chrom
##### $3 - sampleList
##### $4 - bamList
##### $5 - vcf
##### $6 - threshold


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%a-%A.out
#SBATCH --error=%x_%a-%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --time=120:00:00
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

## not included -----
## SBATCH --exclude=zeta[243-262],bigterra152
## SBATCH --gres=gpu:1
## -----



## Load modules
## -----
module load htslib/1.18 bcftools/1.18 samtools/1.18
#unset $PYTHONPATH
export PYTHONPATH=$PYTHONPATH:~/.local/lib/python3.9/site-packages
export PATH=$PATH:~/_softwares/HapCUT2/build
export PATH=$PATH:~/_softwares/HapCUT2/utilities


## Initiate variables
## -----
baseDIR=$1
chrom=$2
sampleList=$3
bamList=$4
vcf=$5
threshold=$6

sample=$(tail +2 $sampleList | cut -f1 | sed -n "${SLURM_ARRAY_TASK_ID}p")
bam=$(cat $bamList | grep $sample)
coverage=$(tail +2 $sampleList | cut -f2 | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [ ! -d $baseDIR/cov-${coverage}_${chrom}_${sample} ]; then mkdir -p $baseDIR/cov-${coverage}_${chrom}_${sample}; fi
sampleDIR=$baseDIR/cov-${coverage}_${chrom}_${sample}

sampleHetVCF=cov-${coverage}_${chrom}_${sample}
sampleHomVCF=cov-${coverage}_${chrom}_${sample}

unLinkedFragments=cov-${coverage}_${chrom}_${sample}
LinkedFragments=cov-${coverage}_${chrom}_${sample}
outPhased=cov-${coverage}_${chrom}_${sample}_thres${threshold}.hapcut2.output
outPhased_BX=cov-${coverage}_${chrom}_${sample}_thres${threshold}.BX.hapcut2.output


cd $sampleDIR
echo sampleDIR: $sampleDIR
echo chrom: $chrom
echo sample: $sample
echo coverage: $coverage
echo BAM: $bam
echo VCF: $vcf
echo sample HET VCF: $sampleHetVCF
echo hapcut2 output: $outPhased



# ## Preprocessing 
# ## -----

# # Extract sample VCF
# #bcftools view -s $sample $vcf | bcftools view -i "AC=1" -Ov -o ${sampleHetVCF}.het.vcf ## method1
# #time bcftools view -s $sample -r Chr6:1-100000 $vcf | awk '/^#/;/Chr/ {OFS="\t"}; !/^#/ && $10~/^0\/1/' > ${sampleHetVCF}.het.vcf ## test for a truncated region
# time bcftools view -s $sample $vcf | awk '/^#/;/Chr/ {OFS="\t"}; !/^#/ && $10~/^0\/1/' > ${sampleHetVCF}.het.vcf



## HapCUT2
## -----

# # Extract unlinked fragments
# time extractHAIRS --10X 1 --bam $bam --VCF ${sampleHetVCF}.het.vcf --region $chrom --out ${unLinkedFragments}.unlinked.fragments
echo -e '\n'
echo -e "No. of unlinked fragments:" $(wc -l ${unLinkedFragments}.unlinked.fragments)

# Remove non-barcodes
#grep -v "NULL" ${unLinkedFragments}.unlinked.fragments > ${unLinkedFragments}.BX.unlinked.fragments
#echo -e "No. of unlinked BX fragments" $(wc -l ${unLinkedFragments}.BX.unlinked.fragments)


# Convert unlinked to linked fragments
#time python3 ~/_softwares/HapCUT2/utilities/LinkFragments.py  --bam $bam --VCF ${sampleHetVCF}.het.vcf --fragments ${unLinkedFragments}.BX.unlinked.fragments --out ${LinkedFragments}.BX.linked.fragments -d 50000
#echo -e "\nNo. of linked fragments" $(wc -l ${LinkedFragments}.BX.linked.fragments)

# Phase with HapCUT2
time HAPCUT2 --fragments ${LinkedFragments}.linked.fragments --VCF ${sampleHetVCF}.het.vcf --out ${outPhased} --nf 1 --threshold $threshold --error_analysis_mode 1 --call_homozygous 1 --outvcf 1  --v 0




## Postprocessing
## -----
