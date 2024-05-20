#!/bin/bash

##### SLURM script to run STITCH imputation (Run 2)
##### author: Arka Pal
##### written: 30.01.2024
##### update: 30.01.2024

##### USAGE: sbatch --array=1-n -J <job-name> job-impute_stitch.sbatch <options>
##### $1 - chrom
##### $2 - buffer
##### $3 - K #75
##### $4 - downsampleToCov #20
##### $5 - use_bx_tag (TRUE/FALSE) #TRUE
##### $6 - ngen #100
##### $7 - niter #40
##### $8 - expRate #0.5
##### $9 - plot (TRUE/FALSE) #TRUE


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%a-stitch-%A.out
#SBATCH --error=%x_%a-stitch-%A.out
#SBATCH --open-mode=append
#SBATCH --partition=defaultp
### #SBATCH --constraint=bookworm
#SBATCH --exclude=bigterra152,zeta[243-262]
#SBATCH --time=240:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=30G

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
module load stitch/1.6.10 bcftools/1.18


## Initiate variables
## -----

# General
vcfDIR=~/snap_hap/variants/vcf_bcftools_Am_all
stitchDIR=~/snap_hap/variants/stitch
chrom=$1
chromSegments=~/snap_hap/ref_genome/chromSegments/${chrom}_segments.txt
#SLURM_ARRAY_TASK_ID=1

# posfile variables
chromStart=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $chromSegments | cut -f1)
chromEnd=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $chromSegments | cut -f2)
buffer=$2
inVCF=$vcfDIR/$chrom/Am_all_bcftools_${chrom}_SnpGap5_biSNPs_filtered-DP_500-7732_QUAL20_MQ30.vcf.gz
posfile=$stitchDIR/${chrom}/posfile_chromSegments/${chrom}_${chromStart}-${chromEnd}_buffer$buffer.pos


# STITCH variables
K=$3
downsampleToCov=$4
use_bx_tag=$5
ngen=$6
niter=$7
expRate=$8
plot=$9
bamlist=~/snap_hap/sample_info/bam_info/bams_Am_all.txt
outputDIR=~/snap_hap/variants/stitch_run2/$chrom/stitch_chromSegments/$(basename ${posfile/.pos})_K${K}_cov${downsampleToCov}_bxTRUE_niter${niter}_ngen${ngen}_r${expRate}_plotFALSE
if [ ! -d  ~/snap_hap/variants/stitch_run2/$chrom/stitch_chromSegments ]; then mkdir -p ~/snap_hap/variants/stitch_run2/$chrom/stitch_chromSegments; fi
if [ ! -d  $outputDIR ]; then mkdir -p $outputDIR; fi


## Create posfile
## ----- 
echo -e "\n$chrom $chromStart $chromEnd"
echo inVCF: $inVCF

# echo Creating posfile: $posfile
# if [ ! -e  $posfile ]; 
# then 
# time bash ~/snap_hap/_scripts/bash/stitch/extract_posfile.sh $chrom $chromStart $chromEnd $buffer $inVCF $posfile
# echo -e '\nposfile created\n'
# else
# echo -e '\nposfile exists\n';
# fi


## Run STITCH
## -----
echo bamlist: $bamlist
echo posfile: $posfile
echo OUTPUT directory: $outputDIR
echo Region: $chrom $chromStart $chromEnd $buffer

srun STITCH.R --chr=$chrom \
    --K=$K \
    --downsampleToCov=$downsampleToCov \
    --use_bx_tag=${use_bx_tag} \
    --niterations=$niter \
    --expRate=0.5  \
    --nGen=$ngen \
    --plotAfterImputation=$plot \
    --plotHapSumDuringIterations=$plot \
    --plot_shuffle_haplotype_attempts=$plot \
    --regionStart=$chromStart \
    --regionEnd=$chromEnd \
    --buffer=$buffer \
    --bamlist=$bamlist \
    --posfile=$posfile \
    --outputdir=$outputDIR \
    --nCores=${SLURM_CPUS_PER_TASK}


## Unused Options 
## --regenerateInput=FALSE \
## --originalRegionName=${chrom}.${chromStart}.${chromEnd} 
