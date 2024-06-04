#!/bin/bash

##### SLURM script to merge BAMs 
##### author: Arka Pal
##### written: 02.05.2022
##### update: 06.07.2023

##### Usage: sbatch --array=1-960 job-merge_bams.sbatch


## Define SLURM variables
#----------------------------------------------------------------

#SBATCH --job-name=mergeBAMs
#SBATCH --output=mergeBAMs_%a_%A.out
#SBATCH --error=mergeBAMs_%a_%A.err

#SBATCH --partition=defaultp
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#----------------------------------------------------------------


## Load modules
module load samtools


## Set variables
WORKDIR=/nfs/scistore18/bartogrp/apal/snap_hap
merged_bamDIR=$WORKDIR/bams/v3.5/bams_merged_2023
bams_to_merge=$WORKDIR/bams/v3.5/sample_list/bams_to_merge.txt

#Read plantID and bamlist
read plantID bamlist < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $bams_to_merge)
if [ $(echo $bamlist | tr " " "\n" | wc -l ) = 3 ];
then
    tmp=$(basename $(echo $bamlist | tr " " "\n" | head -1))
    outBAM=${tmp/.bam/.merged.bam}
else
    tmp=$(basename $(echo $bamlist | tr " " "\n" | tail -1))
    outBAM=${tmp/.bam/.merged.bam}
fi

#Print variables
echo PlantID: $plantID
echo BAMlist: $bamlist
echo outBAM: $outBAM


## Merge BAMs with samtools
time samtools merge -@20 --no-PG -fcp -o ${merged_bamDIR}/$outBAM $bamlist
#time sambamba merge -t20 -p ${merged_bamDIR}/$outBAM $bamlist ## with sambamba


## Fix BAM header
#Extract relevant header from BAM
samtools view -H ${merged_bamDIR}/$outBAM | grep -Ev '@PG|@RG'  > ${merged_bamDIR}/RG_$plantID.txt

#Add read group (@RG)
read_group_ID=${plantID}_$(echo $bamlist | cut -d/ -f10,20,30 | tr "/" "_")
read_group_SM=${tmp/.sorted.BXnum.bam}
read_group_LB=Snapdragon-Haplo_${tmp/_*}
read_group=@RG\\tID:${read_group_ID}\\tSM:${read_group_SM}\\tLB:${read_group_LB}
echo -e $read_group >> ${merged_bamDIR}/RG_${plantID}.txt

#Print BAM header
echo -e "\n BAM header"
cat ${merged_bamDIR}/RG_${plantID}.txt

#Fix BAM header
samtools reheader ${merged_bamDIR}/RG_${plantID}.txt ${merged_bamDIR}/$outBAM > ${merged_bamDIR}/tmp_$plantID.bam

#Rename the re-headered file with outBAM file
mv ${merged_bamDIR}/tmp_$plantID.bam ${merged_bamDIR}/$outBAM

#Remove header file
rm ${merged_bamDIR}/RG_${plantID}.txt


## Index merged BAM file
samtools index -@20 ${merged_bamDIR}/$outBAM