#!/bin/bash

##### EMA Alignment Bash Script
##### written by Arka Pal, last updated 20220317
##### usage: bash ema-align.sh [ref_genome] [input-bam-dir] [input_bam_no]
##### usage (for SLURM job script): srun bash <> $1 $2 $SLURM_ARRAY_TASK_ID
##### usage (for SLURM job submission): sbatch --array=1-247:2 job-HapCall.slurm.sh <ref-genome> <input-bam-dir>  



module load gatk
module load java

ref_genome=$1
input_bam_dir=$2 #eg, ~/Proj-Snap_hap/ema-align/bams-Chr6
input_bam_no=$3 #1-248:2 - put in as a #SLURM_ARRAY_TASK_ID


input_bam=$(ls $input_bam_dir | sed -n "${input_bam_no}p")
output_dir=/nfs/scistore03/bartogrp/apal/Proj-Snap_hap/gatk-HapCall-Chr6


echo Processing file 
gatk HaplotypeCaller -R $ref_genome -I $input_bam_dir/$input_bam -O $output_dir/${input_bam/.bam/.vcf}