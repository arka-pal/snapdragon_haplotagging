#!/bin/bash

##### Bash script to change chromosome annotation in files
##### author: Arka Pal
##### written: 25.04.2023
##### update: 28.04.2023

##### NB: changes GWHBJVT00000001 to Chr 1 (and so on...)

##### USAGE: bash change_chromName.sh $1 $2
##### USAGE Example 1: bash change_chromName.sh /nfs/scistore18/bartogrp/apal/snap_hap/bams/v3.5/ema_align/*/*bam bam

## $1 - filetype (bam, vcf)
## $2 - filepath (NB: absolute path)

# Load modules
module load samtools bcftools

# Read input parameters
filetype=$1
filepath=$2



## Changing BAM files ## 
#NB: check the code again before running!

if [ $filetype = 'bam' ]
then
    module load samtools
    
    prefix=${filepath/.$filetype}
    #echo -e '\n\n'
    echo Filename: $prefix.$filetype
    
    # Change bam
    echo Changing BAM
    time samtools view -h $prefix.$filetype | sed 's/GWHBJVT0000000/Chr/g' | samtools view -Shb - -o ${prefix}_new.$filetype
    #mv ${prefix}_new.$filetype $prefix.$filetype
    
    # Index bam
    echo Indexing BAM
    samtools index -@4 ${prefix}_new.$filetype
    
    # Check for linked read file
    if [ -e $prefix.linked_reads.full.bed ]
    then
        # Change linked read file 
        echo Changing linked read file
        cat $prefix.linked_reads.full.bed | sed 's/GWHBJVT0000000/Chr/g' > $prefix.linked_reads.full_new.bed
        #mv $prefix.linked_reads.full_new.bed $prefix.linked_reads.full.bed
    fi
    
    echo 'Done'
fi



if [ $filetype = 'vcf.gz' ]
then
    module load bcftools
    
    echo vcf Filename: $filepath
    
    # Change & index VCF
    echo Changing VCF
    time bcftools view $filepath | sed 's/GWHBJVT0000000/Chr/g' | bcftools view -Oz -o ${filepath/_old_chromName} -
    
    echo Indexing VCF
    tabix ${filepath/_old_chromName}
    
    echo 'Done'
fi