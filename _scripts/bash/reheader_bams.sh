#!/bin/bash

##### Bash script to reheader BAM files (relevant to only mis-named 10xNEW samples while merging)
##### author: Arka Pal
##### written: 28.04.2023
##### update: 06.07.2023

##### USAGE: bash reheader_bams.sh $1 $2 $3
##### $1 - plantID
##### $2 - old sequencing batch i.e., name to be changed
##### $3 - new sequencing batch i.e., changed name

##### Example usage: bash reheader_bams.sh pb1042 10xNEW 60x2


## Load modules
module load samtools


## Set variables
bamDIR=~/snap_hap/bams/v3.5/bams_merged_2023
plantID=$1
bamfile=($bamDIR/*$plantID*bam)
old_batch=$2 #10xE,10xNEW
new_batch=$3 #10x,10x2,60x2 

echo PlantID: $plantID
echo BAMfile: $bamfile


## Create header
samtools view -H $bamfile | grep -vE '@PG' | sed "s/SM:${old_batch}/SM:${new_batch}/" | sed "s/Haplo_${old_batch}/Haplo_${new_batch}/" > $bamDIR/header_$plantID.txt
echo -e "\n\n New Header"
cat $bamDIR/header_$plantID.txt


## Fix BAM header
echo -e "\n\n Fixing BAM Header"
time samtools reheader $bamDIR/header_$plantID.txt $bamfile > $bamDIR/tmp_$plantID.bam


## Rename and index re-headered BAM
mv $bamDIR/tmp_$plantID.bam $bamfile
mv $bamfile ${bamfile/${old_batch}/${new_batch}}
echo -e "\n\n Indexing BAM"
time samtools index ${bamfile/${old_batch}/${new_batch}}


## Clean-up
rm $bamDIR/header_$plantID.txt
rm $bamfile.bai

echo -e "\nDone!"