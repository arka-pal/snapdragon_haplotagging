#!/bin/bash

##### Bash script to extractn50 from molecular phasing output
##### author: Marek Kucka, Frank Chan, Arka Pal
##### date: 16.08.2023
##### update: 16.08.2023

##### USAGE: bash extract_n50.sh <options>
##### $1 - HapCut2 BED file
##### $2 - threshold number; n50


## Load modules
module load datamash/1.8

## Read input
file=$1
num=$2

## Extract n50
if [ -n "$2" ]; then num=$2; else num=50; fi
sum=`sort -nr -k 8,8 $file | datamash sum 8`
threshold=`echo $sum | awk '{print sprintf("%.2f", $1*'$num'/100)}'`
sort -k 8,8nr $file | \
	awk 'NR==1 && $8 >= '$threshold'{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}; 
		{sum+=$8};
		sum < '$threshold'{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' | \
	sort -k 1,1 -k 2,2n > ${file/.bed/.n$num.bed}


sort -k 8,8nr $file | \
	awk 'NR==1 && $7 >= '$threshold'{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}; 
		{sum+=$7};
		sum < '$threshold' {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'