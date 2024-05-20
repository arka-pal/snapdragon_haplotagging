#!/bin/bash

##### Bash script to convert HapCut2 output to a bed file
##### author: Marek Kucka, Frank Chan, Arka Pal
##### date: 16.08.2023
##### update: 16.08.2023

##### USAGE: bash extract_posfile.sh <options>
##### $1 - HapCut2 output file


## Load modules

## Read input
file=$1

## HapCut2 to BED
grep -A1 BLOCK $file | \
	awk '/BLOCK/ {offset=$3;len=$5;phased=$7;span=$9;frags=$11};
		/Chr/ {print $4"\t"$5"\t"$5+span"\t"$4":"$5"-"$5+span":"phased":"span":"frags"\t"offset"\t"len"\t"phased"\t"span"\t"frags}' > ${file/.output/.bed}