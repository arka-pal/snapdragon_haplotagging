#!/bin/bash

##### EMA Alignment Bash Script for Antirrhinum haplotagging samples.
##### This script only for read-alignment, but NOT ema-preprocessing steps.
##### check ema-align.sh for whole EMA alignment to run.

##### written by Arka Pal, last updated 20221209

##### bash usage:  bash ema-align_only.sh <samp_no> <batch> <ref.fa> <ref.version>
##### SLURM script usage:  srun bash ema-align_only.sh $SLURM_ARRAY_TASK_ID $1 $2 $3
##### SLURM job submission:  sbatch --array=1-n job-EMAalign_only.slurm <sequencing_batch> <ref.fa> <ref.version>

## NOTE: change $ema_path acc. to need.

## NOTE: reference genome must be indexed using BWA and samtools.
## samtools faidx ref.fa
## bwa index ref.fa

##### ------



#### -----
#### LOAD ENVIRONMENT & SET VARIABLES
#### -----


## Load softwares & packages via modules
module load bioconda samtools sambamba bwa parallel compression-tools/20220329
# NB: if installed and loaded separately, add to path
# export $PATH=$PATH:</path>

## User Inputs
samp_no=$1 #nth sample from the list of samples provided in $2
batch=$2 #sequencing batch, e.g., TRIO/n96/2x/10x/60x/2xE/10xE/60xE
ref_genome=$3 #path to reference genome file
ref_gen_ver=$4 #reference genome version used


## Set working directories
WORKDIR='/nfs/scistore18/bartogrp/apal/snap_hap'
ema_path=$WORKDIR/bams/v3.5.SL/ema_dir/$batch ## change according to need
samp_folder=$(sed -n "${samp_no}p" <(ls $ema_path))
prefix=${samp_folder%_*}_${ref_gen_ver}

## Set sequence info
library=Snapdragon-Haplo_${batch} #used in BAM header

## Set EMA base directories (create if necessary)
ema_dir=${WORKDIR}/bams/$ref_gen_ver/ema_dir/$batch
if [ ! -e $ema_dir ]; then mkdir -p $ema_dir; fi
ema_align=${WORKDIR}/bams/$ref_gen_ver/ema_align/$batch
if [ ! -e $ema_align ]; then mkdir -p $ema_align; fi

# ## set input variables
read_group=@RG\\tID:$samp_folder\\tSM:${prefix}\\tLB:$library #BAM @RG tag

## Set mapping directory
if [ ! -e "$ema_dir/$prefix" ]; then mkdir -p $ema_dir/$prefix && cd "$_"; else cd $ema_dir/$prefix; fi
if [ ! -e $ema_dir/$prefix/ema-bin ]; then mkdir -p $ema_dir/$prefix/ema-bin; fi

## Check variables
echo Prefix: $prefix
echo EMA bins input directory: $ema_path/$samp_folder
echo EMA working directory: $(pwd)
echo final BAM directory: $ema_align



#### -----
#### ALIGNMENT
#### -----


## EMA Align (sequentially or parallel)
echo "Running EMA on all ema-bin files"
time parallel --bar -j 5 \
   "ema align -t 5 -d -r $ref_genome -p 10x -s {} | \
   samtools sort -@ 5 -m 5G -O bam -l 0 -o ema-bin/{/}.sorted.bam - " ::: $ema_path/$samp_folder/*ema-bin/ema-bin-???

## BWA align for non-BC reads
echo "Running BWA for non-BC file and sort them"
bwa mem -p -t 40 -M -R $read_group $ref_genome $ema_path/$samp_folder/*ema-bin/ema-nobc | samtools sort -@ 4 -O bam -l 0 -o ema-bin/ema-nobc.sorted.bam

## Mark duplicates in BWA alignment 
echo "Marking duplicates in BWA alignment"
sambamba markdup -t 40 -p -l 0 ema-bin/ema-nobc.sorted.bam ema-bin/ema-nobc-pMarkedup.sorted.bam
rm ema-bin/ema-nobc.sorted.bam


#### -----
#### POSTPROCESSING
#### -----


## Merging all bam files
echo "Merging all BAM files"
sambamba merge -t 40 -p $prefix.sorted.bam ema-bin/*.bam

## Change barcodes from base-pairs to AxxBxxCxxDxx
echo "Reconverting 16 base barcodes to haplotag barcodes"
samtools view -h $prefix.sorted.bam |\
awk 'BEGIN {split("AAAT,AAAG,AAAC,AATA,AATT,AATG,AATC,AAGA,AAGT,AAGG,AAGC,AACA,AACT,AACG,AACC,ATAA,ATAT,ATAG,ATAC,ATTA,ATTT,ATTG,ATTC,ATGA,ATGT,ATGG,ATGC,ATCA,ATCT,ATCG,ATCC,AGAA,AGAT,AGAG,AGAC,AGTA,AGTT,AGTG,AGTC,AGGA,AGGT,AGGG,AGGC,AGCA,AGCT,AGCG,AGCC,ACAA,ACAT,ACAG,ACAC,ACTA,ACTT,ACTG,ACTC,ACGA,ACGT,ACGG,ACGC,ACCA,ACCT,ACCG,ACCC,TAAA,TAAT,TAAG,TAAC,TATA,TATT,TATG,TATC,TAGA,TAGT,TAGG,TAGC,TACA,TACT,TACG,TACC,TTAA,TTAT,TTAG,TTAC,TTTA,TTTT,TTTG,TTTC,TTGA,TTGT,TTGG,TTGC,TTCA,TTCT,TTCG,TTCC,TGAA",val,","); \
for(i=1;i<=96;i++){lookup[val[i]]=sprintf("%02d",i)}};/BX:Z:/ {match($0,"BX:Z");bx=substr($0,RSTART,23);out="BX:Z:A"lookup[substr(bx,6,4)]"C"lookup[substr(bx,10,4)]"B"lookup[substr(bx,14,4)]"D"lookup[substr(bx,18,4)]substr(bx,22,2);gsub(bx,out,$0);print $0}; !/BX:Z/' |\
samtools view -@5 - -O BAM -o $prefix.sorted.BXnum.bam

## Index BAM file
echo "Indexing BAM file"
samtools index -@5 $prefix.sorted.BXnum.bam

## Move final BAM and BAI file to output dir
echo "Moving BAM file (sorted, indexed) to final directory: $ema_align"
mv $prefix.sorted.BXnum.bam* $ema_align/

#### -----
echo "Read Mapping COMPLETED! Let's get a beer?"
#### -----
