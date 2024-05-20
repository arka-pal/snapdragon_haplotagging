#!/bin/bash

##### EMA Alignment Bash Script for Antirrhinum haplotagging samples.

##### written by Arka Pal, last updated 20221206

##### bash usage:  bash ema-align_v3.5.SL.sh <samp_no> <samp_list> <fast_dir> <batch> <loc> <species>
##### SLURM script usage:  srun bash </path/to/ema-align.sh> $SLURM_ARRAY_TASK_ID $1 $2
##### SLURM job submission:  sbatch --array=1-n </path/to/job-align.slurm.sh> <samp_list_filename> <fastq_directory> <sequencing_batch> <sample_location> 

## NOTE: reference genome must be indexed using BWA. 
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
samp_list=$2 #list of samples, .txt file
fastq_dir=$3 #absolute directory path to fastq files
batch=$4 #sequencing batch, e.g., TRIO/n96/2x/10x/60x/2xE/10xE/60xE
loc=$5 #sampling location
species=$6 #sample species

## Set working directories
WORKDIR='/nfs/scistore18/bartogrp/apal/snap_hap'
samp_list=$WORKDIR/bams/samples/$samp_list

## Set sequence info
library=Snapdragon-Haplo_${batch} #used in BAM header
ref_gen_ver='v3.5.SL' #reference genome version used
ref_genome=${WORKDIR}/reference_genome/$ref_gen_ver/Amajus_v3.5_SLocus.fasta #path to reference genome file

## Set EMA base directories (create if necessary)
ema_dir=${WORKDIR}/bams/$ref_gen_ver/ema_dir/$batch
if [ ! -e $ema_dir ]; then mkdir -p $ema_dir; fi
ema_align=${WORKDIR}/bams/$ref_gen_ver/ema_align/$batch
if [ ! -e $ema_align ]; then mkdir -p $ema_align; fi

## read Sample ID from $samp_list
samp_id=$(sed -n "${samp_no}p" ${samp_list})

## set input variables
read1=$(ls $fastq_dir | grep $samp_id | cut -f1 | head -1) #fastq_R1 read
read2=${read1/_R1/_R2} ##fastq_R2 read
prefix=${batch}-${samp_no}_${species}_${loc}_${samp_id}_${ref_gen_ver} #BAM file prefix
read_group=@RG\\tID:$samp_id\\tSM:${prefix}\\tLB:$library #BAM @RG tag

## Set mapping directory
if [ ! -e "$ema_dir/$prefix" ]; then mkdir -p $ema_dir/$prefix && cd "$_"; else cd $ema_dir/$prefix; fi

## Create symbolic links to fastq reads
ln -s -f ${fastq_dir}${read1} $read1
ln -s -f ${fastq_dir}${read2} $read2

## Check variables
echo Sample list: $samp_list
echo Sample ID: $samp_id
echo Read1: $read1
echo Read2: $read2
echo Prefix: $prefix
echo EMA working directory: $(pwd)
echo final BAM directory: $ema_align

# ##################### START Comment
# #### -----
# #### PREPROCESSING
# #### -----


# ## Translate haplotag barcodes to 16base barcodes
# echo "Translating haplotag barcodes to 16 base barcodes"
# zcat $read1 | sed 's/-N...\tRX:/\tRX:/' | 16BaseBCGen ${samp_id} | bgzip -@ 16 > ${read1/.fastq/.16BCgen.fastq}
# ln -s -f $fastq_dir/$read2 ${read2/.fastq/.16BCgen.fastq}
# cut -f 2 ${samp_id}_HaploTag_to_16BaseBCs | tail +2 > ${samp_id}_HaploTag_to_16BaseBCs.ema

# ## EMA Count
# echo "Running EMA count"
# paste <(pigz -c -d ${read1/.fastq/.16BCgen.fastq} | paste - - - - | awk '{print $1"\t"$5"\t"$6"\t"$7}') <(pigz -c -d ${read2/.fastq/.16BCgen.fastq} | paste - - - - | awk '{print $1"\t"$5"\t"$6"\t"$7}' ) |\
# tr "\t" "\n" |\
# ema count -w ${samp_id}_HaploTag_to_16BaseBCs.ema -o $samp_id.16BCgen 2>$samp_id.16BCgen.log
# cat $samp_id.16BCgen.log

# ## EMA Preprocessing
# echo "Running EMA preproc"
# paste <(pigz -c -d ${read1/.fastq/.16BCgen.fastq} | paste - - - - | awk '{print $1"\t"$5"\t"$6"\t"$7}') <(pigz -c -d ${read2/.fastq/.16BCgen.fastq} | paste - - - - | awk '{print $1"\t"$5"\t"$6"\t"$7}' ) |\
# tr "\t" "\n" |\
# ema preproc -w ${samp_id}_HaploTag_to_16BaseBCs.ema -n 500 -t 40 -o ${samp_id}_ema-bin ${samp_id}.16BCgen.ema-ncnt 2>&1 |\
# tee ${samp_id}_preproc.log


# #### -----
# #### ALIGNMENT
# #### -----


# ## EMA Align (sequentially or parallel)
# echo "Running EMA on all ema-bin files"
# #for file in ${samp_id}_ema-bin/ema-bin-???; 
# #do echo "Processing $file"; ema align -t 4 -d -r $ref_genome -R $read_group -p 10x -s $file | samtools sort -@ 4 -O bam -l 0 -m 4G -o ${file}.sorted.bam;
# #done
# parallel --bar -j10 \
#     "ema align -t 4 -d -r $ref_genome -p 10x -s {} | \
#     samtools sort -@ 4 -O bam -l 0 -m 4G -o {}.sorted.bam - " ::: ${samp_id}_ema-bin/ema-bin-???

# ## BWA align for non-BC reads
# echo "Running BWA for non-BC file and sort them"
# bwa mem -p -t 40 -M -R $read_group $ref_genome ${samp_id}_ema-bin/ema-nobc | samtools sort -@ 4 -O bam -l 0 -m 4G -o ${samp_id}_ema-bin/ema-nobc.sorted.bam

# ## Mark duplicates in BWA alignment 
# echo "Marking duplicates in BWA alignment"
# sambamba markdup -t 40 -p -l 0 ${samp_id}_ema-bin/ema-nobc.sorted.bam ${samp_id}_ema-bin/ema-nobc-pMarkedup.sorted.bam
# rm ${samp_id}_ema-bin/ema-nobc.sorted.bam
# ##################### START Comment

#### -----
#### POSTPROCESSING
#### -----


## Merging all bam files
echo "Merging all BAM files"
sambamba merge -t 40 -p $prefix.sorted.bam ${samp_id}_ema-bin/*.bam

## Change barcodes from base-pairs to AxxBxxCxxDxx
echo "Reconverting 16 base barcodes to haplotag barcodes"
samtools view -h $prefix.sorted.bam |\
awk 'BEGIN {split("AAAT,AAAG,AAAC,AATA,AATT,AATG,AATC,AAGA,AAGT,AAGG,AAGC,AACA,AACT,AACG,AACC,ATAA,ATAT,ATAG,ATAC,ATTA,ATTT,ATTG,ATTC,ATGA,ATGT,ATGG,ATGC,ATCA,ATCT,ATCG,ATCC,AGAA,AGAT,AGAG,AGAC,AGTA,AGTT,AGTG,AGTC,AGGA,AGGT,AGGG,AGGC,AGCA,AGCT,AGCG,AGCC,ACAA,ACAT,ACAG,ACAC,ACTA,ACTT,ACTG,ACTC,ACGA,ACGT,ACGG,ACGC,ACCA,ACCT,ACCG,ACCC,TAAA,TAAT,TAAG,TAAC,TATA,TATT,TATG,TATC,TAGA,TAGT,TAGG,TAGC,TACA,TACT,TACG,TACC,TTAA,TTAT,TTAG,TTAC,TTTA,TTTT,TTTG,TTTC,TTGA,TTGT,TTGG,TTGC,TTCA,TTCT,TTCG,TTCC,TGAA",val,","); \
for(i=1;i<=96;i++){lookup[val[i]]=sprintf("%02d",i)}};/BX:Z:/ {match($0,"BX:Z");bx=substr($0,RSTART,23);out="BX:Z:A"lookup[substr(bx,6,4)]"C"lookup[substr(bx,10,4)]"B"lookup[substr(bx,14,4)]"D"lookup[substr(bx,18,4)]substr(bx,22,2);gsub(bx,out,$0);print $0}; !/BX:Z/' |\
samtools view -@ 4 - -O BAM -o $prefix.sorted.BXnum.bam

## Index BAM file
echo "Indexing BAM file"
samtools index -@ 4 $prefix.sorted.BXnum.bam

## Move final BAM and BAI file to output dir
echo "Moving BAM file (sorted, indexed) to final directory: $ema_align"
mv $prefix.sorted.BXnum.bam* $ema_align/


echo "Read Mapping COMPLETED!"