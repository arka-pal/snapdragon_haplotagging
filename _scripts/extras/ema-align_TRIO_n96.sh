#!/bin/bash

##### EMA Alignment Bash Script
##### written by Arka Pal, last updated 20220314
##### usage: bash ema-align.sh [sample_info_file] [line_number_in_samp_info_file]
##### usage (for SLURM job script): srun bash </path/to/ema-align.sh> $1 $SLURM_ARRAY_TASK_ID
##### usage (for SLURM job submission): sbatch --array=1-96 </path/to/job-align.slurm.sh> </path/to/sample_info> 

## load modules
echo Loading required modules
module load bioconda
module load samtools
module load sambamba
module load bwa

## MANUAL: Set working directories
fastq_dir='/nfs/scistore03/bartogrp/apal/snapdragon_genomics/haplotag_n96_2021/fastq_split_cutadapt' #folder n96 fastq files
#fastq_dir='/nfs/scistore03/bartogrp/apal/snapdragon_genomics/haplotag_TRIO_n28_2020/fastq' #folder TRIO n28 fastq files
ema_dir='/nfs/scistore03/bartogrp/apal/Proj-Snap_hap/ema-dir' #ema processing folder for each sample
output_dir='/nfs/scistore03/bartogrp/apal/Proj-Snap_hap/ema-align' #output folder with BAM and BAI files

## MANUAL: Set seq info
library='Snapdragon-Haplo-n96-2020'
#library='Snapdragon-Haplo-TRIO-n28-2020' #sequencing library TRIO n28
seq_platform='Illumina.NovaSeq_S4.2x151'
#seq_platform='Illumina.HiSeq3000.2x150' #sequencing platform TRIO n28
ref_gen_ver='v3' #reference genome version used
ref_genome='/nfs/scistore03/bartogrp/apal/snapdragon_genomics/reference_genome_V3/Amajus.IGDBv3.chr.fasta' #path to reference genome file

## Reading inputs
echo Reading inputs
samp_info=$1 #file: sample ID, barcode, Sampling location, Species
samp_no=$2 #line no in samp_info file

## set variables
echo Setting variables
read1=$(ls $fastq_dir/ | grep $(sed -n "${samp_no}p" ${samp_info} | cut -f 1) | head -1) #eg, x4593_C03_R1_001.fastq.gz
read2=${read1/_R1/_R2} #eg, x4593_C03_R2_001.fastq.gz
samp_loc=$(sed -n "${samp_no}p" ${samp_info} | cut -f 3) #eg, Ave
samp_sp=$(sed -n "${samp_no}p" ${samp_info} | cut -f 4) #eg, Amajus
prefix=$(sed -n "${samp_no}p" ${samp_info} | awk '{print $1"_"$2}') #eg, x4593_C03
#prefix=TRIO_$(sed -n "${samp_no}p" ${samp_info} | awk '{print $1"_"$2}') #eg, TRIO_406_C09  TRIO n28
samp=${prefix}_${samp_loc}_${samp_sp} #eg, x4593_C03_Ave_Amajus OR TRIO_406_C09_GH_Amajus
fbname=${samp}_${ref_gen_ver}.sorted #eg, x4593_C03_Ave_Amajus_v3.sorted
read_group=@RG\\tID:$prefix\\tSM:$samp\\tLB:$library\\tPL:$seq_platform #@RG tag for BAM file withh all above info

## Make relevant directories and symbolic links
mkdir -p $ema_dir/$samp && cd "$_"
ln -s -f $fastq_dir/$read1 $read1
ln -s -f $fastq_dir/$read2 $read2

## Translate haplotag barcodes to 16base barcodes
echo Translating haplotag barcodes to 16 base barcodes
zcat $read1 | sed 's/AC\tRX:/\tRX:/' | 16BaseBCGen ${prefix} | bgzip -@ 16 > ${read1/.fastq/.16BCgen.fastq} #for n96 files
#zcat $read1 | 16BaseBCGen ${prefix} | bgzip -@ 16 > ${read1/.fastq/.16BCgen.fastq} #for TRIO n28 files
ln -s -f $(pwd)/$read2 ${read2/.fastq/.16BCgen.fastq}
cut -f 2 ${prefix}_HaploTag_to_16BaseBCs | tail +2 > ${prefix}_HaploTag_to_16BaseBCs.ema

## EMA Count
echo Running ema count
paste <(pigz -c -d ${read1/.fastq/.16BCgen.fastq} | paste - - - - | awk '{print $1"\t"$5"\t"$6"\t"$7}') <(pigz -c -d ${read2/.fastq/.16BCgen.fastq} | paste - - - - | awk '{print $1"\t"$5"\t"$6"\t"$7}' ) |\
tr "\t" "\n" |\
ema count -w ${prefix}_HaploTag_to_16BaseBCs.ema -o $prefix.16BCgen 2>$prefix.16BCgen.log
cat $prefix.16BCgen.log

## EMA Preprocessing
echo Running ema preproc
paste <(pigz -c -d ${read1/.fastq/.16BCgen.fastq} | paste - - - - | awk '{print $1"\t"$5"\t"$6"\t"$7}') <(pigz -c -d ${read2/.fastq/.16BCgen.fastq} | paste - - - - | awk '{print $1"\t"$5"\t"$6"\t"$7}' ) |\
tr "\t" "\n" |\
ema preproc -w ${prefix}_HaploTag_to_16BaseBCs.ema -n 500 -t 40 -o ${prefix}_ema-bin ${prefix}.16BCgen.ema-ncnt 2>&1 |\
tee ${prefix}_preproc.log

## EMA Align
echo Running ema align 
for file in ${prefix}_ema-bin/ema-bin-???; 
do echo "Processing $file"; ema align -t 4 -d -r $ref_genome -R $read_group -p 10x -s $file | samtools sort -@ 4 -O bam -l 0 -m 4G -o ${file}.sorted.bam;
done

## BWA align for non-BC reads
echo Running BWA for non-BC file
bwa mem -p -t 40 -M -R $read_group $ref_genome ${prefix}_ema-bin/ema-nobc | samtools sort -@ 4 -O bam -l 0 -m 4G -o ${prefix}_ema-bin/ema-nobc.sorted.bam

## Mark duplicates in BWA alignment 
echo Marking duplicates in BWA alignment
sambamba markdup -t 40 -p -l 0 ${prefix}_ema-bin/ema-nobc.sorted.bam ${prefix}_ema-bin/ema-nobc-pMarkedup.sorted.bam
rm ${prefix}_ema-bin/ema-nobc.sorted.bam

## Merging all bam files
echo Merging all BAM files
sambamba merge -t 40 -p $fbname.bam ${prefix}_ema-bin/*.bam

## Change barcodes from base-pairs to AxxBxxCxxDxx
samtools view -h $fbname.bam |\
awk 'BEGIN {split("AAAT,AAAG,AAAC,AATA,AATT,AATG,AATC,AAGA,AAGT,AAGG,AAGC,AACA,AACT,AACG,AACC,ATAA,ATAT,ATAG,ATAC,ATTA,ATTT,ATTG,ATTC,ATGA,ATGT,ATGG,ATGC,ATCA,ATCT,ATCG,ATCC,AGAA,AGAT,AGAG,AGAC,AGTA,AGTT,AGTG,AGTC,AGGA,AGGT,AGGG,AGGC,AGCA,AGCT,AGCG,AGCC,ACAA,ACAT,ACAG,ACAC,ACTA,ACTT,ACTG,ACTC,ACGA,ACGT,ACGG,ACGC,ACCA,ACCT,ACCG,ACCC,TAAA,TAAT,TAAG,TAAC,TATA,TATT,TATG,TATC,TAGA,TAGT,TAGG,TAGC,TACA,TACT,TACG,TACC,TTAA,TTAT,TTAG,TTAC,TTTA,TTTT,TTTG,TTTC,TTGA,TTGT,TTGG,TTGC,TTCA,TTCT,TTCG,TTCC,TGAA",val,","); \
for(i=1;i<=96;i++){lookup[val[i]]=sprintf("%02d",i)}};/BX:Z:/ {match($0,"BX:Z");bx=substr($0,RSTART,23);out="BX:Z:A"lookup[substr(bx,6,4)]"C"lookup[substr(bx,10,4)]"B"lookup[substr(bx,14,4)]"D"lookup[substr(bx,18,4)]substr(bx,22,2);gsub(bx,out,$0);print $0}; !/BX:Z/' |\
samtools view -@ 4 - -O BAM -o $fbname.BXnum.bam

## Index BAM file
samtools index -@ 4 $fbname.BXnum.bam

## Move final BAM and BAI file to output dir
mv $fbname.BXnum.bam* $output_dir/ && cd "$_"
