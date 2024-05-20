#!/bin/bash

##### EMA Alignment Bash Script for 2x/10x/60x coverage files in MEGA batch of samples

##### written by Arka Pal, last updated 20221005
##### usage: bash ema-align.sh <coverage> <plate no.> <line_number_in_samp_info_file>
##### usage (for SLURM job script): srun bash </path/to/ema-align.sh> $SLURM_ARRAY_TASK_ID $1 $2
##### usage (for SLURM job submission): sbatch --array=1-n </path/to/job-align.slurm.sh> <coverage> <plate no.> 



#### LOAD & INITIAITE ENVIRONMENT


## load modules
echo Loading required modules
module load bioconda samtools sambamba bwa compression-tools/20220329


## Inputs
coverage=$1
plate=$2
samp_no=$3


## Set directories
# Parent diretocry
WORKDIR='/nfs/scistore18/bartogrp/apal/snap_hap'
# fastq directory
if [ $coverage = "2x" ];
then fastq_dir=${WORKDIR}/fastq/haplotag_error_MEGA/${coverage}Coverage-e/$plate; #2x;
else fastq_dir=${WORKDIR}/fastq/haplotag_error_MEGA/${coverage}Coverage-e; #10x/60x;
fi
# EMA-mapped reads directory
ema_dir=${WORKDIR}/bams/bams_error_MEGA/ema_dir/${coverage}Coverage/$plate #2x/10x/60x
ema_align=${WORKDIR}/bams/bams_error_MEGA/ema_align/${coverage}Coverage/


## Create directoies if not existing
if [ ! -e $ema_dir ]; then mkdir -p $ema_dir; fi
if [ ! -e $ema_align ]; then mkdir -p $ema_align; fi


## Set sequence info
library='Snapdragon-Haplo-ERROR-n960-2022'
ref_gen_ver='AmajusV3' #reference genome version used
ref_genome=${WORKDIR}/reference_genome/reference_genome_V3/Amajus.IGDBv3.chr.fasta ## path to reference genome file


## Reading inputs
echo Reading inputs
samp_info=$WORKDIR/sample_info/sample_MEGA_n960_2022/samples_${coverage}_${plate/e}.txt #file: sample ID
echo $samp_info
prefix=$(sed -n "${samp_no}p" ${samp_info})
echo $prefix


## set input variables
read1=$(ls $fastq_dir | grep $prefix | cut -f1 | head -1)
read2=${read1/_R1/_R2}
samp=${prefix}_${coverage}_${plate} #eg, pb0112_2x_N701e 
fbname=${samp}_${ref_gen_ver}.sorted #eg, pb0112_2x_N701e_AmajusV3.sorted
read_group=@RG\\tID:$prefix\\tSM:${fbname/.sorted}\\tLB:$library #@RG tag for BAM file with all above info


## Make relevant directories and symbolic links
if [ ! -e "$ema_dir/$samp" ]; then mkdir -p $ema_dir/$samp && cd "$_"; else cd $ema_dir/$samp; fi
ln -s -f $fastq_dir/$read1 $read1
ln -s -f $fastq_dir/$read2 $read2



echo Aligning $samp



#### PREPROCESSING


## Translate haplotag barcodes to 16base barcodes
echo Translating haplotag barcodes to 16 base barcodes
zcat $read1 | sed 's/-N...\tRX:/\tRX:/' | 16BaseBCGen ${samp} | bgzip -@ 16 > ${read1/.fastq/.16BCgen.fastq}
ln -s -f $fastq_dir/$read2 ${read2/.fastq/.16BCgen.fastq}
cut -f 2 ${samp}_HaploTag_to_16BaseBCs | tail +2 > ${samp}_HaploTag_to_16BaseBCs.ema


## EMA Count
echo Running ema count
paste <(pigz -c -d ${read1/.fastq/.16BCgen.fastq} | paste - - - - | awk '{print $1"\t"$5"\t"$6"\t"$7}') <(pigz -c -d ${read2/.fastq/.16BCgen.fastq} | paste - - - - | awk '{print $1"\t"$5"\t"$6"\t"$7}' ) |\
tr "\t" "\n" |\
ema count -w ${samp}_HaploTag_to_16BaseBCs.ema -o $samp.16BCgen 2>$samp.16BCgen.log
cat $samp.16BCgen.log


## EMA Preprocessing
echo Running ema preproc
paste <(pigz -c -d ${read1/.fastq/.16BCgen.fastq} | paste - - - - | awk '{print $1"\t"$5"\t"$6"\t"$7}') <(pigz -c -d ${read2/.fastq/.16BCgen.fastq} | paste - - - - | awk '{print $1"\t"$5"\t"$6"\t"$7}' ) |\
tr "\t" "\n" |\
ema preproc -w ${samp}_HaploTag_to_16BaseBCs.ema -n 500 -t 40 -o ${samp}_ema-bin ${samp}.16BCgen.ema-ncnt 2>&1 |\
tee ${samp}_preproc.log



#### ALIGNMENT


## EMA Align
echo Running ema align 
for file in ${samp}_ema-bin/ema-bin-???; 
do echo "Processing $file"; ema align -t 4 -d -r $ref_genome -R $read_group -p 10x -s $file | samtools sort -@ 4 -O bam -l 0 -m 4G -o ${file}.sorted.bam;
done


## BWA align for non-BC reads
echo Running BWA for non-BC file
bwa mem -p -t 40 -M -R $read_group $ref_genome ${samp}_ema-bin/ema-nobc | samtools sort -@ 4 -O bam -l 0 -m 4G -o ${samp}_ema-bin/ema-nobc.sorted.bam


## Mark duplicates in BWA alignment 
echo Marking duplicates in BWA alignment
sambamba markdup -t 40 -p -l 0 ${samp}_ema-bin/ema-nobc.sorted.bam ${samp}_ema-bin/ema-nobc-pMarkedup.sorted.bam
rm ${samp}_ema-bin/ema-nobc.sorted.bam



#### POSTPROCESSING

## Merging all bam files
echo Merging all BAM files
sambamba merge -t 40 -p $fbname.bam ${samp}_ema-bin/*.bam


## Change barcodes from base-pairs to AxxBxxCxxDxx
samtools view -h $fbname.bam |\
awk 'BEGIN {split("AAAT,AAAG,AAAC,AATA,AATT,AATG,AATC,AAGA,AAGT,AAGG,AAGC,AACA,AACT,AACG,AACC,ATAA,ATAT,ATAG,ATAC,ATTA,ATTT,ATTG,ATTC,ATGA,ATGT,ATGG,ATGC,ATCA,ATCT,ATCG,ATCC,AGAA,AGAT,AGAG,AGAC,AGTA,AGTT,AGTG,AGTC,AGGA,AGGT,AGGG,AGGC,AGCA,AGCT,AGCG,AGCC,ACAA,ACAT,ACAG,ACAC,ACTA,ACTT,ACTG,ACTC,ACGA,ACGT,ACGG,ACGC,ACCA,ACCT,ACCG,ACCC,TAAA,TAAT,TAAG,TAAC,TATA,TATT,TATG,TATC,TAGA,TAGT,TAGG,TAGC,TACA,TACT,TACG,TACC,TTAA,TTAT,TTAG,TTAC,TTTA,TTTT,TTTG,TTTC,TTGA,TTGT,TTGG,TTGC,TTCA,TTCT,TTCG,TTCC,TGAA",val,","); \
for(i=1;i<=96;i++){lookup[val[i]]=sprintf("%02d",i)}};/BX:Z:/ {match($0,"BX:Z");bx=substr($0,RSTART,23);out="BX:Z:A"lookup[substr(bx,6,4)]"C"lookup[substr(bx,10,4)]"B"lookup[substr(bx,14,4)]"D"lookup[substr(bx,18,4)]substr(bx,22,2);gsub(bx,out,$0);print $0}; !/BX:Z/' |\
samtools view -@ 4 - -O BAM -o $fbname.BXnum.bam


## Index BAM file
samtools index -@ 4 $fbname.BXnum.bam


## Move final BAM and BAI file to output dir
mv $fbname.BXnum.bam* $ema_align/ && cd "$_"

echo DONE