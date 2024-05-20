#!/bin/bash
##### Script to add FORMAT/PG and FORMAT/PL fields to STITCH-generated output using awk.
##### written by: Frank Chan, Marek Kucka & Arka Pal
##### modified by: Arka Pal, last Update 2023-09-26

##### Usage: ./add_PG_PL.sh <options>
##### $1: input vcf


i=$1

## Index vcf file
tabix -f $i;

## Create new VCF header 
bcftools view -h $i > ${i/.vcf.gz/.header}; 
awk '/^#/;
	/reference and alternate alleles/ {print "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"\">"};
	/Dosage/ {
		print "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">";
		print "##FORMAT=<ID=PG,Number=1,Type=String,Description=\"Best Guessed Genotype with posterior probability threshold of 0.9\">";}' \
	${i/.vcf.gz/.header}  |\
	awk '/FORMAT=<ID=GT,Number=1,Type=String,Description=/ { sub(" with posterior probability threshold of 0.9", "", $0) } 1' >\
	${i/.vcf.gz/.header}.1; 
mv ${i/.vcf.gz/.header}.1 ${i/.vcf.gz/.header}



## Create new VCF file with best guess genotypes
cp ${i/.vcf.gz/.header} ${i/.vcf.gz/.PL.vcf}; 
bcftools view $i -H | \
	awk 'BEGIN {OFS="\t"};
		{$9="GT:PG:GP:DS:PL"; 
			for (i=10; i<=NF; i++) {
				split($i,field,":");
				pg=field[1];
				split(field[2],gp,",");
				if(gp[1]>gp[2] && gp[1]>gp[3]) {field[1]="0/0"};
				if(gp[2]>gp[1] && gp[2]>gp[3]) {field[1]="0/1"};
				if(gp[3]>gp[1] && gp[3]>gp[2]) {field[1]="1/1"};
				pl1=-log(gp[1])/log(10)*10;
				pl2=-log(gp[2])/log(10)*10;
				pl3=-log(gp[3])/log(10)*10;
				if(pl1 > 255) pl1=255;
				if(pl2 > 255) pl2=255;
				if(pl3 > 255) pl3=255;
				plMIN=pl1; 
				if (pl2 < plMIN) plMIN=pl2; 
				if (pl3 < plMIN) plMIN=pl3;
				pl1_weighted=pl1-plMIN;
				pl2_weighted=pl2-plMIN;
				pl3_weighted=pl3-plMIN;
				$i=field[1]":"pg":"field[2]":"field[3]":"int(pl1_weighted+0.5)","int(pl2_weighted+0.5)","int(pl3_weighted+0.5)
			};		
		print $0}' >> \
	${i/.vcf.gz/.PL.vcf}


## Compress and index
bgzip -f ${i/.vcf.gz/.PL.vcf}; 
tabix -f ${i/.vcf.gz/.PL.vcf.gz}