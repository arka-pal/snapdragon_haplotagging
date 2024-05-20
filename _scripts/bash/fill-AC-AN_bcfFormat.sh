#!/bin/bash

##### Bash script to add AC, AN annotations to a VCF file
##### author: Arka Pal
##### written: 22.01.2024
##### update: 22.01.2024

##### USAGE: bash fill-AC-AN.sh $1

##### $1 - VCF file (.bcf format)

## Change PATH
export PATH=~/.local/bin:$PATH

## Read input BCF file
inputVCF=$1

## Create new header file by adding AC and AN fields
echo "Creating new header file"
bcftools view -h $inputVCF > ${inputVCF/.bcf/.header}; 
awk '/^#/;
	/Dosage/ {
		print "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes for each ALT allele, in the same order as listed\">";
		print "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">";}' \
	${inputVCF/.bcf/.header}  |\
	awk '/FORMAT=<ID=GT,Number=1,Type=String,Description=/ { sub(" with posterior probability threshold of 0.9", "", $0) } 1' >\
	${inputVCF/.bcf/.header}.1; 
mv ${inputVCF/.bcf/.header}.1 ${inputVCF/.bcf/.header}


## Create new VCF file with INFO/AC and INFO/AN
echo "Adding INFO/AC and INFO/AN annotations"
cp ${inputVCF/.bcf/.header} ${inputVCF/.bcf/.tagged.vcf}; 
bcftools +fill-AN-AC $inputVCF | bcftools view -H - >> ${inputVCF/.bcf/.tagged.vcf}

## Compress and index
echo "Compress and index"
bgzip -f ${inputVCF/.bcf/.tagged.vcf}
bcftools tabix -f ${inputVCF/.bcf/.tagged.vcf.gz}