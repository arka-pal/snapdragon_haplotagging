file=$1;
echo $file
fbname=$(basename $file .bam)
echo $fbname
dir=$(dirname $file)
echo $dir
#marks duplicates
java -Xmx8g -XX:ParallelGCThreads=10 -jar /fml/chones/local/picard-2.18.25/picard.jar MarkDuplicates \
       I=$file \
       TMP_DIR=./ \
       O=$dir/$fbname.pMark.bam \
       M=$dir/$fbname.pMark.metrics \
CREATE_INDEX=TRUE READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX VALIDATION_STRINGENCY=LENIENT \
#de-duplicates the bam
samtools view -@ 10 -h $fbname.pMark.bam -F 1024 -O bam -o $fbname.dedup.bam
#BX-tag sorts the deduplicated bam
samtools sort -@ 10 -t BX $fbname.dedup.bam -T ./$fbname.tmpsort -O BAM -o $fbname.dedup.bxsorted.bam
#make linked-read molecules off of reads with same BX-tag within 50kbp from each other
bed_write.full.pl $fbname.dedup.bxsorted.bam
#removes 00 containing molecules
awk '!/A00C|C00B|B00D|D00/' $fbname.dedup.bxsorted.linked_reads.full.bed > $fbname.dedup.bxsorted.linked_reads.full.no00.bed
#makes histogram of reads per molecule, e.g. how many DNA molecules have 50 reads or 10 reads etc
cut -f10 $fbname.dedup.bxsorted.linked_reads.full.no00.bed | datamash -s groupby 1 count 1 | sort -k1nr > $fbname.reads.per.mol.log
#writes histogram of molecule size in 1kb bins, sorted from longest to shortest molecule
awk '{ print $3-$2 }' $fbname.dedup.bxsorted.linked_reads.full.no00.bed | datamash bin:1000 1 | datamash -s groupby 1 count 1 | sort -k1nr > $fbname.1kb-bin.molecule.histogram
reads=$(cut -f10 $fbname.dedup.bxsorted.linked_reads.full.no00.bed | datamash sum 1 | awk '{print $1/2}')
size=$(awk '{ print $3-$2 }' $fbname.dedup.bxsorted.linked_reads.full.no00.bed | datamash sum 1 | awk '{print $1/2}')
#calculates N50_reads_per_molecule
cut -f10 $fbname.dedup.bxsorted.linked_reads.full.no00.bed | sort -k1nr | awk -v val="$reads" '($1+prev)>val{exit} ($1+prev)<=val; {prev+=$1}' |tail -1 > $fbname.N50_reads_per_mol.info
#calculates N50_molecule_length
awk '{ print $3-$2 }' $fbname.dedup.bxsorted.linked_reads.full.no00.bed | sort -k1nr |  awk -v val="$size" '($1+prev)>val{exit} ($1+prev)<=val; {prev+=$1}' | tail -1 > $fbname.N50_molecule_size.info
#writes all size-sorted DNA molecules; showing Mol-size, BX-BC, number of reads/mol, and chr location
awk '{ print $3-$2"\t"$4"\t"$10"\t"$1"\t"$2"\t"$3 }' $fbname.dedup.bxsorted.linked_reads.full.no00.bed | sort -k1nr > $fbname_Mol-size_BC_ReadsPerMol_location.log
rm $fbname.dedup.bam
rm $fbname.dedup.bxsorted.bam