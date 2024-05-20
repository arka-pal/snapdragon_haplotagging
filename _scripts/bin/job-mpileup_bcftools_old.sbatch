#!/bin/bash

##### Usage: sbatch job-mpileup.slurm.sh </path/to/reference_genome> </path/to/input_bams> <output.vcf.gz> <region>
##### $1: reference_genome, $2: input_bams, $3: output.vcf.gz, $4: region


#SBATCH --job-name=mpileup
#SBATCH --output=job-mpileup.slurm-%A.out
#SBATCH --error=job-mpileup.slurm-%A.err

#SBATCH --partition=gpu
#SBATCH --time=96:00:00
#SBATCH -c 1
#SBATCH --mem-per-cpu=50G
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=ALL

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

## load any module you need here
module load bcftools

## run commands on SLURM's srun


#srun bcftools mpileup -a AD,DP,QS,SP -Ob -f ~/snap_hap/ref_genome/v3.5/Amajus_v3.5.fa -b ~/snap_hap/dirtyVCF_for_LouiseA/bamlist.txt | bcftools call -a GQ,GP -m -Oz -o dirtyVCF.vcf.gz
#or you could directly put all input bam files!



srun bcftools mpileup -a AD,DP,QS,SP -r $4 -Ob -f $1 -b $2 | bcftools call -a GQ,GP -m -v -V indels -Oz -o $3

srun bcftools tabix $3

#srun bcftools stats $2 > ${2/.gz/.stats}

#srun bcftools view --max-alleles 2 --exclude-types indels $2 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > ${2/.vcf.gz/.filtered.vcf}

srun bcftools view -m2 -M2 -v snps -Oz -o ${3/.vcf.gz/__biSNPs.vcf.gz} $3

srun bcftools tabix ${3/.vcf.gz/__biSNPs.vcf.gz}