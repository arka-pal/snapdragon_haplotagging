#!/bin/bash

##### SLURM script for bcftools pileup, calling and filtering variants for n96 samples
##### author: Arka Pal
##### written: 26.10.2023
##### update: 30.10.2023

##### Usage: sbatch -J <job-name> job-bcftools_mpileup-n96.sbatch $1 $2 $3 $4 $5 $6
##### $1: chromosome; Chr1
##### $2: reference genome; $WORKDIR/ref_genome/v3.5/Amajus_v3.5.fa
##### $3: bamlist; bamlist=$WORKDIR/sample_info/bam_info/bams_n96.txt
##### $4: outVCF; Am_all_bcftools_Chr1.vcf.gz


## Define SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------
## Defining SLURM variables


## Load modules
module load bcftools/1.18


## Read inputs
chrom=$1
ref=$2
bamlist=$3
outVCF=$4


## Define variables


## Print variables
echo chrom: $chrom
echo Reference Genome: $ref
echo BAM Files: $bamlist
echo outVCF: $outVCF
echo -e "\n"


# ## Run bcftools mpileup and call variants
# echo -e Calling SNPs
# srun time bcftools mpileup --annotate AD,ADF,ADR,DP,QS,SP \
#         -r $chrom \
#         -f $ref \
#         -Ou -d 500 \
#         -b $bamlist \
#         --threads ${SLURM_CPUS_PER_TASK} | \
#     bcftools call --annotate GQ,GP \
#         -m -Oz \
#         -o $outVCF \
#         --threads ${SLURM_CPUS_PER_TASK}
        

# ## Index VCF
# echo -e Indexing VCF file
# srun time bcftools tabix -f $outVCF


## Filter VCF
echo -e Filtering VCF file
srun time bcftools filter --threads ${SLURM_CPUS_PER_TASK} -e "INFO/DP>165 | QUAL<20 | MQ<30 | F_MISSING>0.50 | AC==0 | AC==AN" --SnpGap 5 $outVCF | \
    bcftools view --threads ${SLURM_CPUS_PER_TASK} -m2 -M2 -v snps -Oz -o ${outVCF/.vcf.gz/_biSNPs_filtered.vcf.gz} -
bcftools tabix -f ${outVCF/.vcf.gz/_biSNPs_filtered.vcf.gz}