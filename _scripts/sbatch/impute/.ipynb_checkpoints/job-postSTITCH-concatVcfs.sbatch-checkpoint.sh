#!/bin/bash

##### SLURM script to concat STITCH vcfs
##### author: Arka Pal
##### written: 02.10.2023
##### update: 11.10.2023

##### USAGE: sbatch -J <job-name> job-postSTITCH-concatVcfs.sbatch <options>
##### $1 - chrom
##### $2 - stitchVcfList
##### $3 - outVCF; e.g., Am_all_stitch_${chrom}_SnpGap5_biSNPs_filtered-DP_500-7732_QUAL20_MQ30.vcf.gz


## Defining SLURM variables
#----------------------------------------------------------------
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.out
#SBATCH --open-mode=append

#SBATCH --partition=defaultp
### #SBATCH --partition=debian12
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4G

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=FAIL,END

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#----------------------------------------------------------------


## Load modules
module load bcftools/1.18


## Read inputs
chrom=$1
stitchVcfList=$2
outVCF=$3

echo -e "\n"
echo chrom: $chrom
echo stitchVcfList: $stitchVcfList
echo outVCF: $outVCF
echo -e '\n'


## Merge VCFs
echo -e "\nMerging VCFs\n"
srun time bcftools concat -a --threads ${SLURM_CPUS_PER_TASK} -Oz -o $outVCF -f $stitchVcfList
echo -e "\nIndexing VCFs\n"
srun time bcftools index -f --tbi --threads ${SLURM_CPUS_PER_TASK} $outVCF
echo -e "\nbcftools stats\n"
srun time bcftools stats --threads ${SLURM_CPUS_PER_TASK} -s- --af-bins 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1 --depth 0,25000,25 $outVCF > ${outVCF/.vcf.gz/.stats}