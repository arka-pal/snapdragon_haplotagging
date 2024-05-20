#!/bin/bash

#SBATCH --job-name=download_sequence
#SBATCH --output=download_sequence.out
#SBATCH --error=download_sequence.err

#SBATCH --partition=defaultp
#SBATCH --time=120:00:00
#SBATCH -c 1
#SBATCH --mem-per-cpu=1G
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1

#SBATCH --mail-user=arka.pal@ist.ac.at
#SBATCH --mail-type=ALL

#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

## Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

##################

### # MEGA samples
#srun wget -r -np -nH --cut-dirs=4 -I "fml/ag-chan/grandpaSon/10_Plates_2022/2xCoverage/*e/" 'ftp://ftp.tuebingen.mpg.de/fml/ag-chan/grandpaSon/10_Plates_2022/2xCoverage/'

### 10xNEW samples
#srun wget -r -np -nH --cut-dirs=4 ftp://ftp.tuebingen.mpg.de/fml/ag-chan/grandpaSon/Plate9_Plate10_Haplo4/Plate9/
#srun wget -r -np -nH --cut-dirs=4 ftp://ftp.tuebingen.mpg.de/fml/ag-chan/grandpaSon/Plate9_Plate10_Haplo4/Plate10/
srun wget -r -np -nH --cut-dirs=4 ftp://ftp.tuebingen.mpg.de/fml/ag-chan/grandpaSon/Plate9_Plate10_Haplo4/Plate9/Plate9_C[6-9]*