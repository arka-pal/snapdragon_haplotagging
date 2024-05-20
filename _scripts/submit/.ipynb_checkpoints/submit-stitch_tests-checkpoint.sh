##### Bash script to submit job-impute_stitch jobs to cluster
##### author: Arka Pal
##### written: 10.02.2023
##### update: 12.07.2023

##### USAGE: bash submit-stitch_tests.sh <options>
##### $1 - chr
##### $2 - K
##### $3 - downsampletoCov (default is 50)
##### $4 - use_bx_tag (TRUE or FALSE)
##### $5 - niterations
##### $6 - ngen
##### $7 - recombination rate


WORKDIR='/nfs/scistore18/bartogrp/apal/snap_hap/impute/stitch_tests'

for posfile in $WORKDIR/$1/stitch_regions/posfile_v1/*.pos;
    do echo posfile: $posfile
    chr=$1
    K=$2
    downsampleToCov=$3
    bx_tag=$4
    niter=$5
    ngen=$6
    expRate=$7
    snpID=$(cut -d- -f1 <(echo ${posfile##*/}))
    sbatch -J run2_$chr-${snpID}-K$K-cov$downsampleToCov-bxTRUE-niter${niter}-ngen${ngen}-r${expRate}-plotFALSE \
            ~/snap_hap/_scripts/sbatch/impute/job-stitch_tests.sbatch \
            $chr $K $posfile $downsampleToCov TRUE $niter $ngen $expRate
done
