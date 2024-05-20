##### Bash script to submit EMA-align_only jobs to cluster
##### author: Arka Pal
##### last update: 12.12.2022

##### USAGE: bash file.sh $ref_genome $ref_gen_ver
##### USAGE Example: bash ~/snap_hap/scripts/submit/run-EMAalign.sh ~/snap_hap/ref_genome/v3.5/Amajus_v3.5.fa v3.5

ref_genome=$1
ref_gen_ver=$2

WORKDIR='/nfs/scistore18/bartogrp/apal/snap_hap'

## Set working directory
cd $WORKDIR/bams/$ref_gen_ver/jobs/
#batch_list=(*)
#batch_list=(2x-N701 2x-N702 2x-N703 2x-N704 2x-N705 2x-N706 2x-N707 2x-N708 2x-N709 2x-N710)
#batch_list=(10x)
#batch_list=(60x 60xE)
#batch_list=(n96_Ave n96_Pla n96_misc TRIO 2xE-N705e 2xE-N706e 2xE-N707e 2xE-N708e 2xE-N709e 2xE-N710e)

# echo ${#batch_list[@]} ##no of batches


for batch in ${batch_list[@]};
do 
    ema_path=$WORKDIR/bams/v3.5.SL/ema_dir/$batch;
    # echo $ema_path;
    no_of_samples=$(ls $ema_path | wc -l);
    
    if [ ! -d $WORKDIR/bams/$ref_gen_ver/jobs/$batch ]; 
    then 
        mkdir -p $WORKDIR/bams/$ref_gen_ver/jobs/$batch; 
    fi
    
    cd $WORKDIR/bams/$ref_gen_ver/jobs/$batch
    # echo $(pwd)		
    echo "Submitting $no_of_samples jobs for $batch";
    sbatch -a 1-$no_of_samples \
            -J $batch-$ref_gen_ver \
            $WORKDIR/scripts/sbatch/job-EMAalign_only.sbatch $batch $ref_genome $ref_gen_ver;
done