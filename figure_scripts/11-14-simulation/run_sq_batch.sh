#!/bin/bash

#SBATCH --account=????
#SBATCH --job-name=runSq
#SBATCH --time=1-00
#SBATCH --mail-user=?????
#SBATCH --mail-type=ALL
# Find total number of parameter lines (total number of jobs to be submitted across all arrays)

## Set channel set data #########################################################

# Use your username
uname=????

# S1 fitting
for ch_id_f in {1..156}; do 
	# Check number of jobs 
        # Job limit is 500, to leave n spare jobs for anything else, specify 500-n as the limit
        squeue -u $uname > job_list
        njobs=$(wc -l < job_list)
        echo "there are $njobs jobs"
        while [ $njobs -ge 495 ]; do
                echo "too many jobs, sleeping"
                sleep 30s
                squeue -u $uname > job_list
                njobs=$(wc -l < job_list)
                echo "slept, now there are $njobs jobs"
        done

        echo "array submitting (from line $line)" 
        channel_id=$ch_id_f area_id=1 sbatch ./sq_top10_batch.sh
        echo "submitted"
done

# S2 fitting
for ch_id_s in {1..101}; do
	# Check number of jobs 
        # Job limit is 500, to leave n spare jobs for anything else, specify 500-n as the limit
        squeue -u $uname > job_list
        njobs=$(wc -l < job_list)
        echo "there are $njobs jobs"
        while [ $njobs -ge 495 ]; do
                echo "too many jobs, sleeping"
                sleep 30s
                squeue -u $uname > job_list
                njobs=$(wc -l < job_list)
                echo "slept, now there are $njobs jobs"
        done

        echo "array submitting (from line $line)" 
        channel_id=$ch_id_s area_id=2 sbatch ./sq_top10_batch.sh 
        echo "submitted"
done

rm -rf job_list
