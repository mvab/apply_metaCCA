#!/bin/bash

# input from command line

chr="$1"  # chromosome to process; add ' -F "chr#" ' at the end to submission script
cores=16
data_path="/"


# request resources:
#PBS -l nodes=1:ppn=16
#PBS -l walltime=24:00:00

## not using: #PBS -q himem


# Define working directory
export WORK_DIR=$HOME/metaCCA_big_jobs

# load latest R
module add languages/R-3.5.1-ATLAS-gcc-6.1

# on compute node, change directory to 'submission directory':
cd $WORK_DIR


echo JOB ID: $PBS_JOBID
echo Working Directory: `pwd`
echo Start Time: `date`

# run your program, timing it:

# gene-based analysis
time Rscript 3_run_metaCCA_analysis.R  ${cores} genes $chr ${data_path}
# snp-based analysis
#time Rscript 3_run_metaCCA_analysis.R  ${cores} snps all ${data_path}

echo End Time: `date`


