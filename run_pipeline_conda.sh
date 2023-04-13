#!/usr/bin/bash

#SBATCH --time 24:00:00
#SBATCH --partition exacloud
#SBATCH --job-name run_pipeline 
#SBATCH --output=jobs/run_pipeline_%j.log

# This is a wrapper script for running the cutTag pipeline as a batch job and using conda.

# Make sure to do the following before running this script:
# - conda activate snakemake environment
# - check that slurm profile for snakemake is set up
# - check that "jobs" folder exists

# To run this wrapper, do: sbatch run_pipeline_conda.sh

# set the number of jobs to run at a time (no spaces)
num_jobs=100

# run snakemake pipeline
snakemake -j $num_jobs \
--verbose \
--use-conda \
--profile slurm \
--cluster-config cluster.yaml \
--rerun-incomplete


exit