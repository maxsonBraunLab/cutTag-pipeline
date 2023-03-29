#!/usr/bin/bash

#SBATCH --time 24:00:00
#SBATCH --partition exacloud
#SBATCH --job-name run_pipeline 
#SBATCH --output=jobs/run_pipeline_%j.log

# This is a wrapper script for running the cutTag pipeline as a batch job and using conda.

# conda activate snakemake environment before running script 
# make sure slurm profile for snakemake is set up before running
# make sure "jobs" folder exists before running

# To run this wrapper, do: sbatch run_pipeline_conda.sh

# set the number of jobs to run at a time (no spaces)
num_jobs=100

snakemake -j $num_jobs \
--verbose \
--use-conda \
--profile slurm \
--cluster-config cluster.yaml \
--rerun-incomplete


exit
