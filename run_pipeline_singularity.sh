#!/usr/bin/bash

#SBATCH --time 24:00:00
#SBATCH --partition exacloud
#SBATCH --job-name run_pipeline 
#SBATCH --output=jobs/run_pipeline_%j.log

# This is a wrapper script for running the cutTag pipeline as a batch job and using Singularity.

# Make sure to do the following before running this script:
# - be on an interactive/compute node so that Singularity module can be activated
# - conda activate snakemake environment
# - check that slurm_singularity profile for snakemake is set up
# - check that "jobs" folder exists in the main pipeline directory (if not, mkdir jobs)
# - add correct paths to indices and fastq folders below

# To run this wrapper, do: sbatch run_pipeline_singularity.sh

# set folder paths
indices_folder="/home/groups/MaxsonLab/indices"
conda_folder="${CONDA_PREFIX_1}/envs"
conda_pkgs_folder="${CONDA_PREFIX_1}/pkgs"
fastq_folder="" # add absolute path to folder containing original fastq files (not the symlinks)

# set the number of jobs to run at a time (no spaces)
num_jobs=100

# module load singularity before running snakemake
module load /etc/modulefiles/singularity/current

# run snakemake pipeline
snakemake -j $num_jobs \
--verbose \
--use-conda \
--use-singularity \
--singularity-args "--bind $indices_folder,$conda_folder,$conda_pkgs_folder,$fastq_folder" \
--profile slurm_singularity \
--conda-prefix $conda_folder \
--cluster-config cluster.yaml


exit
