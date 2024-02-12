#!/usr/bin/bash

#SBATCH --time 24:00:00
#SBATCH --partition exacloud
#SBATCH --job-name run_pipeline 
#SBATCH --output=jobs/run_pipeline_%j.log

# This is a wrapper script for running the cutTag pipeline as a batch job and using Singularity.

# Make sure to do the following before running this script:
# - conda activate snakemake environment
# - check that slurm profile for snakemake is set up (use a slurm_singularity profile instead if regular slurm profile has Conda-specific configurations)
# - check that "jobs" folder exists in the main pipeline directory (if not, mkdir jobs)
# - add correct paths to indices and fastq folders below

# To run this wrapper, do: sbatch run_pipeline_singularity.sh

# set folder paths
# add absolute path to folder containing original fastq files (not the symlinks)
indices_folder="/home/groups/MaxsonLab/indices"
fastq_folder="/home/groups/MaxsonLab/nguythai/projects/pipeline_maintenance/cutTag-pipeline-singularity/.test/downsampled_fastqs"

# set the number of jobs to run at a time (no spaces)
num_jobs=100

# module load singularity before running snakemake
module load /etc/modulefiles/singularity/current

# run snakemake pipeline
# Note: if a separate Snakemake slurm profile for Singularity exists (e.g. slurm_singularity), you can use it instead of the default slurm profile
snakemake -j $num_jobs \
--verbose \
--use-singularity \
--singularity-args "--bind $indices_folder,$fastq_folder" \
--profile slurm_singularity \
--cluster-config cluster.yaml


exit
