#!/usr/bin/bash

#SBATCH --time 8:00:00
#SBATCH --mem=8G
#SBATCH --partition exacloud
#SBATCH --job-name singularity_build_remote
#SBATCH --output=jobs/singularity_build_remote/singularity_build_remote_%j.log

# This is a script for building Singularity containers via the sylabs.io remote builder. By default, this script builds container images from the Singularity  definition files (.def) found in a specific folder in the main pipeline directory.
# This script accepts one command line argument specifying the path to a folder in which to store the Singularity image files.
# Run this script from the main pipeline directory (where the Snakefile is)

# Make sure to do the following before running this script:
# - check that "jobs/singularity_build_remote" folder exists in the main pipeline directory (if not, mkdir -p jobs/singularity_build_remote)
# 
# - If you don't already have singularity remote login configured for your terminal, then:
# 		- use srun to get onto an interactive/compute node
# 		- run `module load /etc/modulefiles/singularity/current` from anywhere on a compute node
# 		- get a sylabs access token from sylabs.io (use Github to login to sylabs.io) 
# 		- run `singularity remote login` and input your sylabs token


# module load singularity before running build
module load /etc/modulefiles/singularity/current

# command line inputs
# make sure folder path doesn't have "/" at end 
output_image_folder=$1

echo -e "|--- Output image folder:\n${output_image_folder}\n"

# if output folder doesn't exist, create it
# but if mkdir fails due to empty/invalid input, then echo error message and exit script
if [ ! -d "$output_image_folder" ]
then
	mkdir -p $output_image_folder || { echo "Error: Invalid output folder specified."; exit 1; }
fi

# build container image files
for def_file in $(ls singularity_definition_files/*.def)
do
	image_filepath="${output_image_folder}/$(basename -s .def ${def_file}).sif"
	
	echo -e "=============================================\n"
	echo -e "[ $(date) ]\n"
	echo -e "|--- Definition file:\n${def_file}\n"
	echo -e "|--- Output image file:\n${image_filepath}\n"
	
	singularity build --remote ${image_filepath} ${def_file}

	echo
done

echo -e "[ $(date) ]\n"
echo "Done!"

exit
