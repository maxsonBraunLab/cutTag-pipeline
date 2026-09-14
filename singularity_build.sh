#!/usr/bin/bash

#SBATCH --time 8:00:00
#SBATCH --mem=8G
#SBATCH --partition batch
#SBATCH --job-name singularity_build
#SBATCH --output=jobs/singularity_build/singularity_build_%j.log

# This is a script for building custom Apptainer/Singularity containers. By default, this script builds container images from the Apptainer/Singularity definition files (.def) found in a specific folder in the main pipeline directory.
# This script accepts one command line argument specifying the path to a folder in which to store the Apptainer/Singularity image files.
# Run this script from the main pipeline directory (where the Snakefile is)


# command line inputs
output_image_folder=$1

echo -e "|--- Output image folder:\n${output_image_folder}\n"

# if output folder doesn't exist, create it
# but if mkdir fails due to empty/invalid input, then echo error message and exit script
if [ ! -d "$output_image_folder" ]
then
	mkdir -p $output_image_folder || { echo "Error: Invalid output folder specified."; exit 1; }
fi

# build container image files
for def_file in $(find singularity_definition_files -name "*.def")
do
	image_filepath="${output_image_folder}/$(basename -s .def ${def_file}).sif"
	
	echo -e "=============================================\n"
	echo -e "[ $(date) ]\n"
	echo -e "|--- Definition file:\n${def_file}\n"
	echo -e "|--- Output image file:\n${image_filepath}\n"
	
	apptainer build ${image_filepath} ${def_file}
	
	echo -e "[ $(date) ]\n"
done

exit

