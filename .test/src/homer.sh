#!/usr/bin/bash

while getopts "m:g:s:p:" op
do
	case "$op" in
		m)  mark="$OPTARG";;
		g)  genome="$OPTARG";;
		s)  slurm="$OPTARG";;
		p)  p="$OPTARG";;
		\?) exit 1;;
	esac
done

if [ -f $genome ];
then
	echo "I found the FASTA file!"
else
	echo -e "$g does not exist. Exiting program."
	exit
fi

if [[ $slurm == 1 || $slurm == 0 ]]; then
	true
else
	echo "ERROR: -s must be either 0 or 1 for SLURM submission."
	exit 1
fi

log_file_exists() {

	if [ -f $1 ]; then
		rm $1
		touch $1
	fi

	if [ ! -f $1 ]; then
		touch $1
	fi

}

# detect contrast combinations per mark
contrasts=$(find "data/deseq2/${mark}" -name "*.bed" | cut -d/ -f4 | cut -d- -f1-2 | sort | uniq)

# HOMER analysis per contrast -----------------------------------------------------------

for contrast in $contrasts; do

	echo "assessing contrast $contrast in mark $mark"

	up_05="data/deseq2/${mark}/${contrast}-${mark}-differential-up-05.bed"
	up_01="data/deseq2/${mark}/${contrast}-${mark}-differential-up-01.bed"
	down_05="data/deseq2/${mark}/${contrast}-${mark}-differential-down-05.bed"
	down_01="data/deseq2/${mark}/${contrast}-${mark}-differential-down-01.bed"

	declare -a differential_peaks=($up_05 $up_01 $down_05 $down_01)

	# assess and run HOMER per bed file
	for file in ${differential_peaks[@]}; do

		# check if file exists
		if [ -f "$file" ]; then

			echo "$file detected"
			peak_count=$(cat $file | wc -l)

			# run HOMER if >= 10 peaks
			if [ "$peak_count" -gt 10 ]; then

				# extract dir,sig info for the run.
				echo "Running HOMER for $peak_count peaks in $file"
				directionality=$(basename $file | cut -d- -f5) # "up" or "down"
				sig=$(basename $file | cut -d- -f6 | cut -d. -f1) # e.g. 0.01 or 0.05

				# file I/O
				output_dir="data/homer/${contrast}-${mark}-${directionality}-${sig}"
				log="data/logs/homer_${contrast}-${mark}-${directionality}-${sig}.log"
				log_file_exists $log

				# HOMER local run/submission to SLURM
				if [ $slurm == 0 ]; then
					findMotifsGenome.pl $file $genome $output_dir -size 200 -p $p > $log 2>&1 &
				fi

				if [ $slurm == 1 ]; then
					job_out="jobs/homer-$contrast.out"
					job_err="jobs/homer-$contrast.err"
					sbatch -e $job_err -o $job_out --job-name "HOMER" --time "02:00:00" --mem="8G" --cpus-per-task=$p --wait --wrap="findMotifsGenome.pl $file $genome $output_dir -size 200 -p $p > $log 2>&1" &
				fi
			fi
		fi
	done
done

wait
touch data/homer/$mark.done
