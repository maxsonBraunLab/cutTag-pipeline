#!/usr/bin/bash

while getopts "m:g:" op
do
	case "$op" in
		m)  mark="$OPTARG";;
		g)  genome="$OPTARG";;
		\?) exit 1;;
	esac
done

if [ -f $g ];
then
	echo "I found the FASTA file!"
else
	echo -e "$g does not exist. Exiting program."
	exit
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

# detect deseq2 output per mark ---------------------------------------------------------

# mark already defined
contrast=$(find "data/deseq2/${mark}" -name "*.bed" | cut -d/ -f4 | cut -d- -f1-2 | uniq)
up_05="data/deseq2/${mark}/${contrast}-${mark}-differential-up-05.bed"
up_01="data/deseq2/${mark}/${contrast}-${mark}-differential-up-01.bed"
down_05="data/deseq2/${mark}/${contrast}-${mark}-differential-down-05.bed"
down_01="data/deseq2/${mark}/${contrast}-${mark}-differential-down-01.bed"

declare -a differential_peaks=($up_05 $up_01 $down_05 $down_01)

# because there can be up to 4 homer runs per mark, we will run these 4 in the same node.

# HOMER analysis ------------------------------------------------------------------------

for file in ${differential_peaks[@]}; do

	echo "assessing $file for HOMER analysis"

	# make sure file exists
	if [ -f "$file" ]; then
		echo "$file detected"
	fi

	peak_count=$(wc -l "$file" | cut -d" " -f1)

	# run HOMER if DE peaks is greater than 10 peaks
	if [ "$peak_count" -lt 10 ]; then
		echo "ERROR: $file had less than 10 differential peaks"

	else

		# extract dir,sig info for the run.
		echo "Running HOMER for $peak_count peaks in $file"
		directionality=$(basename $file | cut -d- -f5) # "up" or "down"
		sig=$(basename $file | cut -d- -f6 | cut -d. -f1) # e.g. 0.01 or 0.05

		# file I/O
		output_dir="data/homer/${contrast}-${mark}-${directionality}-${sig}"
		log="data/logs/homer_${contrast}-${mark}-${directionality}-${sig}.log"
		log_file_exists $log

		# HOMER command in {{{ parallel }}}
		cmd="findMotifsGenome.pl $file $genome $output_dir -size 200 > $log 2>&1 &"
		eval $cmd

	fi

done
wait

touch data/homer/$mark.done
