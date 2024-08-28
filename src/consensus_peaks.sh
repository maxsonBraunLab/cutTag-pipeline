#!/bin/bash

# Searches for peaks with a given mark in their filename.
# Loops over every condition and finds peaks that appear in >= n replicates. Export to tmp files.
# Merge all the tmp files and export as consensus peak file per mark.
# then get a counts table for that mark with respective BAM files.

# require bedtools
if ! command -v bedtools &> /dev/null
then
    echo "ERROR: Bedtools could not be found."
    exit 1
fi

# check args
if [[ $# -eq 0 ]]
then
    echo "ERROR: No command line arguments supplied."
    exit 1
fi

# set default variable values
USE_BLACKLIST_BED=false

# read in command line args
while getopts "m:n:o:b" op
do
	case "$op" in
		b)  USE_BLACKLIST_BED=true;;
		m)  MARK="$OPTARG";;
		n)  N_INTERSECTS="$OPTARG";;
		o)  OUTFILE="$OPTARG";;
		\?) exit 1;;
	esac
done


# array of replicates per mark
declare -a mpks=(data/callpeaks/*${MARK}_peaks.bed)

all_samples=$(echo ${mpks[@]} | tr ' ' '\n' | cut -d/ -f3 | cut -d_ -f1-3)
all_conditions=$(echo ${mpks[@]} | tr ' ' '\n' | cut -d/ -f3 | cut -d_ -f1 | sort | uniq)

echo -e "All samples:\n${all_samples}\n"
echo -e "All conditions:\n${all_conditions}\n"

for condition in $all_conditions; do

	echo -e "|-- Condition: ${condition}\n"

    # file I/O
    tmp_output="data/counts/$MARK.$condition.tmp.bed"

    # list all replicates in one condition
	if [ "$USE_BLACKLIST_BED" = true ]
	then
		echo "Blacklist flag enabled."
		all_replicates=$(find data/callpeaks/ -name "$condition*${MARK}_peaks_noBlacklist.bed" | sort | tr '\n' ' ')
	else
		echo "Blacklist flag not enabled."
		all_replicates=$(find data/callpeaks/ -name "$condition*${MARK}_peaks.bed" | sort | tr '\n' ' ')
	fi

	echo -e "All replicates:\n${all_replicates}\n"

    # find widest peak that appear in at least n replicates + export to tmp file.
    cat $all_replicates | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge | \
    bedtools intersect -c -a stdin -b $all_replicates | \
    awk -v n=$N_INTERSECTS '$4 >= n' | awk -v OFS='\t' '{print $1,$2,$3}' > $tmp_output

done

# merge intervals for all tmp files and export as consensus peak.
all_temp_files=$(find data/counts -name "${MARK}.*.tmp.bed" | sort | tr '\n' ' ')

echo -e "All temp files:\n${all_temp_files}\n"

cat $all_temp_files | sort -k1,1 -k2,2n | bedtools merge > data/counts/${MARK}_consensus.bed
rm $all_temp_files

# count the number of reads per union peak
bedtools multicov -bams data/markd/*${MARK}.sorted.markd.bam -bed data/counts/${MARK}_consensus.bed -D > ${OUTFILE}_tmp

# label the counts table
ls data/markd/*${MARK}.sorted.markd.bam  | sed 's!.*/!!' | cut -d. -f1 | xargs |  tr ' ' '\t' | awk '{print "chrom\tstart\tend\t" $0}' | cat - ${OUTFILE}_tmp > ${OUTFILE}

# remove tmp 
rm ${OUTFILE}_tmp
