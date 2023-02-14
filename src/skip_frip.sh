#!/bin/bash

# For samples with peak calls (i.e. bed files with content), perform frip step as normal.
# For samples with no peak calls (i.e. empty bed files), skip frip and just create frip output file with 0 counts.

INPUT0=$1
INPUT1=$2
OUTPUT0=$3
OUTPUT1=$4

if [ -s ${INPUT0} ]
then
	plotEnrichment -b ${INPUT1} --BED ${INPUT0} --regionLabels 'frip' --outRawCounts ${OUTPUT1} -o ${OUTPUT0} > {log} 2>&1
else
	total_read_count=$(samtools view ${INPUT1} | wc -l)
	echo "file	featureType	percent	featureReadCount	totalReadCount" > ${OUTPUT1} #creates output file with column headers
	echo "${INPUT1}	frip	0	0	${total_read_count}" >> ${OUTPUT1} #appends to existing output file with name of sample and 0 counts
	touch ${OUTPUT0}
fi
