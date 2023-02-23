#!/bin/bash

# For samples with peak calls (i.e. bed files with content), perform frip step as normal.
# For samples with no peak calls (i.e. empty bed files), skip frip and just create frip output file with 0 counts.

input_bed=$1
input_bam=$2
output_png=$3
output_tsv=$4
logfile=$5

if [ -s ${input_bed} ]
then
	plotEnrichment -b ${input_bam} --BED ${input_bed} --regionLabels 'frip' --outRawCounts ${output_tsv} -o ${output_png} > ${logfile} 2>&1
else
	total_read_count=$(samtools view ${input_bam} | wc -l)
	echo "file	featureType	percent	featureReadCount	totalReadCount" > ${output_tsv} #creates output file with column headers
	echo "${input_bam}	frip	0	0	${total_read_count}" >> ${output_tsv} #appends to existing output file with name of sample and 0 counts
	touch ${output_png}
fi
