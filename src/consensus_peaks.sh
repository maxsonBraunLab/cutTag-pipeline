#!/bin/bash

# Searches for peaks with a given mark in their filename, and generates a consesnsus peak file for that mark
# a counts table is then generates from bamfiles with that mark 

# require bedtools
if ! command -v bedtools &> /dev/null
then
    echo "bedtools could not be found"
    exit
fi

# check args
if [[ $# -eq 0 ]]
then
    echo "no arguments supplied"
    exit 1
elif [[ $# -ne 4 ]] 
then
    echo "four positional arguments required"
    exit 1
fi

# the search pattern for the 
# epigenetic mark that occurs in the
# peak and bamfile names: e.g. 27Ac
MARK=$1

# directory of peak bedfiles
PEAKDIR=$2

# directory of bamfiles 
BAMDIR=$3

# output counts table filename
OUTFILE=$4

# array of replicates per mark
declare -a mpks=(${PEAKDIR}/*${MARK}*_peaks.bed)

# get union peak (consensus) across marks and conditions
cat ${mpks[@]} | sort -k1,1 -k2,2n \
    | bedtools merge -d -1 \
    | awk -v OFS='\t' '{print $1,$2,$3}' > data/counts/${MARK}_consensus.bed

# count the number of reads per union peak
bedtools multicov -bams ${BAMDIR}/*${MARK}*.bam -bed data/counts/${MARK}_consensus.bed -D > ${OUTFILE}_tmp

# label the counts table
ls ${BAMDIR}/*${MARK}*.bam  | sed 's!.*/!!' | cut -d_ -f1 | xargs | tr ' ' '\t' | awk '{print "chrom\tstart\tend\t" $0}' | cat - ${OUTFILE}_tmp > ${OUTFILE}

# remove tmp 
rm ${OUTFILE}_tmp
