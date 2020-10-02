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

declare -a mpks=(${PEAKDIR}/*${MARK}*_peaks.bed)

# intersect replicate peaks across marks
bedtools intersect -a ${mpks[0]} -b ${mpks[@]:1} -c | awk -v OFS='\t' '$7>1 {print $1,$2,$3,$4,"0","."}' > data/counts/${MARK}_consensus.bed

# count the number of reads permark 
bedtools multicov -bams ${BAMDIR}/*${MARK}*.bam -bed data/counts/${MARK}_consensus.bed -D > ${OUTFILE}_tmp

# label the counts table
ls ${BAMDIR}/*${MARK}*.bam  | cut -d_ -f1 | xargs | tr ' ' '\t' | awk '{print "chrom\tstart\tend\tpeak\tscore\tstrand\t" $0}' | cat - ${OUTFILE}_tmp > ${OUTFILE}

# remove tmp 
rm ${OUTFILE}_tmp
