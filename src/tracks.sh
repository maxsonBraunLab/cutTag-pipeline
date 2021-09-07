#!/bin/bash

while getopts "i:o:c:w:" op
do
    case "$op" in
        i)  input="$OPTARG";;
        o)  output="$OPTARG";;
        c)  chrom_size="$OPTARG";;
        w)  windows="$OPTARG";;
        \?) exit 1;;
    esac
done

# count total reads + get sample name
N=$(samtools view -c -@ 4 $input)
sample_name=$(basename $input | cut -d. -f1)

# genome coverage in CPM
bedtools genomecov -bga -ibam $input | awk -v OFS='\t' -v N=$N '{print $1,$2,$3,$4/N*1e6}' | sort -k1,1 -k2,2n > data/tracks/${sample_name}.bg

# calculate mean of signal by window
bedtools map -a $windows -b data/tracks/${sample_name}.bg -c 4 -o mean > data/tracks/${sample_name}.smooth.bg

# bedgraph to bigwig
bedGraphToBigWig data/tracks/${sample_name}.smooth.bg $chrom_size $output

rm data/tracks/${sample_name}.bg data/tracks/${sample_name}.smooth.bg