#!/bin/bash

#  fragment length distribution
samtools view $1 | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' | sed 's/ /\t/g' | awk -v OFS='\t' '{print $2,$1}' > $2
