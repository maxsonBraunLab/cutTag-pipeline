#!/bin/bash

# Merge replcate bigwigs output from callpeaks

usage(){
	echo "Usage: " 
        echo "  $0 -c chrom_sizes_file -o output_file rep1.bw rep2.bw rep3.bw"
	exit 1
}

file_exits(){
	local f="$1"
	[[ -f "$f" ]] && return 0 || return 1
}

cleanup() {
    rm tmp.bdg &> /dev/null
    rm tmp.sort.bdg &> /dev/null
}

[[ $# -eq 0 ]] && usage

script_args=()
while [ $OPTIND -le "$#" ]
do
    if getopts c:o: option
    then
        case $option
        in
            c) chrom_size="$OPTARG";;
            o) out_name="$OPTARG";;
        esac
    else
        script_args+=("${!OPTIND}")
        ((OPTIND++))
    fi
done


if [ -z "$chrom_size" ] 
then 
    echo "Please supply the chrom.sizes file." 
    usage
else 
    if ( ! file_exits "$chrom_size" )
    then 
        echo "error: file $chrom_size not found"
        usage 
    fi
fi


if [ -z "$out_name" ] 
then 
    echo "Please supply an output file name."
    usage
fi

echo "Using chrom.sizes file: $chrom_size" 
echo "Output file name will be: $out_name" 

# require bigWigMerge
if ! command -v bigWigMerge &> /dev/null
then
    echo "bigWigMerge not found"
    exit
fi

# require bigWigToBedGraph
if ! command -v bigWigToBedGraph &> /dev/null
then 
    echo "bigWigToBedGraph no found"
    exit
fi


echo "merging the following files..."
echo ${script_args[@]}

cmd="bigWigMerge ${script_args[@]} tmp.bdg"
echo $cmd
bigWigMerge ${script_args[@]} tmp.bdg

if [ $? -eq 0 ]; then
    echo "..."
else
    echo "Command bigWigMerge failed to complete..."
    echo "quitting.."
    cleanup
    exit 1
fi

# ensure sorted
sort -k1,1 -k2,2n tmp.bdg > tmp.sort.bdg

if [ $? -eq 0 ]; then
    echo "..."
else
    echo "Failed to sort merged bedgraph..."
    echo "quitting.."
    cleanup
    exit 1
fi

# convert to bigwig
bedGraphToBigWig tmp.sort.bdg $chrom_size $out_name

if [ $? -eq 0 ]; then
    echo "..."
else
    echo "Command bedGraphToBigWig failed..."
    echo "quitting.."
    cleanup
    exit 1
fi

#cleanup tmp files
cleanup
echo "done"
