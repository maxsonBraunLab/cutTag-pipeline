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
    rm ${out_name}.tmp.bdg &> /dev/null
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
if ! command -v wiggletools &> /dev/null
then
    echo "wiggletools not found"
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

cmd="wiggletools mean ${script_args[@]} | grep -v "_" | awk 'NF==4{print}' > ${out_name}.tmp.bdg"
echo $cmd
wiggletools mean ${script_args[@]} | grep -v "_" | awk 'NF==4{print}' > ${out_name}.tmp.bdg

if [ $? -eq 0 ]; then
    echo "..."
else
    echo "Command wiggletools failed to complete..."
    echo "quitting.."
    cleanup
    exit 1
fi

# convert to bigwig
bedGraphToBigWig ${out_name}.tmp.bdg $chrom_size $out_name

if [ $? -eq 0 ]; then
    echo "..."
else
    echo "Command bedGraphToBigWig failed..."
    echo "quitting.."
    cleanup
    exit 1
fi

cleanup
echo "done"

