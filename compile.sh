#!/bin/bash
if [[ ! -e $1 ]]
then
echo "Use: $0 target [flags] [out_name]"
exit
fi
target=$1
name=$(basename $1)".out"
outname="${name%%.*}"
if [[ $3 ]]
then
outname=$3
fi
gcc -std=c99 $2 -Wall -Doff64_t=__off64_t -lz $target -o $outname
