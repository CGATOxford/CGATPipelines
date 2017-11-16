#!/bin/bash
#
# Apply a method (md5sum, wc) on a file.
# 
# If file ends with .gz, uncompress and compute
# from stream.
#
# In both cases, ignore any lines starting with "#"
#
suffix="${1##*.}"

fn=$1
printf "%s\t" "$fn"
shift

if [[ $suffix = 'gz' ]]; then
    zcat $fn | grep -v '^#' | "$@"
elif [[ $suffix = 'bam' ]]; then
    # do not include header, as file path might differ if temp
    # files are used
    samtools view  $fn | "$@"
else
    cat "$fn" | grep -v '^#' | "$@"
fi
