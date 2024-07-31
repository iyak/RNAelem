#!/bin/bash
# usage: ./data_preparation.sh train.positive.fa train.negative.fa test.fa train.gz test.gz

# Check if exactly five arguments are provided
if [ "$#" -ne 5 ]; then
  echo "Five arguments required"
  exit 1
fi

min_len=99

{
    # Process Positive sequences
    while read -r header; do
        read -r seq
        if [[ ! $seq =~ N ]]; then
            seq=${seq//[^ACGT]/A}
            len=${#seq}
            trim_start=$(( (len - min_len) / 2 ))
            trim_end=$(( trim_start + min_len ))
            echo "$header; class:1"
            echo "${seq:$trim_start:$min_len}"
        fi
    done < "$1"

    # Process Negative sequences
    while read -r header; do
        read -r seq
        if [[ ! $seq =~ N ]]; then
            seq=${seq//[^ACGT]/A}
            len=${#seq}
            trim_start=$(( (len - min_len) / 2 ))
            trim_end=$(( trim_start + min_len ))
            echo "$header; class:0"
            echo "${seq:$trim_start:$min_len}"
        fi
    done < "$2"
} | gzip > "$4"
{
    # Process Test sequences
    while read -r header; do
        read -r seq
        if [[ ! $seq =~ N ]]; then
            seq=${seq//[^ACGT]/A}
            len=${#seq}
            trim_start=$(( (len - min_len) / 2 ))
            trim_end=$(( trim_start + min_len ))
            echo "$header; class:0"
            echo "${seq:$trim_start:$min_len}"
        fi
    done < "$3"
} | gzip > "$5"
