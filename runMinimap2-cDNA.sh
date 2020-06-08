#!/bin/bash
genome=$1
cdna=$2
out=$(basename ${genome%.*} |cut -f 1 -d ".")
minimap2 -ax splice -a -uf -C1 -k 12 -t 36 ${genome} --cs ${cdna} > ${out}-${cdna%.*}.sam
