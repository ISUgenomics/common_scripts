#!/bin/bash

cat $1 | perl -pe 's/ /_/g' | paste - - - - | awk '{print ">"$1;print $2}' > ${1%%.*}.fasta
