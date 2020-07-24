#!/bin/bash

perl -pe 's/  +/\t/g' | awk 'BEGIN{FS="\t"} (NR==2) {for (i = 1; i <= NF; i++) printf "--|";print "\n"$0} (NR==3){print ""$0} (NR>3 || NR <2) {print $0}' | perl -pe 's/\t/|/g'  | awk '(NR==2){print "|"$0} (NR!=2) {print "|"$0"|"}'
