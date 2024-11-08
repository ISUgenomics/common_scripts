#!/bin/bash

# runGenomeScope R1Read R2read 21

readR1=$1
#readR2=$2
kmer=$2
outName=$(basename $1 .gz)
outName=${outName}_${kmer}

module load jellyfish/2.2.7-py2-euxjh7e
jellyfish count -C -m ${kmer} -s 3000000000 -t 30 <(zcat ${readR1} ${readR2}) -o $outName.jf
# Create histogram from that data

jellyfish histo -t 36 ${outName}.jf > ${outName}.hist

# Run genomescope with hist or
# once the histo file is created, visit http://qb.cshl.edu/genomescope/ website to upload the histo file

module use /work/gif/modules/software/
module load genomescope2.0
module load r

genomescope.R -i ${outName}.hist --max_kmercov 1000 -k 21 -o ${outName}_genomescope_out -p 2 --fitted_hist
