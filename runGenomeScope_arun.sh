#!/bin/bash
module load GIF2/jellyfish/2.2.5
set +o posix
readR1=$1 
readR2=$2
kmer=$3
baseName=$4
###############################################
#paste <(ls -1 *R1*) <(ls -1 *R2*)  |awk '{print $1,$2,"21",$1}' |sed 's/\./\t/5' |awk '{print "sh runGenomeScope.sh "$1,$2,$3,$4}' >jellyfish.sh
#the above creates the sh script 
#sh runGenomeScope.sh 312_S59_R1.fq.gz 312_S59_R2.fq.gz 21 312_S59_R1
#sh runGenomeScope.sh 314_S67_R1.fq.gz 314_S67_R2.fq.gz 21 314_S67_R1

###############################################


jellyfish count -C -m $kmer -s 1000000000 -t 10 <(zcat ${readR1} ${readR2}) -o ${baseName}_K${kmer}.jf  ;\
jellyfish histo -t 16 ${baseName}_K${kmer}.jf > ${baseName}_K${kmer}.histo ; \



# once the histo file is created, visit http://qb.cshl.edu/genomescope/ website to upload the histo file

