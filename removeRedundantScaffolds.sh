#!/bin/bash

genome=$1
params_asm=${2:-"asm5"}
# Minimap other info
#- asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence

module load cdbfasta
module load bioawk
module load minimap2/r941
module load bedtools2/2.27.1-s2mtpsu

#Version information
cdbfasta -v > Versions.txt
bioawk --version >> Versions.txt
cdbfasta -v >> Versions.txt
bedtools --version >> Versions.txt


# Minimap Header info
# |Col|Type  |Description                               |
# |--:|:----:|:-----------------------------------------|
# |1  |string|Query sequence name                       |
# |2  |int   |Query sequence length                     |
# |3  |int   |Query start (0-based; BED-like; closed)   |
# |4  |int   |Query end (0-based; BED-like; open)       |
# |5  |char  |Relative strand: "+" or "-"               |
# |6  |string|Target sequence name                      |
# |7  |int   |Target sequence length                    |
# |8  |int   |Target start on original strand (0-based) |
# |9  |int   |Target end on original strand (0-based)   |
# |10 |int   |Number of residue matches                 |
# |11 |int   |Alignment block length                    |
# |12 |int   |Mapping quality (0-255; 255 for missing)  |

# Minimap other info
#- asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence


# Create cdbfasta databases

cdbfasta $genome

# Separate out contigs <1M and >1M
bioawk -c fastx '{print $name,length($seq)}' $genome |awk '$2<1000000' |wc
bioawk -c fastx '{print $name,length($seq)}' $genome |awk '$2>1000000' |wc

bioawk -c fastx '{print $name,length($seq)}' $genome |awk '$2<1000000 {print $1}' |cdbyank $genome >LittleContigs.fasta
bioawk -c fastx '{print $name,length($seq)}' $genome  |awk '$2>999999 {print $1}' |cdbyank $genome >BigContigs.fasta

# Align small contigs to large

minimap2 -x $params_asm -t 36 BigContigs.fasta   LittleContigs.fasta > Little2Big_minimap.paf

# create lengths for littleContigs

seqlen.awk LittleContigs.fasta > LittleContigs.length

# How many little contigs in the assembly did not map at all to the big contigs.



grep ">" $genome | perl -pe 's/\>//g' > AllContigIDs.txt
grep ">" LittleContigs.fasta  | perl -pe 's/\>//g' >LittleContigsIDs.txt
seqlen.awk LittleContigs.fasta >LittleContigsIDs.len
grep ">" BigContigs.fasta  | perl -pe 's/\>//g' >BigContigsIDs.txt
seqlen.awk BigContigs.fasta >BigContigsIDs.len
awk '{print $1}' Little2Big_minimap.paf  | sort | uniq  > AlignedContigIDs.txt


wc -l LittleContigsIDs.txt | awk '{print "LittleContigs: "$0}' > Log_${params_asm}.txt
wc -l BigContigsIDs.txt | awk '{print "BigContigsIDs: "$0}' >> Log_${params_asm}.txt
wc -l AllContigIDs.txt | awk '{print "AllContigIDs: "$0}' >> Log_${params_asm}.txt
wc -l AlignedContigIDs.txt | awk '{print "AlignedContigIDs: "$0}' >> Log_${params_asm}.txt

# Grab those scaffold Ids that did not align

cat AlignedContigIDs.txt LittleContigsIDs.txt LittleContigsIDs.txt | sort | uniq -c | awk '$1==2 {print $2}' > NotAlignedContigIDs.txt
cat AlignedContigIDs.txt LittleContigsIDs.txt LittleContigsIDs.txt | sort | uniq -c | awk '$1==2 {print $2}' | cdbyank $genome > NotAlignedContigs.fasta

wc -l NotAlignedContigIDs.txt | awk '{print "NotAlignedContigs: "$0}' >> Log_${params_asm}.txt


# Remove Redundant scaffolds with small gaps (100) with > 90% alignment length

cat Little2Big_minimap.paf |awk '{print $1,$3,$4}' |tr " " "\t" |sort -k1,1V -k2,3n |bedtools merge -d 100 -i - >LittleContigRedundancy100.bed
cat LittleContigRedundancy100.bed | awk '{print $1}' | xargs -I xx awk '$1=="'xx'"' LittleContigs.length > LittleContigRedundancy100.lengths
paste LittleContigRedundancy100.bed LittleContigRedundancy100.lengths |awk '{print $1,$2,$3,$5}' |awk '(($3-$2)>(.9*$4)) {print $0,($3-$2)/$4}' | awk '{print $1}' | sort | uniq |cdbyank $genome > RedundantScaffolds.fasta
paste LittleContigRedundancy100.bed LittleContigRedundancy100.lengths |awk '{print $1,$2,$3,$5}' |awk '(($3-$2)>(.9*$4)) {print $0,($3-$2)/$4}' | awk '{print $1}' | sort | uniq |wc -l
grep ">" RedundantScaffolds.fasta | perl -pe 's/\>//g' > redudantScaffolds.txt
cat redudantScaffolds.txt LittleContigsIDs.txt LittleContigsIDs.txt | sort | uniq -c | awk '$1==2 {print $2}' > Pass1_RemainingScaffolds.txt


wc -l Pass1_RemainingScaffolds.txt | awk '{print "Pass1_RemainingScaffolds: "$0}' >> Log_${params_asm}.txt


# Remove Medium sized or larger scaffolds with medium sized gaps (1000) (scaffoldSize>2999) with > 90% alignment length

cat Pass1_RemainingScaffolds.txt | xargs -I xx grep xx  Little2Big_minimap.paf |awk '{print $1,$3,$4}' |tr " " "\t" |sort -k1,1V -k2,3n |bedtools merge -d 1000 -i - >LittleContigRedundancy1000.bed
cat LittleContigRedundancy1000.bed | awk '{print $1}' | xargs -I xx awk '$1=="'xx'"' LittleContigs.length   > LittleContigRedundancy1000.lengths
paste LittleContigRedundancy1000.bed LittleContigRedundancy1000.lengths |awk '{print $1,$2,$3,$5}' |awk '(($3-$2)>(.9*$4)) && ($4>2999) {print $0,($3-$2)/$4}' |awk '{print $1}' |sort |uniq |cdbyank $genome >> RedundantScaffolds.fasta

# number of scaffolds identified as redundant
# paste LittleContigRedundancy1000.bed LittleContigRedundancy1000.lengths |awk '{print $1,$2,$3,$5}' |awk '(($3-$2)>(.9*$4)) && ($4>2999) {print $0,($3-$2)/$4}' |awk '{print $1}' |sort |uniq | wc | awk '{print "NotAlignedContigs: "$0}'   > Log_${params_asm}.txt

# Grab scaffolds that I couldn't mark as redundant and run another filter
# paste LittleContigRedundancy1000.bed LittleContigRedundancy1000.lengths |awk '{print $1,$2,$3,$5}' |awk '(($3-$2)<(.9*$4))  {print $0,($3-$2)/$4}' |awk '{print $1}' |sort |uniq > Pass2_RemainingScaffolds.txt

grep ">" RedundantScaffolds.fasta | perl -pe 's/\>//g' > redudantScaffolds.txt
cat redudantScaffolds.txt LittleContigsIDs.txt LittleContigsIDs.txt | sort | uniq -c | awk '$1==2 {print $2}' > Pass2_RemainingScaffolds.txt
wc -l Pass2_RemainingScaffolds.txt | awk '{print "Pass2_RemainingScaffolds: "$0}' >> Log_${params_asm}.txt




# Larger sized scaffolds with slightly larger sized gaps (4000) (scaffolds>10000)

cat Pass2_RemainingScaffolds.txt | xargs -I xx grep xx  Little2Big_minimap.paf |awk '{print $1,$3,$4}' |tr " " "\t" |sort -k1,1V -k2,3n |bedtools merge -d 4000 -i - >LittleContigRedundancy4000.bed
cat LittleContigRedundancy4000.bed | awk '{print $1}' | xargs -I xx awk '$1=="'xx'"' LittleContigs.length   > LittleContigRedundancy4000.lengths

paste LittleContigRedundancy4000.bed LittleContigRedundancy4000.lengths |awk '{print $1,$2,$3,$5}' |awk '(($3-$2)>(.9*$4)) && ($4>9999) {print $0,($3-$2)/$4}' |awk '{print $1}' |sort |uniq |cdbyank $genome >> RedundantScaffolds.fasta
# paste LittleContigRedundancy4000.bed LittleContigRedundancy4000.lengths |awk '{print $1,$2,$3,$5}' |awk '(($3-$2)>(.9*$4)) && ($4>2999) {print $0,($3-$2)/$4}' |awk '{print $1}' |sort |uniq | wc

grep ">" RedundantScaffolds.fasta | perl -pe 's/\>//g' > redudantScaffolds.txt
cat redudantScaffolds.txt LittleContigsIDs.txt LittleContigsIDs.txt | sort | uniq -c | awk '$1==2 {print $2}' > Pass3_RemainingScaffolds.txt
wc -l Pass3_RemainingScaffolds.txt | awk '{print "Pass3_RemainingScaffolds: "$0}' >> Log_${params_asm}.txt

cat Pass3_RemainingScaffolds.txt |cdbyank $genome > Pass3_RemainingScaffolds.fasta
seqlen.awk Pass3_RemainingScaffolds.fasta >Pass3_RemainingScaffolds.len

# print log and version to screen
cat Versions.txt
cat Log_${params_asm}.txt
# grab version information
#cdbfasta -v >> Log_${params_asm}.txt
#bedtools --version >> Log_${params_asm}.txt
#bioawk --version  | awk '{print "bio"$0}'>> Log_${params_asm}.txt
#minimap2 --version | awk '{print "minimap2: "$0}'>> Log_${params_asm}.txt


