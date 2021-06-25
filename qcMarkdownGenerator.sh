#!/bin/bash

# use multiqc -p . 
# the -p will generate the required plots
# then add the png folder to the multiqc_data folder and it works.


png=multiqc_data/png

cat <<MDFile > 01b_QC-results.md
# QC
## Report
### Sequence Quality Histograms
The mean quality value across each base position in the read.
![Mean Quality Score]($png/mqc_fastqc_per_sequence_quality_scores_plot_1.png)

### Per Sequence Quality Scores
The number of reads with average quality scores. Shows if a subset of reads has poor quality.
![Per Seq QC]($png/mqc_fastqc_per_sequence_quality_scores_plot_1.png)

### Per Base Sequence Content
The proportion of each base position for which each of the four normal DNA bases has been called.
![Per BAse Seq Content]($png/mqc_fastqc_per_base_sequence_quality_plot_1.png)

### Per Sequence GC Content
The average  GC content of reads. Normal random library typically have a roughly normal distribution of GC content.
![Per Seq GC content]($png/mqc_fastqc_per_sequence_gc_content_plot_Counts.png)

### Per Base N Content
The percentage of bases calles at each position for which an N was called.
![Per base N content]($png/mqc_fastqc_per_base_sequence_quality_plot_1.png)

### Sequence legth distribution

### Sequence Duplication Levels
The relative level of duplication found for every sequence.
![Seq duplicate levels]($png/mqc_fastqc_sequence_duplication_levels_plot_1.png)

### Adapter Content
The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.
![adapter content]($png/mqc_fastqc_adapter_content_plot_1.png)
MDFile

echo "## General Stats" >> 01b_QC-results.md
cat ./multiqc_data/multiqc_general_stats.txt | md  >> 01B-01b_QC-results.md

echo >> 01b_QC-results.md

echo "## FastQC Stats" >> 01b_QC-results.md
cat ./multiqc_data/multiqc_fastqc.txt |tr -s " " "-" | md  >> 01b_QC-results.md
#I had to replace spaces with `-` in order to get the correct number of feilds in awk! I used md (from common_scripts)
