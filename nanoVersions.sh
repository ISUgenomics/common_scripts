#!/bin/bash

# nanoVerions.sh report_FAO68032_20210219_1551_cc50bf20.md 
# you will need to find the markdown FAO report file.

guppy=$(cat $1 | grep guppy | awk '{print $NF}') # guppy version
minknow=$(cat $1 | grep distribution_version| awk '{print $NF}') # MinKNOW version
bream=$(cat $1 | grep protocols_version| awk '{print $NF}') # Bream version
minknowcore=$(cat $1 | grep configuration_version| awk '{print $NF}') # #MinKNOW Core version
flowcell=$(cat $1 | grep exp_script_name | awk '{print $NF}')

cat <<OUTPUT 
| Software | Version | 
| -- | -- | 
| MinKNOW	| $minknow |
| MinKNOW Core | $minknowcore |
| Bream	| $bream |
| Guppy	| $guppy |
| flowcell | $flowcell | 
OUTPUT
