#!/bin/bash

guppy=$(cat $1 | grep guppy | awk '{print $NF}') # guppy version
minknow=$(cat $1 | grep distribution_version| awk '{print $NF}') # MinKNOW version
bream=$(cat $1 | grep protocols_version| awk '{print $NF}') # Bream version
minknowcore=$(cat $1 | grep configuration_version| awk '{print $NF}') # #MinKNOW Core version

cat <<OUTPUT 
| Software | Version | 
| -- | -- | 
| MinKNOW	| $minknow |
| MinKNOW Core | $minknowcore |
| Bream	| $bream |
| Guppy	| $guppy |
OUTPUT
