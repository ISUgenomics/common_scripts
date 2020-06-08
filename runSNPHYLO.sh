#!/bin/bash
file=$1
module load snphylo
xvfb-run /opt/rit/spack-app/linux-rhel7-x86_64/gcc-4.8.5/snphylo-2016-02-04-qftdohupwpclptju5kxtkhlfcq3cfwux/snphylo.sh \
  -v $file
  -l 0.4 \
  -m 0.1 \
  -M  0.1 \
  -P ${file%.*} \
  -a 200000
