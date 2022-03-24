#!/bin/bash

echo "----------------"
echo "USAGE:"
echo "       . ./get_GitHub_folder.sh <URL to the GitHub folder>"
echo "       . ./get_GitHub_folder.sh https://github.com/ISUgenomics/common_scripts/tree/master/get_GitHub_file"
echo "----------------"
echo ""

URL=`echo $1`
folder=`echo $URL | sed 's|tree/master|trunk|g'`

svn export $folder