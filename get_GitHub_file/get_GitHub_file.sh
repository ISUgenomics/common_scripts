#!/bin/bash

echo "----------------"
echo "USAGE:"
echo "       . ./get_GitHub_folder.sh <URL to the GitHub file>"
echo "e.g.,  . ./get_GitHub_folder.sh https://github.com/ISUgenomics/common_scripts/blob/master/get_GitHub_file/get_GitHub_file.sh"
echo "----------------"
echo ""

URL=`echo $1`
folder=`echo $URL | sed 's|blob/master|trunk|g'`

svn export $folder