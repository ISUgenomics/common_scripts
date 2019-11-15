#!/bin/bash
# This script compares the data in 3 columns, priting the column 2 and number of character mismatches the column 2 has
# compared to the first column
# Arun Seetharam <arnstrm@iastate.edu>
# 03/29/2019

function printUsage () {
    cat <<EOF

Synopsis

    $scriptName [-h | --help] tab-delimited-2-column-text file

Description

    This is a awk script that counts the character mis-matches between 2 strings
    It requires a tab-delimited file with 2 columns, with first column having the query string
    and the second column with target string. The script pritns the target string and the number of mis-matches
    the target string has compared to the query string. 


        -h, --help
        Brings up this help page


       	tab-delimited 2-column text file
        A text file with any extension (un-compressed), with 2 columns. First column query and second target.

Author

    Arun Seetharam, Genome Informatics Facilty, Iowa State University
    arnstrm@iastate.edu
    29 March, 2019



EOF
}

if [ $# -lt 1 ] ; then
        printUsage
        exit 0
fi

file=$1

# sanity checks
# see if file exists
if [ ! -f $file ]; then
    echo "\"$i\" file not found!"
    exit 1;
fi

# see if it has uneven number of columns

col1=$(awk '{print NF}' $file | sort | uniq |wc -l)

if [ ${col1} -gt 1 ]; then
    echo "check file format, some rows have more/less than 2 columns!"
    exit 1;

fi

# check if they have exactly 2 columns
col2=$(awk '{print NF}' $file | sort | uniq | head -n 1)


if [ ${col2} -ne 2 ]; then
    echo "check file format, there should be exactly 2 columns!"
    exit 1;

fi


while :
do
    case $1 in
        -h | --help | -\?)
            printUsage
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            printf >&2 'WARNING: Unknown option (ignored): %s\n' "$1"
            shift
            ;;
        *)
            break
            ;;
    esac
done


awk '{
     {mm=0} max=(length($1) >= length($2))? length($1): length($2)
     for(i=1; i <= max; i++)
     {
     	 v1=substr($1, i, 1)
     	 v2=substr($2, i, 1)
     	 if(v1 != v2){ mm++ }
     }
    }
    { print $2"\t"mm }' $file
