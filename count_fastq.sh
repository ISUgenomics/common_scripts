#!/bin/bash
# This is a bash script to count number of reads in fastq files
# 08/04/2014
# Arun Seetharam <arnstrm@iastate.edu>

function printUsage () {
    cat <<EOF

Synopsis

    $scriptName [-h | --help] [-b | --base-count] fastq_file1 [fastq_file2 fastq_file3 ..] || *.fastq[.gz]

Description

    This is a bash script counts the total number of reads in every fastq file
    It also prints the number of reads tabular and easy to read format

        -h, --help
        Brings up this help page

        -b, --base-count
        include base counts along with number of reads

	fastq_file
	A standard fastq file with any extension [but not compressed]

Author

    Arun Seetharam, Genome Informatics Facilty, Iowa State University
    arnstrm@iastate.edu
    08 April, 2014



EOF
}



if [ $# -lt 1 ] ; then
        printUsage
	exit 0
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

filear=${@};
for i in ${filear[@]}
do

if [ ! -f $i ]; then
    echo "\"$i\" file not found!"
    exit 1;
fi

if [[ $i =~ \.gz$ ]]; then
lines=$(zcat $i | wc -l |cut -d " " -f 1)
count=$(($lines / 4))
echo -n -e "\t$i : "
echo "$count"  | \
sed -r '
  :L
  s=([0-9]+)([0-9]{3})=\1,\2=
  t L'
else
lines=$(wc -l $i|cut -d " " -f 1)
count=$(($lines / 4))
echo -n -e "\t$i : "
echo "$count"  | \
sed -r '
  :L
  s=([0-9]+)([0-9]{3})=\1,\2=
  t L'
fi
done

