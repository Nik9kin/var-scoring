#!/usr/bin/bash

print_usage ()
{
        echo "Usage:"
        echo "  ./edges_extract.sh <database dir> [<clusterisation number>] -- extract edges from .selfloop file"
}

dbd=$1	#database directory
num=$2	#clusterisation number or id

if [[ -z $dbd ]]
then
        echo "Database directory should be specified"
        print_usage
        exit 1
elif [[ !( -e $dbd) || !( -d $dbd) ]]
then
        echo "Database directory incorrect"
        print_usage
        exit 1
fi

while IFS="	" read -r v1 v2 other
do
	if [[ v1 -eq "CLUSTERID1" || $v1 -eq $v2 ]]
	then
		continue
	else
		echo "$v1 $v2" >> $dbd/clustered_spectra$num/${dbd}_edges.txt
	fi
done < $dbd/clustered_spectra$num/${dbd}_networkedges.selfloop
