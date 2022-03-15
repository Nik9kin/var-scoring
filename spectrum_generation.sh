# Run only in dir dereplicator/cycloquest_minimal

#!/usr/bin/bash

print_usage ()
{
	echo "Usage:"
        echo "  ../clusterisation/spectrum_generation.sh <path to database> <library info file> -- print mgf spectrum for every molecule from database"
	echo ""
	echo "	Run script only in dir dereplicator/cycloquest_minimal/"
}

path_to_database=$1
lib_info=$2

if [[ -z $path_to_database ]]
then
	echo "Path to Database should be specified"
	print_usage
	exit 1
elif [[ !( -e $path_to_database) || !( -d $path_to_database) ]]
then
	echo "Path to Database incorrect"
	print_usage
	exit 1
fi

if [[ -z $lib_info ]]
then
	echo "Library info file should be specified"
	print_usage
	exit 2
elif [[ !( -e $path_to_database/$lib_info) || !( -f $path_to_database/$lib_info) ]]
then
	echo "Library info file incorrect"
	print_usage
	exit 2
fi

while IFS=" " read -r mol_file name other
do
	`build/bin/print_structure --print_mgf_spectrum --break_double --structure $path_to_database/$mol_file >> $path_to_database/$lib_info.all_spectrums.mgf`
	echo "$mol_file processed, pnp name: $name"
done < $path_to_database/$lib_info



