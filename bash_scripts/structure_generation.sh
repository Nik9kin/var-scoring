# Run only in dir dereplicator/cycloquest_minimal

#!/usr/bin/bash

print_usage ()
{
        echo "Usage:"
        echo "  ../clusterisation/structure_generation.sh <path to database> <library info file> -- print fragmented graph for every molecule from database"
        echo ""
        echo "  Run script only in dir dereplicator/cycloquest_minimal/"
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
	echo "$mol_file $name" >> $path_to_database/${lib_info}_AAGraphs1.txt
        `build/bin/print_structure --print_summary --print_rule_fragmented_graph --print_structure --structure $path_to_database/$mol_file >> $path_to_database/${lib_info}_AAGraphs1.txt`
        echo "$mol_file processed"
done < $path_to_database/$lib_info



