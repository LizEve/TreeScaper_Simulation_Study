#!/bin/bash

# This will iterate through all files with designated suffix. 
# Output will be a text file with all the commands needed to convert nexus files into a format readable by treescaper. 
# This script executes the text file to add taxa translate blocks and create nexus files. 
# Edit directory for Paup* as needed. 

file="paup_to_nexus_list.txt"

if [ -f $file ] ; then
    rm $file
fi

echo "begin paup;
	set nowarnreset autoclose maxtrees = 200000;
" >> paup_to_nexus_list.txt

for f in *.tree
do
base=`basename $f .tree`
out_nexus=$base".nex"
echo "
	execute $f;
	savetrees file=$out_nexus;
" >> paup_to_nexus_list.txt
done

echo "
	quit;
" >> paup_to_nexus_list.txt

#/Applications/paup4a157_osx paup_to_nexus_list.txt > paup.out
/Applications/PAUP_dev_icc.app paup_to_nexus_list.txt > paup.out

mkdir trees
mkdir startingtrees
mv *starting_tree* startingtrees
mv *.tree trees
mv *.log trees

#rm *.tree
