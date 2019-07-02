#!/bin/bash
echo
echo '########################################'
echo '########### Downloading data ###########'
echo '########################################'

for i in NAD TAD NSC TSC
do
	echo
	echo
	echo "------------ Processing $i -------------"
	python py/getcases.py -t $i
	python py/casemirna.py -t $i
	python py/downdata.py -t $i
done

echo
echo "Done; Files saved in Data/"

echo
echo
echo
echo '########################################'
echo '############ Making Plots ##############'
echo '########################################'

echo
echo "Graph showing the amount of miRNAs by RNA for each Tissue Type is being generated..."
python py/getmirnasCounts.py
echo "Done; File saved in Graphs/"

echo
echo "Graph showing the number of retrieved files is being generated..."
python py/file_number.py
echo "Done; File saved in Graphs/"
echo
