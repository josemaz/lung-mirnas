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
	python3 py/getcases.py -t $i
	python3 py/casemirna.py -t $i
	python3 py/downdata.py -t $i
done

echo
echo "Done; Files saved in Data/(Adeno|Squamous)/(NAD|TAD|NSC|TSC)"

echo
echo
echo
echo '########################################'
echo '############ Making Plots ##############'
echo '########################################'

echo
echo "Graph showing the amount of miRNAs by RNA for each Tissue Type is being generated..."
python3 py/getmirnasCounts.py
echo "Done; File saved in Plots/"

echo
echo "Graph showing the number of retrieved files is being generated..."
python3 py/file_number.py
echo "Done; File saved in Plots/"
echo
