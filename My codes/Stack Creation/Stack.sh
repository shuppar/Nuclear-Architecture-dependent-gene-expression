#!/bin/bash

echo "Type scriptname followed by matlab programme name followed by the range of filenames. For example: Stack.sh StacksCreation.m 1 15"
fname=$1
s=$2
f=$3
fold=$4 
mkdir AvgStack      
	for i in `seq $s $f`;
        do
		
		printf "\n\033[1;31m This is image %d.tif\033[0m\r" "$i"  #To make it red.

		sed -i 7s/"1"/"$i"/g $fname
		sed -i 9s/"1"/"$i"/g $fname
		sed -i 18s/"1"/"$i"/g $fname
		sed -i 19s/"1"/"$i"/g $fname
		sed -i 20s/"1"/"$i"/g $fname
		matlab -nodesktop -nosplash < $fname
		#read -p "Press [Enter] key to start backup..."
		sed -i 7s/"$i"/"1"/g $fname
		sed -i 9s/"$i"/"1"/g $fname
		sed -i 18s/"$i"/"1"/g $fname
		sed -i 19s/"$i"/"1"/g $fname
		sed -i 20s/"$i"/"1"/g $fname

        done   
