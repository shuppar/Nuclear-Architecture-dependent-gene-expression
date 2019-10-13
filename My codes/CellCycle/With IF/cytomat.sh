#!/bin/bash

fname=$1
s=$2
f=$3
fold=$4      
	
	
	for i in `seq $s $f`;
        do

		printf "\n\033[1;31m This is image %d.tif\033[0m\r" "$i"  #To make it red.

		sed -i s/"1.tif"/"$i.tif"/g $fname
		sed -i s/"1.jpg"/"$i.jpg"/g $fname
        sed -i s/"1.dat"/"$i.dat"/g $fname
		matlab -nodesktop -nosplash < $fname
		#read -p "Press [Enter] key to start backup..."
		sed -i s/"$i.tif"/"1.tif"/g $fname
		sed -i s/"$i.jpg"/"1.jpg"/g $fname
        sed -i s/"$i.dat"/"1.dat"/g $fname

        done   

		
		matlab -nodesktop -nosplash - nojvm -r 'bimodefit('Intensity.dat')';

