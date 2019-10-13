#!/bin/bash

fname=$1
s=$2
f=$3
fold=$4       
	for i in `seq $s $f`;
        do

		printf "\n\033[1;31m This is image %d.tif\033[0m\r" "$i"  #To make it red.

		sed -i s/"1.tif"/"$i.tif"/g $fname
		matlab -nodesktop -nosplash < $fname
		#read -p "Press [Enter] key to start backup..."
		sed -i s/"$i.tif"/"1.tif"/g $fname

        done   
