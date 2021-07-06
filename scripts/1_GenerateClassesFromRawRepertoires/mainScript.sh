#!/bin/bash
let "k=0"
#for i in 5*/
for i in */
do
	#if [ -f  $i ]
	#then  #only if it is a file, not a folder
	let "k=k+1"
	echo "Processing ${i%/}"
	echo "nr  ${k}"
	#Stupid linux needs a space after if
	if [[ $k == 35 ]]; then
		#echo "10x"
		let "k = 0"
		nice -n 19 Rscript mainScript.R "${i%/}" >& newRes${i%/}.txt
	else
		#echo "11x"
		nice -n 19 Rscript mainScript.R "${i%/}" >& newRes${i%/}.txt &
	fi
done

