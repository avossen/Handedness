#!/bin/bash                                                                                                                                    

#echo $d " counter is " $counter;                                                                                                              
counter=0;
subCounter=0;
dateString=`date +%d%b%Y`



for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
for res in on_resonance continuum
do


myDir=subData_ex$ex\_$res
echo " dir: $myDir " ;
mkdir /pic/projects/belle/voss771/$myDir
mkdir /pic/projects/belle/voss771/handOut/$myDir


done 
done








