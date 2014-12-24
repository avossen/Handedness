#!/bin/bash                                                                                                                                    

#echo $d " counter is " $counter;                                                                                                              
counter=0;
subCounter=0;
dateString=`date +%d%b%Y`

mkdir /pic/projects/belle/voss771/AsymmetriesData/
mkdir /pic/projects/belle/voss771/handOut/
mkdir /pic/projects/belle/voss771/handOut/AsymCompOutData

for ex in 07 09 11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73
do
for res in on_resonance continuum
do


myDir=subData_ex$ex\_$res
myOutDir=/pic/projects/belle/voss771/handOut/AsymCompOut/$myDir
dataDir=/pic/projects/belle/voss771/$myDir

mkdir $myOutDir

echo " dir: $myDir " ;
echo " out dir: $myOutDir " ;

targetShFile=job_DataEx$ex\_$res.sh
#cp batchHead.sh $targetShFile
cp batchHead1.sh $targetShFile
echo "#SBATCH -o /pic/projects/belle/voss771/handOut/AsymCompOut/O_$myDir.out" >> $targetShFile
echo "#SBATCH -e /pic/projects/belle/voss771/handOut/AsymCompOut/O_$myDir.err" >> $targetShFile
echo "#SBATCH -J AsymComp_$myDir"  >> $targetShFile 
cat batchHead2.sh >> $targetShFile 
echo "/people/voss771/handedness/AsymCompCompiledG1T/TwoHadAsymsCMod $dataDir" >>$targetShFile
echo "find . -iname '*$ex*.root' -amin -60 -print0 " >> $targetShFile 
echo "find . -iname '*$ex*.root' -amin -60 -exec cp {} /pic/projects/belle/voss771/AsymmetriesData/ \;" >> $targetShFile 
cat batchEnd.sh >> $targetShFile
#fi
#fi
done 
done








