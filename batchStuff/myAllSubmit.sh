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
mkdir $myDir
echo " dir: $myDir " ;
mkdir /pic/projects/belle/voss771/$myDir
mkdir /pic/projects/belle/voss771/handOut/$myDir


#for d in `cat mc$1_onRes.list`
for d in `cat lists/data$ex\_$res.list`
do
let "counter+=1"

#if [ "$counter" -le "$5" ]
#then
#if [ "$counter" -ge "$4" ]
#then
let "subCounter+=1"
targetShFile=$myDir/job_$subCounter.sh
#cp batchHead.sh $targetShFile
cp batchHead1.sh $targetShFile
echo "#SBATCH -o /pic/projects/belle/voss771/handOut/$myDir/jobId_$subCounter.out" >> $targetShFile
echo "#SBATCH -e /pic/projects/belle/voss771/handOut/$myDir/jobId_$subCounter.err" >> $targetShFile
echo "#SBATCH -J handedness_$subCounter"  >> $targetShFile 
cat batchHead2.sh >> $targetShFile 
echo "module put_parameter handAna rfname\\/pic/projects/belle/voss771/$myDir/job_$subCounter.root" >> $targetShFile
cat batchMiddle.sh >> $targetShFile
echo "process_event $d 0" >> $targetShFile
cat batchEnd.sh >> $targetShFile
#fi
#fi

done
done 
done








