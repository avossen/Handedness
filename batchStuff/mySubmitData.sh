#!/bin/bash                                                                                                                                    

#echo $d " counter is " $counter;                                                                                                              
counter=0;
subCounter=0;
dateString=`date +%d%b%Y`

if [ ! -d sub$dateString ]
then
mkdir sub$dateString
fi

if [ ! -d /pic/projects/belle/voss771/sub$dateString ]
then 
mkdir /pic/projects/belle/voss771/sub$dateString/
fi

mkdir /pic/projects/belle/voss771/handOut/sub$dateString/


for d in `cat dataEx$1_onRes.list`
do
let "counter+=1"

if [ "$counter" -le "$3" ]
then
if [ "$counter" -ge "$2" ]
then
let "subCounter+=1"
targetShFile=sub$dateString/myDataOnResJob_ex$1_$subCounter.sh
cp batchHead.sh $targetShFile
echo "module put_parameter handAna rfname\\/pic/projects/belle/voss771/sub$dateString/data_onResEx$1_$subCounter.root" >> $targetShFile
cat batchMiddle.sh >> $targetShFile
echo "process_event $d 0" >> $targetShFile
cat batchEnd.sh >> $targetShFile
fi
fi

done









