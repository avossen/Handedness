#!/bin/bash                                                                                                                                    

#echo $d " counter is " $counter;                                                                                                              
counter=0;
subCounter=0;
dateString=`date +%d%b%Y`

mkdir sub$dateString
mkdir /pic/projects/belle/voss771/sub$dateString/
mkdir /pic/projects/belle/voss771/handOut/sub$dateString/


#for d in `cat mc$1_onRes.list`
for d in `cat $1_$2_$3.list`
do
let "counter+=1"

if [ "$counter" -le "$5" ]
then
if [ "$counter" -ge "$4" ]
then
let "subCounter+=1"
targetShFile=sub$dateString/$1_$2_$3_$subCounter.sh
cp batchHead.sh $targetShFile
echo "module put_parameter handAna rfname\\/pic/projects/belle/voss771/sub$dateString/$1_$2_$3_$subCounter.root" >> $targetShFile
cat batchMiddle.sh >> $targetShFile
echo "process_event $d 0" >> $targetShFile
cat batchEnd.sh >> $targetShFile
fi
fi

done









