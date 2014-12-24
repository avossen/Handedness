#!/bin/bash                                                                                                                                    
#use as argument e.g. "sub14Apr/data"
rm fpda*
declare -i num_jobs
declare -i max_jobs
declare -i exit_code

max_jobs=5000
num_jobs=${max_jobs}  # immediate reset below

for d in "${1}"*.sh ; do
  if [[ ${num_jobs} -ge ${max_jobs} ]] ; then  # Reset when we reach ${max_jobs}
    num_jobs=`squeue -p shared| grep voss | wc -l` 
  fi

  while [[ ${num_jobs} -ge ${max_jobs} ]] ;  do 
    echo " waiting 20 secs, we have to many jobs in the queue "
    sleep 20;
    num_jobs=`squeue -p shared| grep voss | wc -l`
  done

  sbatch -n 1 -N 1 -p shared ${d}
  exit_code=${?}
  if [[ ${exit_code} -ne 0 ]] ; then 
    echo ${d} "Not submitted. Waiting 20 seconds, before moving on"
    sleep 20
    num_jobs=`squeue -p shared| grep voss | wc -l` # paranoid
  fi

  let num_jobs=num_jobs+1

done

