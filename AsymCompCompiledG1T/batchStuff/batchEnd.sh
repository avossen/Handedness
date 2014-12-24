

dateString=`date +%d%b%Y`

# Grab the exit code of BASF
BASFRET=$?
echo VOSS771_BASF_FINISH `date`


# Mark the file as bad if BASF returned something other than 0
if [ $BASFRET -ne 0 ]; then
  for p in /pic/projects/belle/voss771/handOut/AsymCompOut/*${SLURM_JOBID}*; do
    mv $p $p.badret
  done
fi

