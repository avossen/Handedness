#!/bin/bash
#SBATCH -A belle        # this is the account to which time is charged --
#                       # use `belle' for belle analysis 
#SBATCH -t 12:00:00     # time limit for job in HH:MM:SS
#SBATCH -N 1            # number of CPU cores requested
#SBATCH -o /pic/projects/belle/voss771/handOut/AsymCompOut/O_subMC_ex19_on_resonance_uds.out
#SBATCH -e /pic/projects/belle/voss771/handOut/AsymCompOut/O_subMC_ex19_on_resonance_uds.err
#SBATCH -J AsymComp_subMC_ex19_on_resonance_uds

# `.brofile' is my equivalent to `.bashrc'
echo Zuhause: $HOME
#source $HOME/.bashrc

# This sets up your environment to run code in the Belle paradigm
#source /pic/projects/belle/scripts/general_env.sh

# Add my modifications to the standard Belle environment
#source $HOME/custom_env.sh

module purge
source $HOME/.bashrc

# This is necessary to get code that fills HBOOK ntuples running comfortably
# on Olympus
shopt -s expand_aliases
belleset_hbook32768_libs

echo ANSELM_BASF_START `date`
export USE_GRAND_REPROCESS_DATA=1
export BELLE_MESSAGE_LEVEL=INFO
export LD_LIBRARY_PATH=/people/voss771/handedness/AsymCompiledG1T:`root-config --libdir`:$LD_LIBRARY_PATH

#setenv BASF_NPROCESS 0

echo    $LD_LIBRARY_PATH

/people/voss771/handedness/AsymCompCompiledG1T/TwoHadAsymsCMod /pic/projects/belle/voss771/subMC_ex19_on_resonance_uds mc
cp *19*uds*.root /pic/projects/belle/voss771/AsymmetriesMC/


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

