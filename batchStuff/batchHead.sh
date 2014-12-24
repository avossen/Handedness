#!/bin/bash
#SBATCH -A belle        # this is the account to which time is charged --
#                       # use `belle' for belle analysis 
#SBATCH -t 12:00:00     # time limit for job in HH:MM:SS
#SBATCH -N 1            # number of CPU cores requested
#
# Redirect the job's stdout
#SBATCH -o /pic/projects/belle/voss771/handOut/jobId_%j.out
#
# Redirect the job's stderr
#SBATCH -e /pic/projects/belle/voss771/handOut/jobId_%j.err
#
# Set the job name to be `ypipi-e000053r000077-b20090127_0910.mdst'
#SBATCH -J handedness

# `.brofile' is my equivalent to `.bashrc'
echo Zuhause: $HOME
source $HOME/.bashrc

# This sets up your environment to run code in the Belle paradigm
source /pic/projects/belle/scripts/general_env.sh

# Add my modifications to the standard Belle environment
source $HOME/custom_env.sh

# This is necessary to get code that fills HBOOK ntuples running comfortably
# on Olympus
shopt -s expand_aliases
belleset_hbook32768_libs

echo ANSELM_BASF_START `date`
export USE_GRAND_REPROCESS_DATA=1
export BELLE_MESSAGE_LEVEL=INFO
export LD_LIBRARY_PATH=/people/voss771/handedness/:`root-config --libdir`:$LD_LIBRARY_PATH
export BASF_MODULE_DIR=/people/voss771/handedness/:$BASF_MODULE_DIR
#setenv BASF_NPROCESS 0

echo    $LD_LIBRARY_PATH
echo    $BASF_MODULE_DIR
echo    $PANTHER_TABLE_DIR
MODULE=handAna
HBKFILE=/dev/null
# I put the BASF script inline like so
basf <<EOF
path create main
path create analysis
path add_module main fix_mdst
path add_condition main >:0:analysis
path add_condition main <=:0:KILL
path add_module analysis ${MODULE}
