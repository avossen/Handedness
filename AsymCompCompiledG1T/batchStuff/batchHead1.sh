#!/bin/bash
#SBATCH -A belle        # this is the account to which time is charged --
#                       # use `belle' for belle analysis 
#SBATCH -t 12:00:00     # time limit for job in HH:MM:SS
#SBATCH -N 1            # number of CPU cores requested
