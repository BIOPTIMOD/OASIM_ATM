#! /bin/bash
#
#SBATCH -N 1 -n 1              # 1 tasks on 1 node
#SBATCH --time=8:00:00         # time limits: 20 min
#SBATCH --error myJob.err      # std-error file
#SBATCH --output myJob.out     # std-output file
#SBATCH --account=OGS20_PRACE_P_1  # account number
#SBATCH --partition=m100_usr_prod  # partition to be used


#

cd $SLURM_SUBMIT_DIR

date

./launcher.sh

date
