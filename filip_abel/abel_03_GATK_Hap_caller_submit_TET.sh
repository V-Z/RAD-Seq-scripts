#!/bin/sh
#SBATCH --account=uio
#SBATCH --time=50:0:0
#SBATCH --mem-per-cpu=200M
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

arrayrun 1-188 abel_03_GATK_Hap_caller_worker_TET.sh   #change N of arrays based on N of samples=dirs to process
