#!/bin/sh
#SBATCH --account=uio
#SBATCH --time=50:0:0
#SBATCH --mem-per-cpu=200M
source /cluster/bin/jobsetup
module purge   # clear any inherited modules
set -o errexit # exit on errors

arrayrun 1-97 abel_00_bwa_mem_worker.sh   #change N of arrays
