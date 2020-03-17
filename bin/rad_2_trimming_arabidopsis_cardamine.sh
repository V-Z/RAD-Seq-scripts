#!/bin/bash

# qsub -l walltime=48:0:0 -l select=1:ncpus=4:mem=6gb:scratch_local=300gb -q ibot -m abe ~/bin/rad_2_trimming_arabidopsis_cardamine.sh

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a $SCRATCHDIR $DATADIR/ && clean_scratch' TERM

# Location of data to be trimmed
DATADIR='/storage/ostrava2-archive/tape_tape/backup/VO_cuni_prf_arab/shared/rad'
# Library to process
# run_1_2_dip_tet_harvard run_3_4_dip_tet_2017_01_arenosa_lyrata run_5_dip_tet_2017_11_lyrata run_6_7_dip_tet_2018_01_lyrata run_8_dip_tet_arenosa run_9_dip_tet_lyrata run_10_dip_tet_lyrata run_11_dip_tet_arenosa run_12_dip_tet_arenosa_lyrata
# cardamine_run_01_test cardamine_run_02_test
LIBRARY='cardamine_run_02_test'

# Change working directory
cd "$SCRATCHDIR"/ || exit 1

# Launch it
module add parallel-20160622 || exit 1
module add trimmomatic-0.36 || exit 1

# Prepare the task
cp -a /storage/praha1/home/"$LOGNAME"/rad/2_trimming/* "$SCRATCHDIR"/ || exit 1
cp -a "$DATADIR"/"$LIBRARY"/1_demultiplexed "$SCRATCHDIR"/ || exit 1

./radseq_2_trimming.sh -f 1_demultiplexed -c 3 -o trimmed -a adaptors.fa -m 5g -t "$TRIMMOMATIC_BIN" | tee trimming.log

# Copy results back to storage
cp -a trimmed trimming.log "$DATADIR"/"$LIBRARY"/ || export CLEAN_SCRATCH='false'

exit

