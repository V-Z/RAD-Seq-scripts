#!/bin/bash

# qsub -l walltime=4:0:0 -l select=1:ncpus=1:mem=100mb:scratch_local=50gb -q ibot -m abe ~/bin/rad_3_mapping_stats_arabidopsis_cardamine.sh

# Set data directories
DATADIR="/storage/praha1/home/$LOGNAME"

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a $SCRATCHDIR $DATADIR/ && clean_scratch' TERM

# Change working directory
cd "$SCRATCHDIR"/ || exit 1

# Launch it
module add samtools-1.9 || exit 1

# Prepare the task
cp -a /storage/ostrava2-archive/tape_tape/backup/VO_cuni_prf_arab/shared/for_stats "$SCRATCHDIR"/ || exit 1

# Calculating the statistics
find for_stats/ -name "*.rg.bam" -exec echo '{}' >> mapping_stats.txt \; -exec samtools flagstat '{}' >> mapping_stats.txt \; -exec echo >> mapping_stats.txt \;

# Copy results back to storage
cp mapping_stats.txt /storage/ostrava2-archive/tape_tape/backup/VO_cuni_prf_arab/shared/for_stats/ || export CLEAN_SCRATCH='false'

exit

