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
echo "Going to working directory ${SCRATCHDIR}"
cd "${SCRATCHDIR}"/ || exit 1
echo

# Required modules
echo "Loading modules"
module add trimmomatic-0.36 || exit 1
module add bbmap-38.42 || exit 1
module add fastQC-0.11.5 || exit 1
module add parallel-20160622 || exit 1
echo

# Copy data
echo "Copying..."
echo "Scripts etc. - /storage/praha1/home/${LOGNAME}/rad/"
cp -a /storage/praha1/home/"$LOGNAME"/rad/2_trimming/* "$SCRATCHDIR"/ || exit 1
echo "Data to process - ${DATADIR}/${LIBRARY}"
cp -a "$DATADIR"/"$LIBRARY"/1_demultiplexed "$SCRATCHDIR"/ || exit 1
echo

# Running the task
echo "Preprocessing the FASTQ files..."
./hybseq_1_prep_2_run.sh -f 0_data -c 4 -o 1_trimmed -d 2_dedup -q 3_qual_rep -a adaptors.fa -m 4 -t "${TRIMMOMATIC_BIN}" | tee hybseq_prepare.log # HybSeq
./radseq_2_trimming.sh -f 1_demultiplexed -c 3 -o trimmed -a adaptors.fa -m 5g -t "$TRIMMOMATIC_BIN" | tee trimming.log
echo

# Copy results back to storage
cp -a trimmed trimming.log "$DATADIR"/"$LIBRARY"/ || export CLEAN_SCRATCH='false'

exit

# Remove unneeded file
echo "Removing unneeded files"
rm adaptors.fa hybseq_1_prep_2_run.sh
echo

# Copy results back to storage
echo "Copying results back to ${DATADIR}"
cp -a "${SCRATCHDIR}" "${DATADIR}"/ || export CLEAN_SCRATCH='false'
echo
