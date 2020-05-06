#!/bin/bash

# qsub -l walltime=24:0:0 -l select=1:ncpus=3:mem=4gb:scratch_local=400gb -q ibot -m abe ~/bin/rad_1_demultiplexing_arabidopsis_cardamine.sh

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a $SCRATCHDIR $DATADIR/ && clean_scratch' TERM

# Location of data to be demultiplexed
DATADIR='/auto/pruhonice1-ibot/shared/brassicaceae/rad'
# Library to process
# run_3_4_dip_tet_2017_01_arenosa_lyrata run_5_dip_tet_2017_11_lyrata run_6_7_dip_tet_2018_01_lyrata run_8_dip_tet_arenosa run_9_dip_tet_lyrata run_10_dip_tet_lyrata run_11_dip_tet_arenosa run_12_dip_tet_arenosa_lyrata
# cardamine_run_01_test cardamine_run_02_test
LIBRARY='cardamine_run_02_test'
# Table with demultiplexing key
# samples_demultiplexing_list_run_3_4.tsv samples_demultiplexing_list_run_5_lyrata.tsv samples_demultiplexing_list_run_6_7_lyrata.tsv samples_demultiplexing_list_run_8_arenosa.tsv samples_demultiplexing_list_run_9_lyrata.tsv samples_demultiplexing_list_run_10_lyrata.tsv samples_demultiplexing_list_run_11_arenosa.tsv samples_demultiplexing_list_run_12_arenosa_lyrata.tsv
# samples_demultiplexing_list_cardamine_run_01_test.tsv samples_demultiplexing_list_cardamine_run_02_test.tsv
TABLE='samples_demultiplexing_list_cardamine_run_02_test.tsv'

# Change working directory
echo "Going to working directory ${SCRATCHDIR}"
cd "${SCRATCHDIR}"/ || exit 1
echo

# Required modules
echo "Loading modules"
module add parallel-20160622 || exit 1
module add fastx-0.0.14 || exit 1
echo

# Copy data
echo "Copying..."
echo "Scripts etc. - /storage/praha1/home/${LOGNAME}/rad/"
cp -a /storage/praha1/home/"${LOGNAME}"/rad/1_demultiplexing/* "${SCRATCHDIR}"/ || exit 1
echo "Data to process - ${DATADIR}/${LIBRARY}"
cp -a "${DATADIR}"/"${LIBRARY}"/0_raw_illumina "${SCRATCHDIR}"/ || exit 1
echo

# Running the task
echo "Preprocessing the FASTQ files..."
./radseq_1_demultiplexing.sh -s "${TABLE}" -c 3 -o demultiplexed -f 0_raw_illumina -x fastq.bz2 -u bzcat | tee demultiplexing.log
echo

# Copy results back to storage
echo "Copying results back to ${DATADIR}/${LIBRARY}"
cp -a demultiplexed "${DATADIR}"/"${LIBRARY}"/ || export CLEAN_SCRATCH='false'
echo

exit

