#!/bin/bash

# Author: Vojtěch Zeisek, https://trapa.cz/
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# 

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# qsub -l walltime=24:0:0 -l select=1:ncpus=4:mem=4gb:scratch_local=400gb -q ibot -m abe ~/radseq/bin/rad_1_demultiplexing_1_qsub.sh

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a ${SCRATCHDIR} ${DATADIR}/ && clean_scratch' TERM

# Location of data to be demultiplexed
DATADIR='/auto/pruhonice1-ibot/shared/brassicaceae/rad'
# Library to process
# run_03_04_dip_tet_2017_01_arenosa_lyrata run_05_dip_tet_2017_11_lyrata run_06_07_dip_tet_2018_01_lyrata run_08_dip_tet_2018_08_arenosa run_09_dip_tet_2019_01_lyrata run_10_dip_tet_2019_05_lyrata run_11_dip_tet_2019_06_arenosa run_12_dip_tet_2019_12_arenosa_lyrata run_2018_06_tet_arenosa
# cardamine_run_01_test cardamine_run_02_test
LIBRARY='run_2018_06_tet_arenosa'
# Table with demultiplexing key
# samples_demultiplexing_list_run_03_04_dip_tet_2017_01_arenosa_lyrata.tsv samples_demultiplexing_list_run_05_dip_tet_2017_11_lyrata.tsv samples_demultiplexing_list_run_06_07_dip_tet_2018_01_lyrata.tsv samples_demultiplexing_list_run_08_dip_tet_2018_08_arenosa.tsv samples_demultiplexing_list_run_09_dip_tet_2019_01_lyrata.tsv samples_demultiplexing_list_run_10_dip_tet_2019_05_lyrata.tsv samples_demultiplexing_list_run_11_dip_tet_2019_06_arenosa.tsv samples_demultiplexing_list_run_12_dip_tet_2019_12_arenosa_lyrata.tsv samples_demultiplexing_list_run_miseq_tet_2018_06_arenosa.tsv
# samples_demultiplexing_list_cardamine_run_01_test.tsv samples_demultiplexing_list_cardamine_run_02_test.tsv
TABLE='samples_demultiplexing_list_run_miseq_tet_2018_06_arenosa.tsv'

# Change working directory
echo "Going to working directory ${SCRATCHDIR}"
cd "${SCRATCHDIR}"/ || exit 1
echo

# Required modules
echo "Loading modules"
module add parallel-20200322 || exit 1
module add fastx-0.0.14 || exit 1
echo

# Copy data
echo "Copying..."
echo "Scripts etc. - /storage/praha1/home/${LOGNAME}/radseq/"
cp -a /storage/praha1/home/"${LOGNAME}"/radseq/{bin/rad_1_demultiplexing_2_run.sh,demultiplexing_lists/"${TABLE}"} "${SCRATCHDIR}"/ || exit 1
echo "Data to process - ${DATADIR}/${LIBRARY}"
cp -a "${DATADIR}"/"${LIBRARY}"/0_raw_illumina "${SCRATCHDIR}"/ || exit 1
echo

# Running the task
echo "Preprocessing the FASTQ files..."
./rad_1_demultiplexing_2_run.sh -s "${TABLE}" -c 4 -o 1_demultiplexed -x txt.bz2 -f 0_raw_illumina -u bzcat | tee demultiplexing.log
echo

# Remove unneeded file
echo "Removing unneeded files"
rm -rf 0_raw_illumina rad_1_demultiplexing_2_run.sh "${TABLE}"
echo

# Copy results back to storage
echo "Copying results back to ${DATADIR}/${LIBRARY}"
cp -a "${SCRATCHDIR}" "${DATADIR}"/"${LIBRARY}"/ || export CLEAN_SCRATCH='false'
echo

exit

