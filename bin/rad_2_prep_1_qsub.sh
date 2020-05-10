#!/bin/bash

# Author: Vojtěch Zeisek, https://trapa.cz/
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# Process all files in the working directory (named *.R[12].fq[.bz2]) to be ready for mapping and calling haplotypes (trimming, deduplication, statistics, quality checks).

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# qsub -l walltime=24:0:0 -l select=1:ncpus=4:mem=8gb:scratch_local=400gb -q ibot -m abe ~/radseq/bin/rad_2_prep_1_qsub.sh

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a ${SCRATCHDIR} ${DATADIR}/ && clean_scratch' TERM

# Location of data to be trimmed
DATADIR='/auto/pruhonice1-ibot/shared/brassicaceae/rad'
# Library to process
# run_1_2_dip_tet_harvard run_3_4_dip_tet_2017_01_arenosa_lyrata run_5_dip_tet_2017_11_lyrata run_6_7_dip_tet_2018_01_lyrata run_8_dip_tet_arenosa run_9_dip_tet_lyrata run_10_dip_tet_lyrata run_11_dip_tet_arenosa run_12_dip_tet_arenosa_lyrata
# cardamine_run_01_test cardamine_run_02_test
LIBRARY='run_8_dip_tet_arenosa'

# Change working directory
echo "Going to working directory ${SCRATCHDIR}"
cd "${SCRATCHDIR}"/ || exit 1
echo

# Required modules
echo "Loading modules"
module add trimmomatic-0.36 || exit 1
module add bbmap-38.42 || exit 1
module add fastQC-0.11.5 || exit 1
module add parallel-20200322 || exit 1
echo

# Copy data
echo "Copying..."
echo "Scripts etc. - /storage/praha1/home/${LOGNAME}/radseq/"
cp -a /storage/praha1/home/"${LOGNAME}"/radseq/{adaptors.fa,bin/rad_2_prep_2_run.sh} "${SCRATCHDIR}"/ || exit 1
echo "Data to process - ${DATADIR}/${LIBRARY}"
cp -a "${DATADIR}"/"${LIBRARY}"/1_demultiplexed "${SCRATCHDIR}"/ || exit 1
echo

# Running the task
echo "Preprocessing the FASTQ files..."
./rad_2_prep_2_run.sh -f 1_demultiplexed -c 4 -o 2_trimmed -d 2_dedup -q 3_qual_rep -a adaptors.fa -m 2 -t "${TRIMMOMATIC_BIN}" | tee radseq_prepare.log
echo

# Remove unneeded file
echo "Removing unneeded files"
rm -rf 1_demultiplexed adaptors.fa rad_2_prep_2_run.sh
echo

# Copy results back to storage
echo "Copying results back to ${DATADIR}"
cp -a "${SCRATCHDIR}" "${DATADIR}"/"${LIBRARY}"/ || export CLEAN_SCRATCH='false'
echo

exit

