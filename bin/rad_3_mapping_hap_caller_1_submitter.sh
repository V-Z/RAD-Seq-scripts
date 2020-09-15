#!/bin/bash

# Author: Vojtěch Zeisek, https://trapa.cz/
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# 

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# Set data directories
WORKDIR="/storage/pruhonice1-ibot/home/${LOGNAME}/radseq"

# Location of data to map and call variants
DATADIR='/storage/pruhonice1-ibot/shared/brassicaceae/rad'
# Library to process
# run_01_02_brian run_01_02_dip_tet_harvard run_03_04_dip_tet_2017_01_arenosa_lyrata run_05_dip_tet_2017_11_lyrata run_06_07_dip_tet_2018_01_lyrata run_2018_06_tet_arenosa run_08_dip_tet_2018_08_arenosa run_09_dip_tet_2019_01_lyrata run_10_dip_tet_2019_05_lyrata run_10_dip_tet_2019_05_lyrata_merged run_11_dip_tet_2019_06_arenosa run_12_dip_tet_2019_12_arenosa_lyrata run_13_dip_tet_2020_08_arenosa
# cardamine_run_01_test cardamine_run_02_test
# minuartia_run_01_2017 minuartia_run_02_2019 minuartia_run_02_2019_merged
LIBRARY='run_2018_06_tet_arenosa'

# Reference
# ref/arabidopsis/alygenomes.fasta ref/cardamine/pseudohap_Camara_90M_10kb.fasta ref/minuartia/minuartia_combined_ref.fasta
REF='ref/arabidopsis/alygenomes.fasta'

# Submitting individual tasks

# Go to working directory
echo "Switching to ${DATADIR}/${LIBRARY}"
cd "${DATADIR}"/"${LIBRARY}"/ || exit 1
echo

# Make output directory
echo "Making output directory"
mkdir mapped
echo

# Processing all samples
echo "Processing all samples at $(date)..."
echo
for ALN in $(find 3_dedup/ -type d | tail -n+2 | sort); do
	ALNB="$(basename "${ALN}")"
	echo "Processing ${ALNB}"
	qsub -l walltime=48:0:0 -l select=1:ncpus=2:mem=16gb:scratch_local=10gb -q ibot -N RADSeq_mapping_hapcaller."${ALNB}" -v WORKDIR="${WORKDIR}",DATADIR="${DATADIR}",LIBRARY="${LIBRARY}",REF="${REF}",ALNF="${ALNB}" ~/radseq/bin/rad_3_mapping_hap_caller_2_qsub.sh || exit 1
	echo
	done

echo "All jobs submitted..."
echo

exit

