#!/bin/bash

# Author: Vojtěch Zeisek, https://trapa.cz/
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# 

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# qsub -l walltime=48:0:0 -l select=1:ncpus=8:mem=64gb:scratch_local=100gb -q ibot -m abe ~/radseq/bin/rad_4_genotype_vcf_1_qsub.sh

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a "${SCRATCHDIR}" "${DATADIR}"/ && clean_scratch' TERM

# Location of data to merge
DATADIR='/storage/pruhonice1-ibot/shared/brassicaceae/rad_vcf/0_for_join/arenosa'
# DATADIR='/storage/pruhonice1-ibot/shared/brassicaceae/rad_vcf/0_for_join/lyrata'

# Reference
# ref/arabidopsis/alygenomes.fasta ref/cardamine/pseudohap_Camara_90M_10kb.fasta ref/minuartia/minuartia_combined_ref.fasta
REF='ref/arabidopsis/alygenomes.fasta'

# Change working directory
echo "Going to working directory ${SCRATCHDIR}"
cd "${SCRATCHDIR}"/ || exit 1
echo

# Required modules
echo "Loading modules"
module add jdk-8 || exit 1
module add parallel-20200322 || exit 1
echo

# Copy data
echo "Copying..."
echo "Scripts etc. - /storage/pruhonice1-ibot/home/${LOGNAME}/radseq/"
cp -a /storage/pruhonice1-ibot/home/"${LOGNAME}"/radseq/{${REF%.*}*,bin/rad_4_genotype_vcf_2_run.sh} "${SCRATCHDIR}"/ || exit 1
echo "Data to process - ${DATADIR}"
cp -a "${DATADIR}" "${SCRATCHDIR}"/ || exit 1
echo

# Reference base name
echo "Obtaining basename of reference file ${REF}"
REFB="$(basename "${REF}")" || exit 1
echo

# Data dir base name
echo "Obtaining basename of data directory ${DATADIR}"
DATADIRB="$(basename "${DATADIR}")" || exit 1
echo

# Temp directory
echo "Creating temporal directory"
mkdir tmp || exit 1
echo

# Running the task
echo "Preprocessing the FASTQ files..."
# ./rad_4_genotype_vcf_2_run.sh -w "raw.g.vcf" -u ".gz" -x ".join.raw.vcf.gz" -f "${DATADIRB}" -c 7 -o "${DATADIRB}"_vcf -n "${DATADIRB}"_var -a "${REFB}" -j /packages/run/jdk-8/current/bin/java -m 63g -g /storage/pruhonice1-ibot/home/"${LOGNAME}"/bin/GenomeAnalysisTK.jar | tee "${DATADIRB}"_var_joining_genotype_vcf.log
./rad_4_genotype_vcf_2_run.sh -w "raw.g.vcf" -u ".gz" -x ".join.raw.vcf.gz" -f "${DATADIRB}" -c 7 -o "${DATADIRB}"_vcf -n "${DATADIRB}"_all -a "${REFB}" -j /packages/run/jdk-8/current/bin/java -m 63g -g /storage/pruhonice1-ibot/home/"${LOGNAME}"/bin/GenomeAnalysisTK.jar -i | tee "${DATADIRB}"_all_joining_genotype_vcf.log
echo

# Remove unneeded file
echo "Removing unneeded files"
rm -rf tmp ${REFB%.*}* "${DATADIRB}" || { export CLEAN_SCRATCH='false'; exit 1; }
echo

# Copy results back to storage
echo "Copying results back to ${DATADIR}"
cp -a "${SCRATCHDIR}" "${DATADIR}"/ || { export CLEAN_SCRATCH='false'; exit 1; }
echo

exit

