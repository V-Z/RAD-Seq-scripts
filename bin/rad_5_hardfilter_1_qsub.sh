#!/bin/bash

# Author: Vojtěch Zeisek, https://trapa.cz/
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# 

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# qsub -l walltime=96:0:0 -l select=1:ncpus=2:mem=96gb:scratch_local=100gb -q ibot -m abe ~/radseq/bin/rad_5_hardfilter_1_qsub.sh

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a "${SCRATCHDIR}" "${DATADIR}"/ && clean_scratch' TERM

# Location of data to filter
DATADIR='/storage/pruhonice1-ibot/shared/brassicaceae/rad_vcf/1_joined_vcf/arenosa'
# DATADIR='/storage/pruhonice1-ibot/shared/brassicaceae/rad_vcf/1_joined_vcf/lyrata'
# DATADIR='/storage/pruhonice1-ibot/shared/brassicaceae/rad_vcf/1_joined_vcf/minuartia'

# Sample to process
# arenosa_all.join.vcf.gz arenosa_var.join.vcf.gz lyrata_all.join.vcf.gz lyrata_var.join.vcf.gz minuartia_all.join.vcf.gz minuartia_var.join.vcf.gz
SAMPLE='arenosa_var.join.vcf.gz'

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
module add bcftools-1.10.2 || exit 1
module add R-4.0.0-gcc || exit 1
echo

# Copy data
echo "Copying..."
echo "Scripts etc. - /storage/pruhonice1-ibot/home/${LOGNAME}/radseq/"
cp -a /storage/pruhonice1-ibot/home/"${LOGNAME}"/radseq/{blacklisty.intervals,${REF%.*}*,rpackages,bin/rad_5_hardfilter_2_run.sh,bin/rad_5_hardfilter_2_stats.r} "${SCRATCHDIR}"/ || exit 1
echo "Data to process - ${DATADIR}"
cp "${DATADIR}"/"${SAMPLE}" "${DATADIR}"/"${SAMPLE}".tbi "${SCRATCHDIR}"/ || exit 1
echo

# Reference base name
echo "Obtaining basename of reference file ${REF}"
REFB="$(basename "${REF}")" || exit 1
echo

# Temp directory
echo "Creating temporal directory"
mkdir tmp || exit 1
echo

# Running the task
echo "Preprocessing the FASTQ files..."
./rad_5_hardfilter_2_run.sh -f "${SAMPLE}" -n ${SAMPLE%.*}.filtered -a "${REFB}" -m 95g -g /storage/pruhonice1-ibot/home/"${LOGNAME}"/bin/GenomeAnalysisTK.jar -l 0.5 -w 4 -y 4 | tee ${SAMPLE%.*}_hardfilter.log
# ./rad_5_hardfilter_2_run.sh -f "${SAMPLE}" -n ${SAMPLE%.*}.filtered -a "${REFB}" -m 95g -g /storage/pruhonice1-ibot/home/"${LOGNAME}"/bin/GenomeAnalysisTK.jar -l 0.5 -w 4 -y 4 -i | tee ${SAMPLE%.*}_hardfilter.log
echo

# Remove unneeded file
echo "Removing unneeded files"
rm -rf tmp ${REFB%.*}* blacklisty.intervals "${SAMPLE}" "${SAMPLE}".tbi rad_5_hardfilter_2_run.sh rad_5_hardfilter_2_stats.r rpackages || { export CLEAN_SCRATCH='false'; exit 1; }
echo

# Copy results back to storage
echo "Copying results back to ${DATADIR}"
cp -a "${SCRATCHDIR}" "${DATADIR}"/ || { export CLEAN_SCRATCH='false'; exit 1; }
echo

exit

