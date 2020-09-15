#!/bin/bash

# Author: Vojtěch Zeisek, https://trapa.cz/
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# 

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# qsub -l walltime=24:0:0 -l select=1:ncpus=1:mem=48gb:scratch_local=100gb -q ibot -m abe ~/radseq/bin/rad_5_hardfilter_2_run_stats_pcoa.sh

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a "${SCRATCHDIR}" "${DATADIR}"/ && clean_scratch' TERM

# Data location
DATADIR='/storage/pruhonice1-ibot/shared/brassicaceae/rad_vcf/2_filtered_vcf'

# Reference
# ref/arabidopsis/alygenomes.fasta ref/cardamine/pseudohap_Camara_90M_10kb.fasta ref/minuartia/minuartia_combined_ref.fasta
REF='ref/arabidopsis/alygenomes.fasta'

# Change working directory
echo "Going to working directory ${SCRATCHDIR}"
cd "${SCRATCHDIR}"/ || exit 1
echo

# Required modules
echo "Loading modules"
module add bcftools-1.10.2 || exit 1
module add R-4.0.0-gcc || exit 1
echo

# Copy data
echo "Copying..."
echo "Scripts etc. - /storage/pruhonice1-ibot/home/${LOGNAME}/radseq/"
cp -a /storage/pruhonice1-ibot/home/"${LOGNAME}"/radseq/{${REF%.*}*,rpackages,bin/rad_5_hardfilter_2_stats.r} "${SCRATCHDIR}"/ || exit 1
echo "Data to process"
cp "${DATADIR}"/arenosa/* "${SCRATCHDIR}"/ || exit 1
# cp "${DATADIR}"/lyrata/* "${SCRATCHDIR}"/ || exit 1
echo

# Reference base name
echo "Obtaining basename of reference file ${REF}"
REFB="$(basename "${REF}")" || exit 1
echo

# Statistics using BCFtools

echo "Calculating statistics using BCFtools"
echo
for VCFGZ in *.vcf.gz; do
	echo "Processing ${VCFGZ}"
	bcftools stats -F "${REFB}" "${VCFGZ}" > "${VCFGZ%.vcf.gz}".stats.txt || { export CLEAN_SCRATCH='false'; exit 1; }
	echo
	done

# Statistics using R script

# Do the calculations
echo "Calculating statistics, PCAs and distances using R"
echo
for VCFGZ in *.vcf.gz; do
	echo "Processing ${VCFGZ}"
	echo
	# Create output directory
	mkdir "${VCFGZ%.vcf.gz}" || exit 1
	# Go to output directory
	cd "${VCFGZ%.vcf.gz}" || exit 1
	# Copy R script to working directory, R packages, processed file
	cp -a ../{rad_5_hardfilter_2_stats.r,rpackages,"${VCFGZ}","${VCFGZ}".tbi} . || exit 1
	# Prepare variable storing filename for R to read input tree
	export VCFR="${VCFGZ}" || exit 1
	# Do the calculations
	R CMD BATCH --no-save --no-restore rad_5_hardfilter_2_stats.r "${VCFGZ%.vcf.gz}".log
	# Discard the variable
	unset VCFR || { export CLEAN_SCRATCH='false'; exit 1; }
	# Cleanup
	rm -rf rad_5_hardfilter_2_stats.r rpackages "${VCFGZ}" "${VCFGZ}".tbi || { export CLEAN_SCRATCH='false'; exit 1; }
	# Go back
	cd ../ || { export CLEAN_SCRATCH='false'; exit 1; }
	echo
	done

# Remove unneeded file
echo "Removing unneeded files"
rm -rf rad_5_hardfilter_2_stats.r ${REFB%.*}* rpackages ./*.vcf.gz ./*.vcf.gz.tbi || { export CLEAN_SCRATCH='false'; exit 1; }
echo

# Copy results back to storage
echo "Copying results back to ${DATADIR}"
cp -a "${SCRATCHDIR}" "${DATADIR}"/ || export CLEAN_SCRATCH='false'
echo

exit

