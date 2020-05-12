#!/bin/bash

# qsub -l walltime=48:0:0 -l select=1:ncpus=2:mem=16gb:scratch_local=50gb -q ibot -m abe -N RADSeq_mapping_hapcaller."${ALNB}" -v WORKDIR="${WORKDIR}",DATADIR="${DATADIR}",LIBRARY="${LIBRARY}",ALNF="${ALNB}",REF="${REF}" ~/hybseq/bin/rad_3_mapping_hap_caller_2_qsub.sh

# Checking if all required variables are provided
if [ -z "${ALNF}" ]; then
	echo "Error! Sample name not provided!"
	exit 1
	fi
if [ -z "${WORKDIR}" ]; then
	echo "Error! Data and scripts for HybSeq not provided!"
	exit 1
	fi
if [ -z "${DATADIR}" ]; then
	echo "Error! Directory with data to process not provided!"
	exit 1
	fi
if [ -z "${LIBRARY}" ]; then
	echo "Error! Library with data to process not provided!"
	exit 1
	fi
if [ -z "${REF}" ]; then
	echo "Error! Reference not provided!"
	exit 1
	fi

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a ${SCRATCHDIR} ${DATADIR}/ && clean_scratch' TERM

# Change working directory
echo "Going to working directory ${SCRATCHDIR}"
cd "${SCRATCHDIR}"/ || exit 1
echo

# Launch it
echo "Loading modules"
module add parallel-20200322 || exit 1
module add bwa-0.7.17 || exit 1
module add samtools-1.10 || exit 1
module add picard-2.22.1 || exit 1
module add jdk-8 || exit 1
echo

# Copy data
echo "Copying..."
echo "Scripts etc. - /storage/praha1/home/${LOGNAME}/radseq/"
cp -a /storage/praha1/home/"${LOGNAME}"/radseq/{"${REF%.*}*",bin/rad_3_mapping_hap_caller_3_run.sh} "${SCRATCHDIR}"/ || exit 1
echo "Data to process - ${DATADIR}/${LIBRARY}/2_trimmed/${ALNF}"
cp -a "${DATADIR}"/"${LIBRARY}"/2_trimmed/"${ALNF}"/* "${SCRATCHDIR}"/ || exit 1
echo "and ${DATADIR}/${LIBRARY}/3_dedup/${ALNF}"
cp -a "${DATADIR}"/"${LIBRARY}"/3_dedup/"${ALNF}"/* "${SCRATCHDIR}"/ || exit 1
echo

# Reference base name
echo "Obtaining basename of reference file ${ALNF}"
REFB="$(basename "${REF}")" || exit 1
echo

# Running the task
echo "Preprocessing the FASTQ files..."
./rad_3_mapping_hap_caller_3_run.sh -f "${ALNF}" -a "${REFB}" -j /packages/run/jdk-8/current/bin/java -m 16000m -p /software/picard/2.22.1/build/libs/picard.jar -g /auto/pruhonice1-ibot/home/gunnera/bin/GenomeAnalysisTK.jar | tee mapping_hap_caller."${ALNF}".log
echo

# Remove unneeded file
echo "Removing unneeded files"
rm "${REFB%.*}*" rad_3_mapping_hap_caller_3_run.sh "${ALNF}"*.f*q* || { export CLEAN_SCRATCH='false'; exit 1; }
echo

# Copy results back to storage
echo "Copying results back to ${DATADIR}"
cp -a "${SCRATCHDIR}"/* "${DATADIR}"/"${LIBRARY}"/mapped/ || { export CLEAN_SCRATCH='false'; exit 1; }
echo

exit

