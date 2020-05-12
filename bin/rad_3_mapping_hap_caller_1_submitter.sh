#!/bin/bash

# Set data directories
WORKDIR="/storage/praha1/home/${LOGNAME}/radseq"

# Location of data to be trimmed
DATADIR='/auto/pruhonice1-ibot/shared/brassicaceae/rad'
# Library to process
# run_01_02_brian run_01_02_dip_tet_harvard run_03_04_dip_tet_2017_01_arenosa_lyrata run_05_dip_tet_2017_11_lyrata run_06_07_dip_tet_2018_01_lyrata run_2018_06 run_08_dip_tet_2018_08_arenosa run_09_dip_tet_2019_01_lyrata run_10_dip_tet_2019_05_lyrata run_11_dip_tet_2019_06_arenosa run_12_dip_tet_2019_12_arenosa_lyrata
# cardamine_run_01_test cardamine_run_02_test
LIBRARY='run_08_dip_tet_2018_08_arenosa'

# Reference
# ref/arabidopsis/alygenomes.fasta ref/cardamine/pseudohap_Camara_90M_10kb.fasta
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
	qsub -l walltime=48:0:0 -l select=1:ncpus=2:mem=16gb:scratch_local=50gb -q ibot -m abe -N RADSeq_mapping_hapcaller."${ALNB}" -v WORKDIR="${WORKDIR}",DATADIR="${DATADIR}",LIBRARY="${LIBRARY}",ALNF="${ALNB}",REF="${REF}" ~/radseq/bin/rad_3_mapping_hap_caller_2_qsub.sh || exit 1
	echo
	done

echo "All jobs submitted..."
echo

exit

