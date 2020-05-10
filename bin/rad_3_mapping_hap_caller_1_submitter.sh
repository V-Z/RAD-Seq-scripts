#!/bin/bash

# qsub -l walltime=168:0:0 -l select=1:ncpus=8:mem=128gb:scratch_local=400gb -q ibot -m abe ~/bin/rad_3_mapping_hap_caller_arabidopsis_cardamine.sh

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a $SCRATCHDIR $DATADIR/ && clean_scratch' TERM

# Location of data to be trimmed
DATADIR='/auto/pruhonice1-ibot/shared/brassicaceae/rad'

# Library to process
# run_1_2_dip_tet_harvard run_3_4_dip_tet_2017_01_arenosa_lyrata run_5_dip_tet_2017_11_lyrata run_6_7_dip_tet_2018_01_lyrata run_8_dip_tet_arenosa run_9_dip_tet_lyrata run_10_dip_tet_lyrata run_11_dip_tet_arenosa run_12_dip_tet_arenosa_lyrata
# cardamine_run_01_test cardamine_run_02_test
LIBRARY='cardamine_run_02_test'

# Reference
# alygenomes.fasta pseudohap_Camara_90M_10kb.fasta
REF='pseudohap_Camara_90M_10kb.fasta'

# Change working directory
cd "$SCRATCHDIR"/ || exit 1

# Launch it
module add parallel-20160622 || exit 1
module add bwa-0.7.3a || exit 1
module add samtools-1.9 || exit 1
module add picard-2.9.0 || exit 1
module add gatk-3.8-0 || exit 1
module rm openjdk-8 || exit 1
module add jdk-8 || exit 1

# Prepare the task
cp -a /storage/praha1/home/"$LOGNAME"/rad/3_mapping_hap_caller/* "$SCRATCHDIR"/ || exit 1
cp -a /storage/praha1/home/"$LOGNAME"/rad/ref/* "$SCRATCHDIR"/ || exit 1
cp -a "$DATADIR"/"$LIBRARY"/2_trimmed "$SCRATCHDIR"/ || exit 1

./radseq_3_mapping_hap_caller.sh -f 2_trimmed -c 8 -o mapped -a "$REF" -j /packages/run/jdk-8/current/bin/java -m 16000m -p /software/picard/2.9.0/build/libs/picard.jar -g "$GATK"/GenomeAnalysisTK.jar | tee mapping.log

# Copy results back to storage
cp -a "$SCRATCHDIR" "$DATADIR"/"$LIBRARY"/ || export CLEAN_SCRATCH='false'

exit

# Go to working directory
echo "Switching to ${DATADIR}"
cd "${DATADIR}"/ || exit 1
echo

# Processing all samples
echo "Processing all samples at $(date)..."
echo
while read -r SAMPLE; do
	echo "Processing ${SAMPLE}"
	qsub -l walltime=48:0:0 -l select=1:ncpus="${NCPU}":mem=8gb:scratch_local=10gb -q ibot -m abe -N HybPiper."${SAMPLE}" -v HYBPIPDIR="${HYBPIPDIR}",WORKDIR="${WORKDIR}",DATADIR="${DATADIR}",BAITFILE="${BAITFILE}",NCPU="${NCPU}",SAMPLE="${SAMPLE}" ~/hybseq/bin/hybseq_2_hybpiper_2_qsub.sh || exit 1
	echo
	done < "${SAMPLES}"

echo "All jobs submitted..."
echo

exit

