#!/bin/bash

# qsub -l walltime=24:0:0 -l select=1:ncpus=2:mem=8gb:scratch_local=50gb -m abe bin/arabidopsis_5_hardfilter.sh

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a $SCRATCHDIR $DATADIR/ && clean_scratch' TERM

# Prepare the task
cp -a /storage/praha1/home/$LOGNAME/rad/5_hardfilter/* $SCRATCHDIR/ || exit 1
cp -a /storage/praha1/home/$LOGNAME/rad/ref/* $SCRATCHDIR/ || exit 1
# cp -a /storage/ostrava2-archive/tape_tape/backup/VO_cuni_prf_arab/shared/rad_vcf/joined_vcf/arenosa $SCRATCHDIR/ || exit 1
cp -a /storage/ostrava2-archive/tape_tape/backup/VO_cuni_prf_arab/shared/rad_vcf/joined_vcf/lyrata $SCRATCHDIR/ || exit 1

# Change working directory
cd $SCRATCHDIR/ || exit 2

# Launch it
module add gatk-3.8-0
module add parallel-20160622
module add bcftools-1.9
module add R-3.4.3-gcc

# ./radseq_5_hardfilter.sh -f arenosa/alpine_var.join.raw.vcf.gz -n alpine_var_filtered -a alygenomes.fasta -e blacklisty.intervals -m 7800m -g $GATK/GenomeAnalysisTK.jar -l 0.25 -w 4 -y 4 | tee alpine_var_hardfilter.log
# ./radseq_5_hardfilter.sh -f arenosa/alpine_var.join.raw.vcf.gz -n alpine_var_filtered -a alygenomes.fasta -e blacklisty.intervals -m 7800m -g $GATK/GenomeAnalysisTK.jar -l 0.5 -w 4 -y 4 | tee -a alpine_var_hardfilter.log
# ./radseq_5_hardfilter.sh -f lyrata/lyrata_var.join.raw.vcf.gz -n lyrata_var_filtered -a alygenomes.fasta -e blacklisty.intervals -m 7800m -g $GATK/GenomeAnalysisTK.jar -l 0.25 -w 4 -y 4 | tee lyrata_var_hardfilter.log
./radseq_5_hardfilter.sh -f lyrata/lyrata_var.join.raw.vcf.gz -n lyrata_var_filtered -a alygenomes.fasta -e blacklisty.intervals -m 7800m -g $GATK/GenomeAnalysisTK.jar -l 0.5 -w 4 -y 4 | tee -a lyrata_var_hardfilter.log

# Copy results back to storage
cp -a $SCRATCHDIR /storage/ostrava2-archive/tape_tape/backup/VO_cuni_prf_arab/shared/rad_vcf/filtered_vcf/ || export CLEAN_SCRATCH=false

# Clean-up of SCRATCH
if [ "$CLEAN_SCRATCH" != "false" ]; then rm -rf $SCRATCHDIR/*; fi

exit

