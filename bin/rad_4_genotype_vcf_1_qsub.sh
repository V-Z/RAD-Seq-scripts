#!/bin/bash

# qsub -l walltime=48:0:0 -l select=1:ncpus=6:mem=16gb:scratch_local=100gb -m abe bin/arabidopsis_4_genotype_vcf.sh

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a $SCRATCHDIR $DATADIR/ && clean_scratch' TERM

# Prepare the task
cp -a /storage/praha1/home/$LOGNAME/rad/4_genotype_vcf/* $SCRATCHDIR/ || exit 1
cp -a /storage/praha1/home/$LOGNAME/rad/ref/* $SCRATCHDIR/ || exit 1
# cp -a /auto/pruhonice1-ibot/shared/brassicaceae/rad_vcf/for_join/arenosa $SCRATCHDIR/ || exit 1
cp -a /auto/pruhonice1-ibot/shared/brassicaceae/rad_vcf/for_join/lyrata $SCRATCHDIR/ || exit 1

# Change working directory
cd $SCRATCHDIR/ || exit 2

# Launch it
module add parallel-20160622
module add gatk-3.8-0

# ./radseq_4_genotype_vcf.sh -w "raw.g.vcf" -u ".gz" -x ".join.raw.vcf.gz" -f arenosa -c 5 -o joined_vcf -n alpine_var -a alygenomes.fasta -j /packages/run/jdk-8/current/bin/java -m 15500m -g $GATK/GenomeAnalysisTK.jar | tee alpine_var_joining_genotype_vcf.log
# ./radseq_4_genotype_vcf.sh -w "raw.g.vcf" -u ".gz" -x ".join.raw.vcf.gz" -f arenosa -c 5 -o joined_vcf -n alpine_all -a alygenomes.fasta -j /packages/run/jdk-8/current/bin/java -m 15500m -g $GATK/GenomeAnalysisTK.jar -i | tee alpine_all_joining_genotype_vcf.log
# ./radseq_4_genotype_vcf.sh -w "raw.g.vcf" -u ".gz" -x ".join.raw.vcf.gz" -f lyrata -c 5 -o joined_vcf -n lyrata_var -a alygenomes.fasta -j /packages/run/jdk-8/current/bin/java -m 15500m -g $GATK/GenomeAnalysisTK.jar | tee lyrata_var_joining_genotype_vcf.log
./radseq_4_genotype_vcf.sh -w "raw.g.vcf" -u ".gz" -x ".join.raw.vcf.gz" -f lyrata -c 5 -o joined_vcf -n lyrata_all -a alygenomes.fasta -j /packages/run/jdk-8/current/bin/java -m 15500m -g $GATK/GenomeAnalysisTK.jar -i | tee lyrata_all_joining_genotype_vcf.log

# Copy results back to storage
cp -a $SCRATCHDIR /auto/pruhonice1-ibot/shared/brassicaceae/rad_vcf/joined_vcf/ || export CLEAN_SCRATCH=false

# Clean-up of SCRATCH
if [ "$CLEAN_SCRATCH" != "false" ]; then rm -rf $SCRATCHDIR/*; fi

exit

