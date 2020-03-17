#!/bin/bash

# qsub -l walltime=48:0:0 -l select=1:ncpus=2:mem=24gb:scratch_local=250gb -m abe bin/rad_5_hardfilter_invar_arabidopsis_cardamine.sh

# Clean-up of SCRATCH
trap 'clean_scratch' TERM EXIT
trap 'cp -a $SCRATCHDIR $DATADIR/ && clean_scratch' TERM

# Prepare the task
cp -a /storage/praha1/home/$LOGNAME/rad/5_hardfilter/* $SCRATCHDIR/ || exit 1
cp -a /storage/praha1/home/$LOGNAME/rad/ref/* $SCRATCHDIR/ || exit 1
# cp -a /storage/ostrava2-archive/tape_tape/backup/VO_cuni_prf_arab/shared/arabidopsis/joined_vcf/arenosa/*all*vcf.gz* $SCRATCHDIR/ || exit 1
cp -a /storage/ostrava2-archive/tape_tape/backup/VO_cuni_prf_arab/shared/arabidopsis/joined_vcf/lyrata/*all*vcf.gz* $SCRATCHDIR/ || exit 1

# Change working directory
cd $SCRATCHDIR/ || exit 2

# Launch it
module add gatk-3.8-0

# Create a filtering criterion
java -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R alygenomes.fasta -V lyrata_all.join.raw.vcf.gz -o lyrata_all.raw.vcf.gz --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR < 3.0" --filterName "Rad_filter1" || exit 1

# Select only non-variant sites
java -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V lyrata_all.raw.vcf.gz -selectType NO_VARIATION -o lyrata_all.raw.hardfilter.vcf.gz -L scaffold_1 -L scaffold_2 -L scaffold_3 -L scaffold_4 -L scaffold_5 -L scaffold_6 -L scaffold_7 -L scaffold_8 --excludeIntervals blacklisty.intervals --excludeFiltered || exit 1

# Set the filtering for required minimal DP
java -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R alygenomes.fasta -V lyrata_all.raw.hardfilter.vcf.gz -o lyrata_all.raw.hardfilter.dp50.vcf.gz --genotypeFilterExpression "DP < 50" --genotypeFilterName DP-50 --setFilteredGtToNocall || exit 1

# Set the filtering for required minimal average coverage for each called allele
java -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration -R alygenomes.fasta -V lyrata_all.raw.hardfilter.dp50.vcf.gz -o lyrata_all.raw.hardfilter.dp50.dpan5.vcf.gz --filterExpression "DP / AN < 5" --filterName DP-AN-5 || exit 1

# Select variants based on this interval list (NB variants with < defined coverage will be still present in VCF)
java -Xmx24g -jar $GATK/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V lyrata_all.raw.hardfilter.dp50.dpan5.vcf.gz -o lyrata_all.raw.hardfilter.dp50.dpan5.percmiss05.vcf.gz --maxNOCALLfraction 0.5 || exit 1

# Copy results back to storage
cp -a $SCRATCHDIR /storage/ostrava2-archive/tape_tape/backup/VO_cuni_prf_arab/shared/arabidopsis/filtered_vcf/ || export CLEAN_SCRATCH=false

# Clean-up of SCRATCH
if [ "$CLEAN_SCRATCH" != "false" ]; then rm -rf $SCRATCHDIR/*; fi

exit

