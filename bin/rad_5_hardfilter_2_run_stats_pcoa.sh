#!/bin/bash

# qsub -l walltime=48:0:0 -l select=1:ncpus=2:mem=48gb:scratch_local=50gb -m abe bin/rad_5_hardfilter_rstats_arabidopsis_cardamine.sh

# Change working directory
# cd /auto/pruhonice1-ibot/shared/brassicaceae/rad_vcf/joined_vcf/arenosa/ || exit 1
cd /auto/pruhonice1-ibot/shared/brassicaceae/rad_vcf/joined_vcf/lyrata/ || exit 1

# Copy R packages and script
cp -a /storage/praha1/home/gunnera/arabidopsis/5_hardfilter/rpkgs . || exit 1
cp /storage/praha1/home/gunnera/arabidopsis/5_hardfilter/stat_pca_rad.r . || exit 1

# Launch it
module add R-3.4.3-gcc || exit 1

# PCA and more statistics using R script

# Create directory for R outputs
mkdir rstats || exit 1

# Do the calculations
echo "Calculating statistics, PCAs and distances using R"
for VCFGZ in *.vcf.gz; do
	echo "Processing $VCFGZ"
	# Create output directory
	mkdir rstats/$VCFGZ || exit 1
	# Go to output directory
	cd rstats/$VCFGZ || exit 1
	# Copy R script to working directory
	cp ../../stat_pca_rad.r . || exit 1
	# Copy R packages
	cp -a ../../rpkgs .
	# Copy processed file
	cp ../../$VCFGZ . || exit 1
	# Prepare variable storing filename for R to read input tree
	export VCFR="$VCFGZ"
	# Do the calculations
	R CMD BATCH stat_pca_rad.r
	# Discard the variable
	unset VCFR
	# Cleanup
	rm -rf rpkgs stat_pca_rad.r $VCFGZ || exit 1
	# Go back
	cd ../..
	done

# Cleanup
echo
rm -rf rpkgs stat_pca_rad.r || exit 1

exit

