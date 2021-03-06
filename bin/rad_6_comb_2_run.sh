#!/bin/bash

# Author: Vojtěch Zeisek, https://trapa.cz/
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# 

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# Initialize variables
VCFFILE='' # Input joint VCF file
VCFFILESNP='' # Output file: SNPs
VCFFILESNPHARD='' # Output file: Hardfiltering SNPs
VCFFILESNPPASS='' # Output file: Bi- and multiallelic SNPs
VCFFILESNPBIAL='' # Output file: Biallelic SNPs
VCFFILESNPFILT='' # Output file: Filtered non-variant SNPs
REF='' # Reference sequence
EXCLUDEPARALOGS='' # Paralogs to exclude
JAVA='' # PATH to custom Java binary
JAVAMEM='' # Memory limit for Picard and GATK
JAVAMEMTEST='^[0-9]+[kmgt]$' # Testing if provided value is a number with k, m, g or t
GATK='' # Path to directory containing GATK JAR file
OUTNAME='' # Base name of the output file
OUTNAMETEST='^[a-zA-Z0-9_.]+$' # Testing if base name of output file contains only valid characters
DECNUMTEST='^[0-9]+.?[0-9]*$' # Testing if provided value is decimal number
MAXFRACTFILTGENOT='' # Maximum fraction of samples filtered at the genotype level
NUMTEST='^[0-9]+$' # Testing if provided value is an integer
GENOTFILTDP='' # Genotype filtering level
MINDPANCOV='' # Minimal average coverage for each called allele
echo

# Parse initial arguments
while getopts "hrvf:n:a:e:j:m:g:l:w:y:" INITARGS; do
	case "${INITARGS}" in
		h) # Help and exit
			echo "Usage options:"
			echo -e "\t-h\tPrint this help and exit."
			echo -e "\t-r\tPrint references to used software and exit."
			echo -e "\t-v\tPrint script version and exit."
			echo -e "\t-f\tInput joint VCF file to be processed."
			echo -e "\t-n\tBase name of output files. Allowed characters are letters, numbers, underscore or dot. If not provided, output files will start with name of input file."
			echo -e "\t-a\tReference FASTA file."
			echo -e "\t-e\tParalogs to exclude."
			echo -e "\t-j\tOptional path to custom Java binary (default is output of \`command -v java\`; GATK requires Oracle Java)."
			echo -e "\t-m\tMaximal memory consumption allowed to GATK. Input as common for 'jar -Xmx', e.g. 24g for '-Xmx24g'. Default is 24g."
			echo -e "\t-g\tPath to GATK JAR file."
			echo -e "\t-l\tMaximum fraction of samples filtered at the genotype level. Provide value from 0 to 1. Default is 0.5."
			echo -e "\t-w\tMinimal required coverage. Provide an integer. Default is 4."
			echo -e "\t-y\tMinimal average coverage for each called allele. Provide an integer. Default is 4."
			echo -e "\tOutput files will be in same directory as input file (-f)."
			echo
			exit
			;;
		r) # References to cite and exit
			echo "Software to cite:"
			echo "* Bcftools, https://samtools.github.io/bcftools/"
			echo "* GATK, https://software.broadinstitute.org/gatk/"
			echo "* Java, https://java.com/"
			echo
			exit
			;;
		v) # Print script version and exit
			echo "Version: 1.0"
			echo
			exit
			;;
		f) # Input joint VCF file to be processed
			if [ -r "${OPTARG}" ]; then
				VCFFILE="${OPTARG}"
				echo "Input joint VCF file: ${VCFFILE}"
				echo
				else
					echo "Error! You did not provide path to input joint VCF file (-f) \"${OPTARG}\"!"
					echo
					exit 1
					fi
			;;
		n) # Base name of the output file
			if [[ ${OPTARG} =~ ${OUTNAMETEST} ]]; then
				OUTNAME="${OPTARG}"
				# Names of output files are derived from the base name provided by the user
				VCFFILESNP="${OUTNAME}".raw.snp.vcf.gz # SNPs
				VCFFILESNPHARD="${OUTNAME}".raw.hardfilter.snp.vcf.gz # Hardfiltering SNPs
				VCFFILESNPPASS="${OUTNAME}".raw.hardfilter.snp.pass.vcf.gz # Bi- and multiallelic SNPs
				VCFFILESNPBIAL="${OUTNAME}".raw.hardfilter.snp.pass.bial.vcf.gz # Biallelic SNPs
				echo "Name of the output files will start with: ${OUTNAME}"
				echo
				else
					echo "Error! You did not provide valid base name of output files (-n) \"${OPTARG}\"!"
					echo
					exit 1
					fi
			;;
		a) # Reference FASTA file
			if [ -r "${OPTARG}" ]; then
				REF="${OPTARG}"
				echo "Reference FASTA file: ${REF}"
				echo
				else
					echo "Error! You did not provide path to reference FASTA file (-a) \"${OPTARG}\"!"
					echo
					exit 1
					fi
			;;
		e) # Paralogs to exclude
			if [ -r "${OPTARG}" ]; then
				EXCLUDEPARALOGS="${OPTARG}"
				echo "Paralogs to exclude: ${EXCLUDEPARALOGS}"
				echo
				else
					echo "Error! You did not provide path to file with paralogs to exclude (-e) \"${OPTARG}\"!"
					echo
					exit 1
					fi
			;;
		j) # Path to custom Java binary
			if [ -x "${OPTARG}" ]; then
			JAVA="${OPTARG}"
			echo "Custom Java binary: ${JAVA}"
			echo
			else
				echo "Error! You did not provide path to custom Java binary (-j) \"${OPTARG}\"!"
				echo
				exit 1
				fi
			;;
		m) # Maximal Java memory consumption
			if [[ ${OPTARG} =~ ${JAVAMEMTEST} ]]; then
			JAVAMEM="${OPTARG}"
			echo "Maximal memory consumption by Java (GATK): ${JAVAMEM}"
			echo
			else
				echo "Error! You did not provide correct maximal memory consumption (-m), e.g. 8g or 6000m, \"${OPTARG}\"!"
				echo
				exit 1
				fi
			;;
		g) # Path to GATK JAR file
			if [ -r "${OPTARG}" ]; then
			GATK="${OPTARG}"
			echo "GATK JAR file: ${GATK}"
			echo
			else
				echo "Error! You did not provide path to GATK JAR file (-g) \"${OPTARG}\"!"
				echo
				exit 1
				fi
			;;
		l) # Maximum fraction of samples filtered at the genotype level
			if [[ ${OPTARG} =~ ${DECNUMTEST} ]]; then
			MAXFRACTFILTGENOT="${OPTARG}"
			echo "Maximum fraction of samples filtered at the genotype level: ${MAXFRACTFILTGENOT}"
			echo
			else
				echo "Error! You did not provide correct number for maximum fraction of samples filtered at the genotype level (-l) \"${OPTARG}\"!"
				echo
				exit 1
				fi
			;;
		w) # Minimal required coverage
			if [[ ${OPTARG} =~ ${NUMTEST} ]]; then
			GENOTFILTDP="${OPTARG}"
			echo "Minimal required coverage: ${GENOTFILTDP}"
			echo
			else
				echo "Error! You did not provide correct number for minimal required coverage (-w) \"${OPTARG}\"!"
				echo
				exit 1
				fi
			;;
		y) # Minimal average coverage for each called allele
		if [[ ${OPTARG} =~ ${NUMTEST} ]]; then
		MINDPANCOV="${OPTARG}"
		echo "Minimal average coverage for each called allele: ${GENOTFILTDP}"
		echo
		else
				echo "Error! You did not provide correct number for minimal average coverage for each called allele (-y) \"${OPTARG}\"!"
				echo
				exit 1
				fi
		;;
		*)
			echo "Error! Unknown option!"
			echo "See usage options: \"${0} -h\""
			echo
			exit 1
			;;
		esac
	done

# Check if all required binaries are available
function toolcheck {
	command -v "${1}" >/dev/null 2>&1 || {
		echo >&2 "Error! ${1} is required but not installed. Aborting. Please, install it."
		echo
		exit 1
		}
	}

toolcheck bcftools
toolcheck parallel

# Checking if all required parameters are provided
if [ -z "${VCFFILE}" ]; then
	echo "Error! Input joint VCF file (-f) was not specified!"
	echo "See usage options: \"${0} -h\""
	echo
	exit 1
	fi

if [ -z "${OUTNAME}" ]; then
	echo "Base name of output join VCF (-n) was not set. Names of output files will start with name of input file."
	echo
	# Names of output files are derived from the input name
	VCFFILESNP=$(echo "${VCFFILEB}" | sed 's/\(^.*\)vcf/\1snp.vcf/') # SNPs
	VCFFILESNPHARD=$(echo "${VCFFILESNP}" | sed 's/\(^.*\)snp/\1hardfilter.snp/') # Hardfiltering SNPs
	VCFFILESNPPASS=$(echo "${VCFFILESNPHARD}" | sed 's/\(^.*\)snp/\1snp.pass/') # Bi- and multiallelic SNPs
	VCFFILESNPBIAL=$(echo "${VCFFILESNPPASS}" | sed 's/\(^.*\)snp\.pass/\1snp.pass.bial/') # Biallelic SNPs
	VCFFILESNPFILT=$(echo "${VCFFILESNPPASS}" | sed 's/\(^.*\)snp\.pass/\1snp.pass.filt/') # Filtered non-variant SNPs
	fi

if [ -z "${REF}" ]; then
	echo "Error! Reference FASTA file (-a) was not specified!"
	echo "See usage options: \"${0} -h\""
	echo
	exit 1
	fi

if [ -z "${EXCLUDEPARALOGS}" ]; then
	echo "Error! File with paralogs to exclude (-e) was not specified!"
	echo "See usage options: \"${0} -h\""
	echo
	exit 1
	fi

if [ -z "${JAVA}" ]; then
	echo "Path to custom Java executable (-j) was not specified. Using default $(command -v java)"
	JAVA="$(command -v java)"
	toolcheck java
	echo
	fi

if [ -z "${JAVAMEM}" ]; then
	echo "Java memory consumption for GATK (-m) was not set. Using default value of 24g."
	JAVAMEM="24g"
	echo
	fi

if [ -z "${GATK}" ]; then
	echo "Error! Path to GATK JAR file (-g) was not specified!"
	echo "See usage options: \"${0} -h\""
	echo
	exit 1
	fi

if [ -z "${MAXFRACTFILTGENOT}" ]; then
	echo "Maximum fraction of samples filtered at the genotype level (-l) was not set. Using default value of 0.5."
	echo
	MAXFRACTFILTGENOT='0.5'
	fi

if [ -z "${GENOTFILTDP}" ]; then
	echo "Minimal required coverage (-w) was not set. Using default value of 4."
	echo
	GENOTFILTDP='4'
	fi

if [ -z "${MINDPANCOV}" ]; then
	echo "Minimal average coverage for each called allele (-y) was not set. Using default value of 4."
	echo
	MINDPANCOV='4'
	fi

# Exit on error
function operationfailed {
	echo "Error! Operation failed!"
	echo
	echo "See previous message(s) to be able to trace the problem."
	echo
	# Do not clean SCRATCHDIR, but copy content back to DATADIR
	export CLEAN_SCRATCH='false'
	exit 1
	}

echo "Start: $(date)"
echo

# Set the filtering for required minimal DP
echo "Marking filtered sites in the joint VCF ${VCFFILESNPBIAL}..."
echo
"${JAVA}" -Xmx"${JAVAMEM}" -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar "${GATK}" -T VariantFiltration -R "${REF}" -V "${VCFFILESNPBIAL}" -o "${VCFFILESNPBIAL%.vcf.gz}".dp"${GENOTFILTDP}".vcf.gz --genotypeFilterExpression "DP < ${GENOTFILTDP}" --genotypeFilterName DP-"${GENOTFILTDP}" --setFilteredGtToNocall || operationfailed
echo
echo "Marked filtered sites were saved as ${VCFFILESNPBIAL%.vcf.gz}.dp${GENOTFILTDP}.vcf.gz"
echo

# Set the filtering for required minimal average coverage for each called allele
echo "Checking if each allele that is called is covered by at least ${MINDPANCOV} reads on average"
echo
"${JAVA}" -Xmx"${JAVAMEM}" -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar "${GATK}" -T VariantFiltration -R "${REF}" -V "${VCFFILESNPBIAL%.vcf.gz}".dp"${GENOTFILTDP}".vcf.gz -o "${VCFFILESNPBIAL%.vcf.gz}".dp"${GENOTFILTDP}".dpan"${MINDPANCOV}".vcf.gz --filterExpression "DP / AN < ${MINDPANCOV}" --filterName DP-AN-"${MINDPANCOV}" || operationfailed
echo
echo "Marked filtered sites were saved as ${VCFFILESNPBIAL%.vcf.gz}.dp${GENOTFILTDP}.dpan${MINDPANCOV}.vcf.gz"
echo

# Select variants based on this interval list (NB variants with < defined coverage will be still present in VCF)
echo "Selecting variants based on presence in ${GENOTFILTDP} of indivs in ${VCFFILESNPBIAL%.vcf.gz}.dp${GENOTFILTDP}.vcf.gz joint VCF..."
echo
"${JAVA}" -Xmx"${JAVAMEM}" -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar "${GATK}" -T SelectVariants -R "${REF}" -V "${VCFFILESNPBIAL%.vcf.gz}".dp"${GENOTFILTDP}".dpan"${MINDPANCOV}".vcf.gz -o "${VCFFILESNPBIAL%.vcf.gz}".dp"${GENOTFILTDP}".dpan"${MINDPANCOV}".percmiss"${MAXFRACTFILTGENOT}".vcf.gz --maxNOCALLfraction "${MAXFRACTFILTGENOT}" || operationfailed
echo
echo "Final selected variants were saved as ${VCFFILESNPBIAL%.vcf.gz}.dp${GENOTFILTDP}.dpan${MINDPANCOV}.percmiss${MAXFRACTFILTGENOT}.vcf.gz"
echo

echo "All output files are in directory \"$(pwd)\" and names start with \"${OUTNAME}*\""
echo

# Statistics using BCFtools
echo "Statistics of SNPs in VCF files using bcftools"
find . -name "*.vcf.gz" | parallel -j 2 "echo '{/}' && bcftools stats -F ${REF} '{}' > '{.}'.stats.txt || operationfailed"
echo

# Statistics using R script
echo "Calculating statistics, PCAs and distances using R"
echo
for VCFGZ in *.vcf.gz; do
	echo "Processing ${VCFGZ}"
	echo
	# Create output directory
	mkdir "${VCFGZ%.vcf.gz}" || operationfailed
	# Go to output directory
	cd "${VCFGZ%.vcf.gz}" || operationfailed
	# Copy R script to working directory, R packages, processed file
	cp -a ../{rad_5_hardfilter_2_stats.r,rpackages,"${VCFGZ}","${VCFGZ}".tbi} . || operationfailed
	# Prepare variable storing filename for R to read input tree
	export VCFR="${VCFGZ}" || operationfailed
	# Do the calculations
	R CMD BATCH --no-save --no-restore rad_5_hardfilter_2_stats.r "${VCFGZ%.vcf.gz}".log
	# Discard the variable
	unset VCFR || operationfailed
	# Cleanup
	rm -rf rad_5_hardfilter_2_stats.r rpackages "${VCFGZ}" "${VCFGZ}".tbi || operationfailed
	# Go back
	cd ../ || operationfailed
	echo
	done

echo "End: $(date)"
echo

exit



cd $SCRATCHDIR
mkdir tmp
cp -a ~/radseq/ref/arabidopsis/* .
module add jdk-8
cp /storage/pruhonice1-ibot/shared/brassicaceae/rad_vcf/2_filtered_vcf/arenosa/{arenosa_var.join.vcf.filtered.snp.hardfilter.pass.bial.dp4.dpan4.percmiss0.7.vcf*,arenosa_all.join.vcf.filtered.invar.qual.pass.vcf*} .



java -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar /storage/pruhonice1-ibot/home/gunnera/bin/GenomeAnalysisTK.jar -T CombineVariants -R alygenomes.fasta -V arenosa_all.join.vcf.filtered.invar.qual.pass.vcf.gz -V arenosa_var.join.vcf.filtered.snp.hardfilter.pass.bial.dp4.dpan4.percmiss0.7.vcf.gz -o arenosa_comb.vcf.gz -genotypeMergeOptions UNSORTED -filteredRecordsMergeType KEEP_UNCONDITIONAL | tee -a arenosa_comb.log

java -Xmx19g -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar /storage/pruhonice1-ibot/home/gunnera/bin/GenomeAnalysisTK.jar -T VariantFiltration -R alygenomes.fasta -V arenosa_comb.vcf.gz -o arenosa_comb.dp4.vcf.gz --genotypeFilterExpression "DP < 4" --genotypeFilterName DP-4 --setFilteredGtToNocall | tee -a arenosa_comb.log

java -Xmx19g -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar /storage/pruhonice1-ibot/home/gunnera/bin/GenomeAnalysisTK.jar -T VariantFiltration -R alygenomes.fasta -V arenosa_comb.dp4.vcf.gz -o arenosa_comb.dp4.dpan4.vcf.gz --filterExpression "DP / AN < 4" --filterName DP-AN-4 | tee -a arenosa_comb.log

java -Xmx19g -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar /storage/pruhonice1-ibot/home/gunnera/bin/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V arenosa_comb.dp4.dpan4.vcf.gz -o arenosa_comb.dp4.dpan4.percmiss0.7.vcf.gz --maxNOCALLfraction 0.7 | tee -a arenosa_comb.log



java -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar /storage/pruhonice1-ibot/home/gunnera/bin/GenomeAnalysisTK.jar -T CombineVariants -R alygenomes.fasta -V arenosa_600_U_rad.merged.vcf.gz -V arenosa_comb.dp4.dpan4.percmiss0.7.vcf.gz -o arenosa_comb_600.vcf.gz -genotypeMergeOptions UNIQUIFY

java -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar /storage/pruhonice1-ibot/home/gunnera/bin/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V arenosa_comb_600.vcf.gz -o arenosa_comb_600.vcf.gz -selectType SNP -L scaffold_1 -L scaffold_2 -L scaffold_3 -L scaffold_4 -L scaffold_5 -L scaffold_6 -L scaffold_7 -L scaffold_8 --excludeIntervals blacklisty.intervals

java -Xmx19g -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar /storage/pruhonice1-ibot/home/gunnera/bin/GenomeAnalysisTK.jar -T VariantFiltration -R alygenomes.fasta -V arenosa_comb_600.vcf.gz -o arenosa_comb_600.dp8.vcf.gz --genotypeFilterExpression "DP < 8" --genotypeFilterName DP-8 --setFilteredGtToNocall

java -Xmx19g -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar /storage/pruhonice1-ibot/home/gunnera/bin/GenomeAnalysisTK.jar -T VariantFiltration -R alygenomes.fasta -V arenosa_comb_600.dp8.vcf.gz -o arenosa_comb_600.dp8.dpan8.vcf.gz --filterExpression "DP / AN < 8" --filterName DP-AN-8

java -Xmx19g -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar /storage/pruhonice1-ibot/home/gunnera/bin/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V arenosa_comb_600.dp8.dpan8.vcf.gz -o arenosa_comb_600.dp8.dpan8.percmiss05.vcf.gz --maxNOCALLfraction 0.5

java -Xmx19g -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar /storage/pruhonice1-ibot/home/gunnera/bin/GenomeAnalysisTK.jar -T SelectVariants -R alygenomes.fasta -V arenosa_comb_600.dp8.dpan8.vcf.gz -o arenosa_comb_600.dp8.dpan8.percmiss08.vcf.gz --maxNOCALLfraction 0.8



cp -a $SCRATCHDIR /storage/pruhonice1-ibot/shared/brassicaceae/rad_vcf/3_combined_vcf/


