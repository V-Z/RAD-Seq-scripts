#!/bin/bash

# Initialize variables
SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) # Determine script's directory
VCFFILE='' # Input joint VCF file
VCFFILEB='' # Only filename of the joint VCF file
VCFDIR='' # Directory with input joint VCF file
VCFFILESNP='' # Output file: SNPs
VCFFILESNPHARD='' # Output file: Hardfiltering SNPs
VCFFILESNPPASS='' # Output file: Bi- and multiallelic SNPs
VCFFILESNPBIAL='' # Output file: Biallelic SNPs
STARTDIR=$(pwd) # Current working directory
REF='' # Reference sequence
REFB='' # Only filename of the reference sequence
EXCLUDEPARALOGS='' # Paralogs to exclude
EXCLUDEPARALOGSB='' # Only file name of paralogs to exclude
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
echo

# Parse initial arguments
while getopts "hrvf:n:a:e:j:m:g:l:w:" INITARGS; do
	case "$INITARGS" in
		h) # Help and exit
			echo "Usage options:"
			echo -e "\t-h\tPrint this help and exit."
			echo -e "\t-r\tPrint references to used software and exit."
			echo -e "\t-v\tPrint script version and exit."
			echo -e "\t-f\tInput joint VCF file to be processed."
			echo -e "\t-n\tBase name of output files. Allowed characters are letters, numbers, underscore or dot. If not provided, output files will start with name of input file."
			echo -e "\t-a\tReference FASTA file."
			echo -e "\t-e\tParalogs to exclude."
			echo -e "\t-j\tOptional path to custom Java binary (default is output of \`which java\`; GATK requires Oracle Java)."
			echo -e "\t-m\tMaximal memory consumption allowed to GATK. Input as common for 'jar -Xmx', e.g. 12g for '-Xmx12g'. Default is 12g."
			echo -e "\t-g\tPath to GATK JAR file."
			echo -e "\t-l\tMaximum fraction of samples filtered at the genotype level. Provide value from 0 to 1. Default is 0.1."
			echo -e "\t-w\tMinimal required coverage. Provide an integer. Default is 8."
			echo -e "\tOutput files will be in same directory as input file (-f)."
			echo
			exit
			;;
		r) # References to cite and exit
			echo "Software to cite:"
			echo "* Bcftools, https://samtools.github.io/bcftools/"
			echo "* GATK, https://software.broadinstitute.org/gatk/"
			echo "* Java, https://java.com/ or http://openjdk.java.net/"
			echo "* R, https://www.r-project.org/"
			echo
			exit
			;;
		v) # Print script version and exit
			echo "Version: 1.0"
			echo
			exit
			;;
		f) # Input joint VCF file to be processed
			if [ -r "$OPTARG" ]; then
				VCFFILE="$OPTARG"
				VCFFILEB=$(basename $VCFFILE) # Only filename of the joint VCF file
				VCFDIR=$(dirname $VCFFILE) # Only directory name of the VCF file
				if [ ! -w "$VCFDIR" ]; then
					echo "Error! Directory $VCFDIR containing $VCFFILEB is not writable!"
					echo
					exit 1
					fi
				echo "Input joint VCF file: $VCFFILE"
				echo
				else
					echo "Error! You did not provide path to input joint VCF file (-f) \"$OPTARG\"!"
					echo
					exit 1
					fi
			;;
		n) # Base name of the output file
			if [[ $OPTARG =~ $OUTNAMETEST ]]; then
				OUTNAME="$OPTARG"
				# Names of output files are derived from the base name provided by the user
				VCFFILESNP=$OUTNAME.raw.snp.vcf.gz # SNPs
				VCFFILESNPHARD=$OUTNAME.raw.hardfilter.snp.vcf.gz # Hardfiltering SNPs
				VCFFILESNPPASS=$OUTNAME.raw.hardfilter.snp.pass.vcf.gz # Bi- and multiallelic SNPs
				VCFFILESNPBIAL=$OUTNAME.raw.hardfilter.snp.pass.bial.vcf.gz # Biallelic SNPs
				echo "Name of the output files will start with: $OUTNAME"
				echo
				else
					echo "Error! You did not provide valid base name of output files (-n) \"$OPTARG\"!"
					echo
					exit 1
					fi
			;;
		a) # Reference FASTA file
			if [ -r "$OPTARG" ]; then
				REF="$OPTARG"
				REFB=$(basename $REF) # Only filename of the reference sequence
				echo "Reference FASTA file: $REF"
				echo
				else
					echo "Error! You did not provide path to reference FASTA file (-a) \"$OPTARG\"!"
					echo
					exit 1
					fi
			;;
		e) # Paralogs to exclude
			if [ -r "$OPTARG" ]; then
				EXCLUDEPARALOGS="$OPTARG"
				EXCLUDEPARALOGSB=$(basename $EXCLUDEPARALOGS) # Only filename of the reference sequence
				echo "Paralogs to exclude: $EXCLUDEPARALOGS"
				echo
				else
					echo "Error! You did not provide path to file with paralogs to exclude (-e) \"$OPTARG\"!"
					echo
					exit 1
					fi
			;;
		j) # Path to custom Java binary
			if [ -x "$OPTARG" ]; then
			JAVA="$OPTARG"
			echo "Custom Java binary: $JAVA"
			echo
			else
				echo "Error! You did not provide path to custom Java binary (-j) \"$OPTARG\"!"
				echo
				exit 1
				fi
			;;
		m) # Maximal Java memory consumption
			if [[ $OPTARG =~ $JAVAMEMTEST ]]; then
			JAVAMEM="$OPTARG"
			echo "Maximal memory consumption by Java (GATK): $JAVAMEM"
			echo
			else
				echo "Error! You did not provide correct maximal memory consumption (-m), e.g. 8g or 6000m, \"$OPTARG\"!"
				echo
				exit 1
				fi
			;;
		g) # Path to GATK JAR file
			if [ -r "$OPTARG" ]; then
			GATK="$OPTARG"
			echo "GATK JAR file: $GATK"
			echo
			else
				echo "Error! You did not provide path to GATK JAR file (-g) \"$OPTARG\"!"
				echo
				exit 1
				fi
			;;
		l) # Maximum fraction of samples filtered at the genotype level
			if [[ $OPTARG =~ $DECNUMTEST ]]; then
			MAXFRACTFILTGENOT="$OPTARG"
			echo "Maximum fraction of samples filtered at the genotype level: $MAXFRACTFILTGENOT"
			echo
			else
				echo "Error! You did not provide correct number for maximum fraction of samples filtered at the genotype level (-l) \"$OPTARG\"!"
				echo
				exit 1
				fi
			;;
		w) # Minimal required coverage
			if [[ $OPTARG =~ $NUMTEST ]]; then
			GENOTFILTDP="$OPTARG"
			echo "Minimal required coverage: $GENOTFILTDP"
			echo
			else
				echo "Error! You did not provide correct number for minimal required coverage (-w) \"$OPTARG\"!"
				echo
				exit 1
				fi
			;;
		*)
			echo "Error! Unknown option!"
			echo "See usage options: \"$0 -h\""
			echo
			exit 1
			;;
		esac
	done

# Check if all required binaries are available
function toolcheck {
	command -v $1 >/dev/null 2>&1 || {
		echo >&2 "Error! $1 is required but not installed. Aborting. Please, install it."
		echo
		exit 1
		}
	}

toolcheck dirname
toolcheck bcftools
toolcheck R

# Checking if all required parameters are provided
if [ -z "$VCFFILE" ]; then
	echo "Error! Input joint VCF file (-f) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "$OUTNAME" ]; then
	echo "Base name of output join VCF (-n) was not set. Names of output files will start with name of input file."
	echo
	# Names of output files are derived from the input name
	VCFFILESNP=$(echo "$VCFFILEB" | sed 's/\(^.*\)vcf/\1snp.vcf/') # SNPs
	VCFFILESNPHARD=$(echo "$VCFFILESNP" | sed 's/\(^.*\)snp/\1hardfilter.snp/') # Hardfiltering SNPs
	VCFFILESNPPASS=$(echo "$VCFFILESNPHARD" | sed 's/\(^.*\)snp/\1snp.pass/') # Bi- and multiallelic SNPs
	VCFFILESNPBIAL=$(echo "$VCFFILESNPPASS" | sed 's/\(^.*\)snp\.pass/\1snp.pass.bial/') # Biallelic SNPs
	fi

if [ -z "$REF" ]; then
	echo "Error! Reference FASTA file (-a) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "$EXCLUDEPARALOGS" ]; then
	echo "Error! File with paralogs to exclude (-e) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "$JAVA" ]; then
	toolcheck java
	echo "Path to custom Java executable (-j) was not specified. Using default `which java`"
	JAVA=$(which java)
	echo
	fi

if [ -z "$JAVAMEM" ]; then
	echo "Java memory consumption for GATK (-m) was not set. Using default value of 12g."
	JAVAMEM="12g"
	echo
	fi

if [ -z "$GATK" ]; then
	echo "Error! Path to GATK JAR file (-g) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "$MAXFRACTFILTGENOT" ]; then
	echo "Maximum fraction of samples filtered at the genotype level (-l) was not set. Using default value of 0.1."
	echo
	MAXFRACTFILTGENOT='0.1'
	fi

if [ -z "$GENOTFILTDP" ]; then
	echo "Minimal required coverage (-w) was not set. Using default value of 8."
	echo
	GENOTFILTDP='8'
	fi

# Exit on error
function operationfailed {
	echo "Error! Operation failed!"
	echo
	echo "See previous message(s) to be able to trace the problem."
	echo
	exit 1
	}

echo "Start: `date`"
echo

# Copy reference sequence to working directory
echo "Copying reference sequence to working directory"
cp $REF* $VCFDIR/ && cp ${REF%.*}.dict $VCFDIR/ || operationfailed
echo

# Copy paralogs to exclude to working directory
echo "Copying paralogs to exclude to working directory"
cp $EXCLUDEPARALOGS $VCFDIR/ || operationfailed
echo

# Go to working directory containing joint VCF file
cd $VCFDIR

# Hardfilter SNPs from this joint file using recommended parameters for SNPS and exclude SNPs from set of putatively paralogue loci
echo "Extracting non-paralogous SNPs from joint VCF $VCFFILEB..."
echo
$JAVA -Xmx$JAVAMEM -jar $GATK -T SelectVariants -R $REFB -V $VCFFILEB -selectType SNP -o $VCFFILESNP -L scaffold_1 -L scaffold_2 -L scaffold_3 -L scaffold_4 -L scaffold_5 -L scaffold_6 -L scaffold_7 -L scaffold_8 --excludeIntervals $EXCLUDEPARALOGSB --excludeNonVariants --selectTypeToExclude INDEL || operationfailed
echo
echo "File with extracted non-paralogous SNPs was saved as $VCFFILESNP"
echo

# Create a filtering criterion
echo "Hard filtering SNPs from joint VCF $VCFFILESNP..."
echo
$JAVA -Xmx$JAVAMEM -jar $GATK -T VariantFiltration -R $REFB -V $VCFFILESNP -o $VCFFILESNPHARD --filterExpression "QD < 2.0" --filterName "Rad_filter1" --filterExpression "FS > 60.0" --filterName "Rad_filter2" --filterExpression "MQ < 40.0" --filterName "Rad_filter3" --filterExpression "MQRankSum < -12.5" --filterName "Rad_filter4" --filterExpression "ReadPosRankSum < -8.0" --filterName "Rad_filter5" --filterExpression "SOR < 3.0" --filterName "Rad_filter6" || operationfailed
echo
echo "Hard filtered SNPs were saved as $VCFFILESNPHARD"
echo

# Extract SNPs (both bi- and multiallelic) passing the AA_filter1 criterion
echo "Extracting only passing SNPs from joint VCF $VCFFILESNPHARD..."
echo
$JAVA -Xmx$JAVAMEM -jar $GATK -T SelectVariants -R $REFB -V $VCFFILESNPHARD -o $VCFFILESNPPASS --excludeFiltered || operationfailed
echo
echo "Extracted passed SNPs from previous filtering were saved as $VCFFILESNPPASS"
echo

# Extract BIALLELIC SNPs
echo "Extracting only BIALLELIC SNPs from joint VCF $VCFFILESNPPASS..."
echo
$JAVA -Xmx$JAVAMEM -jar $GATK -T SelectVariants -R $REFB -V $VCFFILESNPPASS -o $VCFFILESNPBIAL --excludeFiltered --restrictAllelesTo BIALLELIC || operationfailed
echo
echo "Biallelic SNPs were saved as $VCFFILESNPBIAL"
echo

# Set the filtering
echo "Marking filtered sites in the joint VCF $VCFFILESNPBIAL..."
echo
$JAVA -Xmx$JAVAMEM -jar $GATK -T VariantFiltration -R $REFB -V $VCFFILESNPBIAL -o $VCFFILESNPBIAL.dp$GENOTFILTDP.vcf.gz --genotypeFilterExpression "DP < $GENOTFILTDP" --genotypeFilterName DP-$GENOTFILTDP --setFilteredGtToNocall || operationfailed #
echo
echo "Marked filtered sites were saved as $VCFFILESNPBIAL.dp$GENOTFILTDP.vcf.gz"
echo

# Select variants based on this interval list (NB variants with < defined coverage will be still present in VCF)
echo "Selecting variants based on presence in $GENOTFILTDP of indivs in ${VCFFILESNPBIAL%.vcf.gz}.dp$GENOTFILTDP.vcf.gz joint VCF..."
echo
$JAVA -Xmx$JAVAMEM -jar $GATK -T SelectVariants -R $REFB -V $VCFFILESNPBIAL.dp$GENOTFILTDP.vcf.gz -o $VCFFILESNPBIAL.dp$GENOTFILTDP.perc$MAXFRACTFILTGENOT.vcf.gz --maxNOCALLfraction $MAXFRACTFILTGENOT || operationfailed # --maxFractionFilteredGenotypes $MAXFRACTFILTGENOT
echo
echo "Final selected variants were saved as $VCFFILESNPBIAL.dp$GENOTFILTDP.perc$MAXFRACTFILTGENOT.vcf.gz"
echo

echo "All output files are in directory \"$VCFDIR\" and names start with \"$OUTNAME*\""
echo

# Delete unneeded files
echo "Final cleanup"
rm $REFB* ${REFB%.*}.dict $EXCLUDEPARALOGSB || operationfailed
echo

# Statistics
echo "Statistics of SNPs in $VCFFILESNPBIAL.dp$GENOTFILTDP.perc$MAXFRACTFILTGENOT.vcf.gz using bcftools"
bcftools stats -s - $VCFFILESNPBIAL.dp$GENOTFILTDP.perc$MAXFRACTFILTGENOT.vcf.gz > $VCFFILESNPBIAL.dp$GENOTFILTDP.perc$MAXFRACTFILTGENOT.stats.txt || operationfailed
echo

# PCA and more statistics using R script
# Variable storing filename for R
export VCFR="$VCFFILESNPBIAL.dp$GENOTFILTDP.perc$MAXFRACTFILTGENOT.vcf.gz"
# Copy R script to working directory
cp $SCRIPTDIR/stat_pca_rad.r . || operationfailed
mkdir rpkgs
# Do the calculations in R
echo "Calculating statistics, PCAs and distances in $VCFFILESNPBIAL.dp$GENOTFILTDP.perc$MAXFRACTFILTGENOT.vcf.gz using R"
R CMD BATCH stat_pca_rad.r || operationfailed
echo
# Cleanup
rm -rf rpkgs stat_pca_rad.r || operationfailed
unset VCFR

cd $STARTDIR

echo "End: `date`"
echo

exit
