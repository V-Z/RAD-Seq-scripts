#!/bin/bash

# Initialize variables
VCFNAME='' # Default name (suffix) of input VCF files (e.g. "raw.g.vcf")
VCFCOMPRESSION='' # Default compression of input VCF files; if any, add "." on the beginning (e.g. ".gz")
VCFOUTSUFIX='' # Default suffix of output VCF file, add "." on the beginning (e.g. ".raw.vcf.gz")
STRINGTEST='^[a-zA-Z0-9_.]+$' # Testing if provided string contains only valid characters
VCFDIR='' # Input directory with VCF files
STARTDIR=$(pwd) # Current working directory
REF='' # Reference sequence
REFB='' # Only filename of the reference sequence
JAVAMEM='' # Memory limit for Picard and GATK
JAVAMEMTEST='^[0-9]+[kmgt]$' # Testing if provided value is a number with k, m, g or t
GATK='' # Path to directory containing GATK script file
OUTDIR='' # Output directory
JOINTNAME='' # Base name of the output file
TEMPDIR='' # Temporal directory for trimming
echo

# Parse initial arguments
while getopts "hrvw:u:x:f:o:n:a:m:g:" INITARGS; do
	case "$INITARGS" in
		h) # Help and exit
			echo "Usage options:"
			echo -e "\t-h\tPrint this help and exit."
			echo -e "\t-r\tPrint references to used software and exit."
			echo -e "\t-v\tPrint script version and exit."
			echo -e "\t-w\tName (suffix) of input VCF files (e.g. \"raw.g.vcf\")."
			echo -e "\t-u\tCompression of input VCF files; if any, add \".\" on the beginning (e.g. \".gz\")."
			echo -e "\t-x\tSuffix of output VCF file, add \".\" on the beginning (e.g. \".raw.vcf.gz\")."
			echo -e "\t-f\tInput directory with VCF files saved as e.g. \"*.raw.g.vcf.gz\"."
			echo -e "\t-o\tOutput directory. It should be empty."
			echo -e "\t-n\tBase name of output join VCF. Suffix (e.g. \".raw.vcf.gz\", see previous options) will be added. Allowed characters are letters, numbers, underscore or dot. If not provided, default \"join\" will be used (e.g. \"join.raw.vcf.gz\")."
			echo -e "\t-a\tReference FASTA file."
			echo -e "\t-m\tMaximal memory consumption allowed to GATK. Input as common for 'jar -Xmx', e.g. 12g for '-Xmx12g'. Default is 24g."
			echo -e "\t-g\tPath to GATK script file (GATK 4 and above is required)."
			echo "This script requires 2 CPU threads to run properly."
			echo
			exit
			;;
		r) # References to cite and exit
			echo "Software to cite:"
			echo "* GATK, https://software.broadinstitute.org/gatk/"
			echo "* GNU Parallel, https://www.gnu.org/software/parallel/"
			echo "* Java, https://java.com/ or http://openjdk.java.net/"
			echo "* Python, https://www.python.org/"
			echo
			exit
			;;
		v) # Print script version and exit
			echo "Version: 1.0"
			echo
			exit
			;;
		w) # Name (suffix) of input VCF files (default is "raw.g.vcf")
			if [[ $OPTARG =~ $STRINGTEST ]]; then
				VCFNAME="$OPTARG"
				echo "Name (suffix) of input VCF files: $VCFNAME"
				echo
				else
					echo "Error! As name (suffix) of input VCF files (-w) \"$OPTARG\" you did not provide a valid string containing only letters, numbers, dots and underscores (e.g. \"raw.g.vcf\")!"
					echo
					exit 1
					fi
			;;
		u) # Compression of input VCF files (default is ".gz")
			if [[ $OPTARG =~ $STRINGTEST ]]; then
				VCFCOMPRESSION="$OPTARG"
				echo "Compression of input VCF files: $VCFCOMPRESSION"
				echo
				else
					echo "Error! As compression of input VCF files (-u) \"$OPTARG\" you did not provide a valid string containing only letters, numbers, dots and underscores (e.g. \".gz\")!"
					echo
					exit 1
					fi
			;;
		x) # Sufix of output VCF file (default is ".raw.vcf.gz")
			if [[ $OPTARG =~ $STRINGTEST ]]; then
				VCFOUTSUFIX="$OPTARG"
				echo "Sufix of output VCF file: $VCFOUTSUFIX"
				echo
				else
					echo "Error! As suffix of output VCF file (-x) \"$OPTARG\" you did not provide a valid string containing only letters, numbers, dots and underscore (e.g. \".raw.vcf.gz\")!"
					echo
					exit 1
					fi
			;;
		f) # Input directory with compressed VCF files to be processed
			if [ -d "$OPTARG" ]; then
				COUNTVCF=`find $OPTARG -name "*.$VCFNAME*" -print 2>/dev/null | wc -l`
				if [ "$COUNTVCF" != 0 ]; then
					VCFDIR="$OPTARG"
					echo "Input directory: $VCFDIR"
					echo
					else
						echo "Error! Given input directory does not contain any VCF files in *.fq*!"
						echo
						exit 1
						fi
				else
					echo "Error! You did not provide path to input directory with VCF files (-f) \"$OPTARG\"!"
					echo
					exit 1
					fi
			;;
		o) # Output directory
			if [ -d "$OPTARG" ]; then
				OUTDIR="$OPTARG"
				echo "Output directory: $OUTDIR"
				echo
				else
					echo "Error! You did not provide path to output directory (-o) \"$OPTARG\"!"
					echo
					exit 1
					fi
			;;
		n) # Base name of the output file
			if [[ $OPTARG =~ $STRINGTEST ]]; then
				JOINTNAME="$OPTARG"
				echo "Name of the output file will be: $JOINTNAME$VCFOUTSUFIX"
				echo
				else
					echo "Error! You did not provide valid base name of output joint file (-n) \"$OPTARG\"!"
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
		g) # Path to GATK script file
			if [ -r "$OPTARG" ]; then
			GATK="$OPTARG"
			echo "GATK script file: $GATK"
			echo
			else
				echo "Error! You did not provide path to GATK script file (-g) \"$OPTARG\"!"
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

toolcheck basename
toolcheck java
toolcheck parallel
toolcheck python

# Checking if all required parameters are provided
if [ -z "$VCFDIR" ]; then
	echo "Error! Input directory with VCF files (-f) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "$OUTDIR" ]; then
	echo "Error! Output directory (-o) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "$JOINTNAME" ]; then
	echo "Base name of output join VCF (-n) was not set. Using default value \"join\"."
	echo
	JOINTNAME='join'
	fi

if [ -z "$REF" ]; then
	echo "Error! Reference FASTA file (-a) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "$JAVAMEM" ]; then
	echo "Java memory consumption for GATK (-m) was not set. Using default value of 24g."
	JAVAMEM='24g'
	echo
	fi

if [ -z "$GATK" ]; then
	echo "Error! Path to GATK script file (-g) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "$VCFNAME" ]; then
	echo "Name (suffix) of input VCF files (-w) was not set. Using default 'raw.g.vcf'."
	VCFNAME='raw.g.vcf'
	echo
	fi

if [ -z "$VCFCOMPRESSION" ]; then
	echo "Compression of input VCF files (-u) was not set. Using default '.gz'."
	VCFCOMPRESSION='.gz'
	echo
	fi

if [ -z "$VCFOUTSUFIX" ]; then
	echo "Sufix of output VCF file (-x) was not set. Using default '.raw.vcf.gz'."
	VCFOUTSUFIX='.raw.vcf.gz'
	echo
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

# Make temporal directory
echo "Making temporal directory in \"`pwd`/$OUTDIR\""
TEMPDIR=`mktemp -d $OUTDIR/VCFJOINXXXX`
echo

# Copy all VCF files to output directory
echo "Copying all individual VCF files into temporal directory \"`pwd`/$TEMPDIR\""
find $VCFDIR -name "*.$VCFNAME*" -print | parallel -j 2 "cp '{}' $TEMPDIR/" || operationfailed
echo

# Copy reference sequence to temporal directory
echo "Copying reference sequence to temporal directory"
cp $REF* $TEMPDIR/ && cp ${REF%.*}.dict $TEMPDIR/ || operationfailed
echo

# Go to output directory with VCF files to be processed
cd $TEMPDIR

# Create list of VCF files to be processed
SAMPLELIST=$(ls -1 *.$VCFNAME$VCFCOMPRESSION | sed 's/^/-V /' | sed 's/$/ /' | tr "\n" " ")

echo "Joining genotyping samples $SAMPLELIST at `date`"
echo
# Running Combine GVCFs
$GATK --java-options "-Xmx$JAVAMEM" CombineGVCFs -O $JOINTNAME.comb$VCFOUTSUFIX -R $REFB $SAMPLELIST --create-output-variant-index || operationfailed
echo

# Running Genotype GVCFs
$GATK --java-options "-Xmx$JAVAMEM" GenotypeGVCFs -O ../$JOINTNAME$VCFOUTSUFIX -R $REFB -V $JOINTNAME.comb$VCFOUTSUFIX --create-output-variant-index || operationfailed
echo
cd $STARTDIR

echo "Final cleanup"
rm -rf $TEMPDIR
echo

echo "End: `date`"
echo

exit
