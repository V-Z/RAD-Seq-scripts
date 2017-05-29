#!/bin/bash

# Initialize variables
FQDIR='' # Input directory - it MUST contain directories named according to individuals and containing F and R FASTQ files for respective individual
NCPU='' # Number of CPU threads for parallel operations
NCPUTEST='^[0-9]+$' # Testing if provided value is an integer
OUTDIR='' # Output directory
COUNTFASTQ='' # Test if input directory contains FASTQ files
FQDIRD='' # List of directories with FASTQ names
FQDIRDO='' # Ensure to have only the directory name, no path
TEMPDIR='' # Temporal directory for trimming
BASEIN1='' # Input file name 1 (F)
BASEIN2='' # Input file name 2 (R)
ADAPTOR='' # Path to FASTA file containing adaptor(s)
ADAPTORP='' # Absolute path to the adaptor file
ADAPTORF='' # Only filename of the adaptor file
JAVA='' # PATH to custom Java binary
JAVAMEM='' # Memory limit for Picard and GATK
JAVAMEMTEST='^[0-9]+[kmgt]$' # Testing if provided value is a number with k, m, g or t
TRIMMOMATIC='' # Path to Trimmomatic JAR file
STARTDIR=$(pwd) # Determine script's directory
echo

# Parse initial arguments
while getopts "hrvf:c:o:a:j:m:t:" INITARGS; do
	case "$INITARGS" in
		h) # Help and exit
			echo "Usage options:"
			echo -e "\t-h\tPrint this help and exit."
			echo -e "\t-r\tPrint references to used software and exit."
			echo -e "\t-v\tPrint script version and exit."
			echo -e "\t-f\tInput directory with FASTQ files saved as \"*.fq*\"."
			echo -e "\t-c\tNumber of CPU threads to use for parallel operations (parameter \"-j\" of GNU Parallel). If not provided, default is 2."
			echo -e "\t-o\tOutput directory. It should be empty."
			echo -e "\t-a\tFASTA file with adaptors."
			echo -e "\t-j\tOptional path to custom Java binary (default is output of \`which java\`)."
			echo -e "\t-m\tMaximal memory consumption allowed to Picard and GATK. Input as common for 'jar -Xmx', e.g. 12g for '-Xmx12g'. Default is 2g."
			echo -e "\t-t\tTrimmomatic JAR file."
			echo
			exit
			;;
		r) # References to cite and exit
			echo "Software to cite:"
			echo "* GNU Parallel, https://www.gnu.org/software/parallel/"
			echo "* Java, https://java.com/ or http://openjdk.java.net/"
			echo "* Trimmomatic, http://www.usadellab.org/cms/?page=trimmomatic"
			echo
			exit
			;;
		v) # Print script version and exit
			echo "Version: 1.0"
			echo
			exit
			;;
		f) # Input directory with compressed FASTQ files to be processed
			if [ -d "$OPTARG" ]; then
				COUNTFASTQ=$(ls -1 $OPTARG/*/*.fq* 2>/dev/null | wc -l)
				if [ "$COUNT" != 0 ]; then
					FQDIR="$OPTARG"
					echo "Input directory: $FQDIR"
					echo
					else
						echo "Error! Given input directory does not contain any FASTQ files in *.fq*!"
						echo
						exit 1
						fi
				else
					echo "Error! You did not provide path to input directory with FASTQ files (-f) \"$OPTARG\"!"
					echo
					exit 1
					fi
			;;
		c) # Number of CPU threads for parallel processing
			if [[ $OPTARG =~ $NCPUTEST ]]; then
				NCPU="$OPTARG"
				echo "Number of CPU threads: $NCPU"
				echo
				else
					echo "Error! As number of CPU threads (-c) \"$OPTARG\" you did not provide a number!"
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
		a) # FASTA file containing adaptor(s)
			if [ -r "$OPTARG" ]; then
				ADAPTOR="$OPTARG"
				# Decompose path to adaptor FASTA file to absolute path to directory and filename
				ADAPTORP=$(realpath $ADAPTOR) # Absolute path to the adaptor file
				ADAPTORF=$(basename $ADAPTOR) # Only filename of the adaptor file
				echo "Adaptor(s) FASTA file: $ADAPTOR"
				echo
				else
					echo "Error! You did not provide path to FASTA file with adaptor(s) (-a) \"$OPTARG\"!"
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
			echo "Maximal memory consumption by Java (Trimmomatic): $JAVAMEM"
			echo
			else
				echo "Error! You did not provide correct maximal memory consumption (-m), e.g. 8g or 6000m, \"$OPTARG\"!"
				echo
				exit 1
				fi
			;;
		t) # Path to Trimmomatic JAR file
			if [ -r "$OPTARG" ]; then
			TRIMMOMATIC="$OPTARG"
			echo "Trimmomatic JAR file: $TRIMMOMATIC"
			echo
			else
				echo "Error! You did not provide path to Trimmomatic JAR file (-t) \"$OPTARG\"!"
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

toolcheck parallel
toolcheck xargs

# Checking if all required parameters are provided
if [ -z "$FQDIR" ]; then
	echo "Error! Input directory with FASTQ files (-f) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "$NCPU" ]; then
	echo "Number of CPU threads (-c) for parallel operations was not set. Using default value of 2."
	echo
	NCPU='2'
	fi

if [ -z "$OUTDIR" ]; then
	echo "Error! Output directory (-o) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "$ADAPTOR" ]; then
	echo "Error! Adaptor file (-a) was not specified!"
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
	echo "Java memory consumption for Trimmomatic (-m) was not set. Using default value of 2g."
	JAVAMEM='2g'
	echo
	fi

if [ -z "$TRIMMOMATIC" ]; then
	echo "Error! Path to Trimmomatic JAR file (-t) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
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

# Make output directories
echo "Making output directories"
find $FQDIR -name '*.fq*' -printf '%h\n' | sort -u | parallel -j $NCPU "echo '{}' && mkdir -p $OUTDIR/{/}" || operationfailed
echo

# Initialize file with statistics
echo "Initializing file with statistics"
printf "indiv%s\t%sR1_trm_i.e.N_read_PAIRS%s\t%sR1_unpaired_reads%s\t%sR2_unpaired_reads%s\n%s" > $OUTDIR/trimming_info.tsv || operationfailed
echo

# Trim the sequences
# Every directory is named according to individual, every directory contains only files for respective individual
# Remove first 8 bases at the beginning of each read, keep reads using default threshold 20 QC score, remove reads shorter than 100 bp, trim some low quality ends
echo "Trimming. This may take longer time..."
echo
for FQDIRD in `find $FQDIR -name '*.fq*' -printf '%h\n' | sort -u`; do
	# Ensure to get name of the directory containing pair of sequences only (not full path)
	FQDIRDO=$(basename $FQDIRD)
	# Go to working directory
	cd $FQDIRD
	# Make temporal output directory
	echo "Making temporal directory"
	TEMPDIR=$(mktemp -d TRMXXX) || operationfailed
	echo
	# Get first and second FASTQ in the respective input directory as input files
	BASEIN1=$(ls -1 *.fq* | head -1)
	BASEIN2=$(ls -1 *.fq* | head -2 | tail -1)
	# Copy the adaptor
	echo "Copying adaptor"
	cp $ADAPTORP . || operationfailed
	echo

	echo "Processing $BASEIN1 and $BASEIN2 at `date`"
	echo "Trimming"
	# Do the trimming
	$JAVA -Xmx$JAVAMEM -jar $TRIMMOMATIC PE -threads $((NCPU-1)) -phred33 -trimlog $TEMPDIR/$FQDIRDO"_trim.log" $BASEIN1 $BASEIN2 $TEMPDIR/$FQDIRDO"_trm_R1.fq" $TEMPDIR/$FQDIRDO"_unp_R1.fq" $TEMPDIR/$FQDIRDO"_trm_R2.fq" $TEMPDIR/$FQDIRDO"_unp_R2.fq" HEADCROP:8 ILLUMINACLIP:$ADAPTORF:2:30:10 SLIDINGWINDOW:4:15 MINLEN:100 || operationfailed
	# Remove the adaptor
	rm $ADAPTORF

	# Statistics
	echo "Statistics"
	printf '%s\t%s' $TEMPDIR/$FQDIRDO >> $STARTDIR/$OUTDIR/trimming_info.tsv  # print sample name
	cat $TEMPDIR/$FQDIRDO"_trm_R1.fq" | echo $((`wc -l`/4)) | xargs printf '%s\t%s' >> $STARTDIR/$OUTDIR/trimming_info.tsv # No of R1 reads (N of lines /4)
	cat $TEMPDIR/$FQDIRDO"_unp_R1.fq" | echo $((`wc -l`/4)) | xargs printf '%s\t%s' >> $STARTDIR/$OUTDIR/trimming_info.tsv # No of R1 unpaired reads (N of lines /4)
	cat $TEMPDIR/$FQDIRDO"_unp_R2.fq" | echo $((`wc -l`/4)) | xargs printf '%s\t%s' >> $STARTDIR/$OUTDIR/trimming_info.tsv # No of R2 unpaired reads (N of lines /4)
	printf '%s\n%s' >> $STARTDIR/$OUTDIR/trimming_info.tsv
	echo

	# Compress output files with bzip2
	echo "Compressing output files"
	parallel -j $NCPU -X bzip2 -9 ::: $TEMPDIR/* || operationfailed
	echo

	# Moving files to final output directory
	echo "Moving files to final destination"
	mv $TEMPDIR/* $STARTDIR/$OUTDIR/$FQDIRDO/ || operationfailed
	rmdir $TEMPDIR
	cd $STARTDIR
	echo
	done

# Remove names of temporal directories and whitespaces on the end of lines from the statistics file
sed -i 's/^.*TRM.\{3\}\///' $STARTDIR/$OUTDIR/trimming_info.tsv
sed -i 's/[[:blank:]]\+$//' $STARTDIR/$OUTDIR/trimming_info.tsv

echo "End: `date`"
echo

exit
