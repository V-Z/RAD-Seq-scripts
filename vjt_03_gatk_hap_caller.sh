#!/bin/bash

# Initialize variables
SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) # Determine script's directory
FQDIR='' # Input directory
STARTDIR=$(pwd) # Current working directory
REF='' # Reference sequence
REFB='' # Only filename of the reference sequence
JAVA='' # PATH to custom Java binary
JAVAMEM='' # Memory limit for Picard and GATK
JAVAMEMTEST='^[0-9]+[kmgt]$' # Testing if provided value is a number with k, m, g or t
PICARD='' # Path to directory containing Picard JAR file
PLATFORM='illumina' # Sequencing platform
GATK='' # Path to directory containing GATK JAR file
OUTDIR='' # Output directory
BAMLIST='' # List of BAM files in each directory
BAMFILE='' # Currently processed BAM file
BAMFILEBASE='' # Base name of currently processed sample
BAMFILERUNBASE='' # Base name of currently processed sample, run number
RUNNUMBER='' # Number of currently processed run
RGLB='' # RGLB parameter for Picard
RGPU='' # RGPU parameter for Picard
BAMHCFILE='' # Currently processed BAM file
BAMHCFILEBASE='' # Base name of currently processed BAM file
echo

# Parse initial arguments
while getopts "hrvf:o:a:j:m:p:g:" INITARGS; do
	case "$INITARGS" in
		h) # Help and exit
			echo "Usage options:"
			echo -e "\t-h\tPrint this help and exit."
			echo -e "\t-r\tPrint references to used software and exit."
			echo -e "\t-v\tPrint script version and exit."
			echo -e "\t-f\tInput directory with FASTQ files saved as \"*.fq*\"."
			echo -e "\t-o\tOutput directory. It should be empty."
			echo -e "\t-a\tReference FASTA file."
			echo -e "\t-j\tOptional path to custom Java binary (default is output of \`which java\`; GATK requires Oracle Java)."
			echo -e "\t-m\tMaximal memory consumption allowed to Picard and GATK. Input as common for 'jar -Xmx', e.g. 12g for '-Xmx12g'. Default is 4g."
			echo -e "\t-p\tPath to Picard JAR file."
			echo -e "\t-g\tPath to GATK JAR file."
			echo
			exit
			;;
		r) # References to cite and exit
			echo "Software to cite:"
			echo "* BWA, http://bio-bwa.sourceforge.net/"
			echo "* GATK, https://software.broadinstitute.org/gatk/"
			echo "* GNU Parallel, https://www.gnu.org/software/parallel/"
			echo "* Java, https://java.com/ or http://openjdk.java.net/"
			echo "* Picard, https://broadinstitute.github.io/picard/"
			echo "* Samtools, http://samtools.sourceforge.net/"
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
				COUNTFASTQ=`ls -1 $OPTARG/*/*.fq* 2>/dev/null | wc -l`
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
			echo "Maximal memory consumption by Java (Picard, GATK): $JAVAMEM"
			echo
			else
				echo "Error! You did not provide correct maximal memory consumption (-m), e.g. 8g or 6000m, \"$OPTARG\"!"
				echo
				exit 1
				fi
			;;
		p) # Path to Picard JAR file
			if [ -r "$OPTARG" ]; then
			PICARD="$OPTARG"
			echo "Picard JAR file: $PICARD"
			echo
			else
				echo "Error! You did not provide path to Picard JAR file (-p) \"$OPTARG\"!"
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

toolcheck bwa
toolcheck samtools
toolcheck R

# Checking if all required parameters are provided
if [ -z "$FQDIR" ]; then
	echo "Error! Input directory with FASTQ files (-f) was not specified!"
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

if [ -z "$REF" ]; then
	echo "Error! Reference FASTA file (-a) was not specified!"
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
	echo "Java memory consumption for Picard and GATK (-m) was not set. Using default value of 4g."
	JAVAMEM='4g'
	echo
	fi

if [ -z "$PICARD" ]; then
	echo "Error! Path to Picard JAR file (-p) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "$GATK" ]; then
	echo "Error! Path to GATK JAR file (-g) was not specified!"
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

# Copy input files to the output destination for processing
echo "Copying files into final destination at `date`"
cp -a $FQDIR/* $OUTDIR/ || operationfailed

# If input files are compressed by BZIP2, uncompress them first
find $OUTDIR/ -name *.fq.bz2 -print | egrep '.*' > /dev/null && toolcheck parallel && echo && echo "Input files are compressed by BZIP2. Decompressing at `date`" && find $OUTDIR/ -name *.fq.bz2 -print | parallel -j 3 "echo '{}' && bunzip2 '{}'" || operationfailed

echo
echo "Mapping of paired and orphaned reads and postprocessing of output BAM files"

# Mapping of paired and orphaned reads
for OUTDIRD in `find $OUTDIR -name '*.fq*' -printf '%h\n' | sort -u`; do
	echo
	echo "Processing $OUTDIRD at `date`"
	echo
	# Copy reference sequence
	cp $REF* $OUTDIRD/ && cp ${REF%.*}.dict $OUTDIRD/ || operationfailed
	# Go to currently processed directory
	cd $OUTDIRD/
	# Do the mapping separately paired files and the orphaned reads (reads without a mate)
	echo "Starting mapping of paired and orphaned reads"
	echo
	bwa-0.7.3a mem $REFB *trm_R1* *trm_R2* | samtools view -bu | samtools sort -l 9 -o $(basename $OUTDIRD)"_paired.bam" || operationfailed
	bwa-0.7.3a mem $REFB *unp_R1* | samtools view -bu | samtools sort -l 9 -o $(basename $OUTDIRD)"_unpaired_R1.bam" || operationfailed
	bwa-0.7.3a mem $REFB *unp_R2* | samtools view -bu | samtools sort -l 9 -o $(basename $OUTDIRD)"_unpaired_R2.bam" || operationfailed
	# Remove original inputs from final destination (and trimming logs, if presented)
	rm *trm_R{1,2}* *unp_R{1,2}* *trim.log* 2>/dev/null
	echo
	echo "Finished mapping of paired and orphaned reads"
	echo
	echo "Merging paired and unpaired BAM files, adding RG headers"
	echo
	# First, merge the PE and SE files for each run and add the read group headers to each file
	# Required headers are RGID, RGLB, RGPL, RGSM
	BAMLIST=(*paired.bam)
	# Process all BAM files
	for ((BAMF=0; BAMF<${#BAMLIST[@]}; ++BAMF)); do
		BAMFILE=${BAMLIST[$BAMF]} # AA016ac_paired_run2.bam
		BAMFILEBASE=`basename ${BAMFILE%_run?_*_paired.bam}` # the file base should be e.g. AA016ac
		BAMFILERUNBASE=${BAMFILE%.bam} # the run base should be e.g. AA016ac_run2_paired
		# Merge all BAM files (should be 3) with the same file base
		echo "Merging all BAM files"
		echo
		$JAVA -Xmx$JAVAMEM -jar $PICARD MergeSamFiles $(printf 'INPUT=%s ' $BAMFILEBASE*.bam) OUTPUT=$BAMFILEBASE.mergeRun.bam || operationfailed
		echo
		# Extract run number for RGID assignment and add RG info headers
		RUNNUMBER=`echo $BAMFILERUNBASE | grep -o "run[[:digit:]]"` # Get the run number from AA016ac_run2_paired, e.g. run2
		RGLB=`basename $OUTDIRD.lib1`
		RGPU=$RUNNUMBER.unit1
		echo "Modifying read groups"
		echo
		$JAVA -Xmx$JAVAMEM -jar $PICARD AddOrReplaceReadGroups INPUT=$BAMFILEBASE.mergeRun.bam OUTPUT=$BAMFILEBASE.rg.bam RGID=$RUNNUMBER RGLB=$RGLB RGPL=$PLATFORM RGSM=$OUTDIRD RGPU=$RGPU || operationfailed
		# Index this BAM
		echo
		echo "Indexing the BAM"
		echo
		$JAVA -Xmx$JAVAMEM -jar $PICARD BuildBamIndex INPUT=$BAMFILEBASE.rg.bam || operationfailed
		done
	# Clean up all those temporary files
	rm *.mergeRun.bam
	cd $STARTDIR
	done
echo
echo "Mapping ended at `date`"
echo

# Produce genomic VCF
echo "Running HaplotypeCaller to produce genomic VCF files"
echo

# Run HaplotypeCaller on each DIPLOID sample to produce genomic VCF
for OUTDIRD in `find $OUTDIR -type d -print | grep 'dip'`; do
	cd $OUTDIRD
	BAMHCFILE=(*_rg.bam)
	BAMHCFILEBASE=${BAMHCFILE%_rg.bam}
	echo "Processing diploid $BAMHCFILEBASE at `date`"
	echo
	$JAVA -Xmx$JAVAMEM -jar $GATK -R $REFB -T HaplotypeCaller -I $BAMHCFILEBASE"_rg.bam" -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 --output_mode EMIT_VARIANTS_ONLY -ploidy 2 -o $BAMHCFILEBASE".raw.g.vcf.gz" || operationfailed
	rm $REFB* ${REFB%.*}.dict
	cd $STARTDIR
	echo
	done

# Run HaplotypeCaller on each TETRAPLOID sample to produce genomic VCF
for OUTDIRD in `find $OUTDIR -type d -print | grep 'tet'`; do
	cd $OUTDIRD
	BAMHCFILE=(*_rg.bam)
	BAMHCFILEBASE=${BAMHCFILE%_rg.bam}
	echo "Processing tetraploid $BAMHCFILEBASE at `date`"
	echo
	$JAVA -Xmx$JAVAMEM -jar $GATK -R $REFB -T HaplotypeCaller -I $BAMHCFILEBASE"_rg.bam" -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 --output_mode EMIT_VARIANTS_ONLY -ploidy 4 -o $BAMHCFILEBASE".raw.g.vcf.gz" || operationfailed
	rm $REFB* ${REFB%.*}.dict
	cd $STARTDIR
	echo
	done

# Calculating depth of coverage for each *_rg.bam file
echo "Calculating statistics of depth of coverage - they will be in file \"$OUTDIR/depth_of_coverage.tsv\""
echo
echo -e "Sample\tMin\t1stQu\tMedian\tMean\t3rdQu\tMax" > $OUTDIR/depth_of_coverage.tsv
for OUTDIRD in `find $OUTDIR -name *_rg.bam -print`; do
	echo "Processing $OUTDIRD"
	echo "$OUTDIRD" | sed -i 's/^.\+\///' >> $OUTDIR/depth_of_coverage.txt || operationfailed
	samtools depth $OUTDIRD | sed 's/[[:blank:]]\+/ /g' | cut -d ' ' -f 3 | Rscript $SCRIPTDIR/stat_depth_of_coverage.r >> $OUTDIR/depth_of_coverage.txt || operationfailed
	done
echo

# Clean-up of the statistics
echo "Postprocessing file with statistics"
sed -i ':a;N;$!ba;s/\_rg\.bam\n//' $OUTDIR/depth_of_coverage.txt
cat $OUTDIR/depth_of_coverage.txt | sed 's/[[:blank:]]\+/ /g' | tr " " "\t" | grep -v "^.*Min.*1st.*Qu.*Median.*Mean.*3rd.*Qu.*Max.*$" >> $OUTDIR/depth_of_coverage.tsv
rm $OUTDIR/depth_of_coverage.txt
echo

echo "End: `date`"
echo

exit
