#!/bin/bash

# Initialize variables
SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) # Determine script's directory
FQDIR='' # Input directory
NCPU='' # Number of CPU threads for parallel operations
NCPUTEST='^[0-9]+$' # Testing if provided value is an integer
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
# BAMHCFILE='' # Currently processed BAM file
# BAMHCFILEBASE='' # Base name of currently processed BAM file
echo

# Parse initial arguments
while getopts "hrvf:c:o:a:j:m:p:g:" INITARGS; do
	case "$INITARGS" in
		h) # Help and exit
			echo "Usage options:"
			echo -e "\t-h\tPrint this help and exit."
			echo -e "\t-r\tPrint references to used software and exit."
			echo -e "\t-v\tPrint script version and exit."
			echo -e "\t-f\tInput directory with FASTQ files saved as \"*.fq*\"."
			echo -e "\t-c\tNumber of CPU threads to use for parallel operations (parameter \"-j\" of GNU Parallel). If not provided, default is 3."
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

if [ -z "$NCPU" ]; then
	echo "Number of CPU threads (-c) for parallel operations was not set. Using default value of 3."
	echo
	NCPU='3'
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
find $OUTDIR/ -name "*.fq.bz2" -print | egrep '.*' > /dev/null && toolcheck parallel && echo && echo "Input files are compressed by BZIP2. Decompressing at `date`" && find $OUTDIR/ -name "*.fq.bz2" -print | parallel -j $NCPU "echo '{}' && bunzip2 '{}' || operationfailed"

echo
echo "Mapping of paired and orphaned reads and postprocessing of output BAM files"

echo
echo "Copying references to working directories at `date`"
find $OUTDIR -name "*trm_R1*" -print | parallel -j $NCPU "cp $REF* '{//}'/ && cp ${REF%.*}.dict '{//}'/ || operationfailed"

# Do the mapping separately paired files and the orphaned reads (reads without a mate)
echo
echo "Starting mapping of paired reads at `date`"
find $OUTDIR -name "*trm_R1*" -print | parallel -j $((NCPU-1)) "echo && echo '{}' && echo && bwa-0.7.3a mem '{//}'/$REFB '{//}'/*trm_R1* '{//}'/*trm_R2* | samtools view -bu | samtools sort -l 9 -o {= s:trm.+$:_paired.bam: =}' || operationfailed"

echo "Starting mapping of orphaned reads (R1) at `date`"
find $OUTDIR -name "*unp_R1*" -print | parallel -j $((NCPU-1)) "echo && echo '{}' && echo && bwa-0.7.3a mem '{//}'/$REFB '{}' | samtools view -bu | samtools sort -l 9 -o '{= s:unp.+$:_unpaired_R1.bam: =}' || operationfailed"

echo "Starting mapping of orphaned reads (R2) at `date`"
find $OUTDIR -name "*unp_R2*" -print | parallel -j $((NCPU-1)) "echo && echo '{}' && echo && bwa-0.7.3a mem '{//}'/$REFB '{}' | samtools view -bu | samtools sort -l 9 -o '{= s:unp.+$:_unpaired_R2.bam: =}' || operationfailed"
echo
echo "Finished mapping of paired and orphaned reads at `date`"
echo

# Remove original inputs from final destination (and trimming logs, if presented)
echo "Deleting unneeded files"
find $OUTDIR -name "*trm_R[12]*" -print | parallel -j $NCPU "rm '{}'"
find $OUTDIR -name "*unp_R[12]*" -print | parallel -j $NCPU "rm '{}'"
find $OUTDIR -name "*trim.log*" -print | parallel -j $NCPU "rm '{}'"
echo

echo "Merging paired and unpaired BAM files, adding RG headers"
for OUTDIRD in `find $OUTDIR -name '*.bam' -printf '%h\n' | sort -u`; do
	echo
	echo "Processing $OUTDIRD at `date`"
	echo
	cd $OUTDIRD/
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

# Diploids
echo "Processing diploids at `date`"
echo
find $OUTDIR -name "*.rg.bam" -print | grep 'dip' | parallel -j $((NCPU-1)) "echo && echo '{}' && echo && $JAVA -Xmx$JAVAMEM -jar $GATK -R {//}/$REFB -T HaplotypeCaller -I '{}' -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 --output_mode EMIT_VARIANTS_ONLY -ploidy 2 -o '{.}'.raw.g.vcf.gz || operationfailed"

# Tetraploids
echo "Processing tetraploids at `date`"
echo
find $OUTDIR -name "*.rg.bam" -print | grep 'tet' | parallel -j $((NCPU-1)) "echo && echo '{}' && echo && $JAVA -Xmx$JAVAMEM -jar $GATK -R {//}/$REFB -T HaplotypeCaller -I '{}' -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 --output_mode EMIT_VARIANTS_ONLY -ploidy 2 -o '{.}'.raw.g.vcf.gz || operationfailed"

# Deleting references in working directories
echo
echo "Deleting references from working directories"
find $OUTDIR -name "$REFB" -print | parallel -j $NCPU "rm '{//}'/$REFB* '{//}'/${REF%.*}.dict || operationfailed"
echo

# Calculating depth of coverage for each *.rg.bam file
echo "Calculating statistics of depth of coverage - they will be in file \"$OUTDIR/depth_of_coverage.tsv\""
echo
echo -e "Sample\tMin\t1stQu\tMedian\tMean\t3rdQu\tMax" > $OUTDIR/depth_of_coverage.tsv
for OUTDIRD in `find $OUTDIR -name "*.rg.bam" -print`; do
	echo "Processing $OUTDIRD"
	echo "$OUTDIRD" | sed 's/^.\+\///' >> $OUTDIR/depth_of_coverage.txt || operationfailed
	samtools depth $OUTDIRD | sed 's/[[:blank:]]\+/ /g' | cut -d ' ' -f 3 | Rscript $SCRIPTDIR/stat_depth_of_coverage.r | grep -v "^.*Min.*1st.*Qu.*Median.*Mean.*3rd.*Qu.*Max.*$" >> $OUTDIR/depth_of_coverage.txt || operationfailed
	done
echo

# Clean-up of the statistics
echo "Postprocessing file with statistics"
for L in `wc -l $OUTDIR/depth_of_coverage.txt`; do sed -i ':a;N;$!ba;s/\.rg\.bam\n//' $OUTDIR/depth_of_coverage.txt; done
cat $OUTDIR/depth_of_coverage.txt | sed 's/[[:blank:]]\+/ /g' | tr " " "\t" >> $OUTDIR/depth_of_coverage.tsv
rm $OUTDIR/depth_of_coverage.txt
echo

echo "End: `date`"
echo

exit
