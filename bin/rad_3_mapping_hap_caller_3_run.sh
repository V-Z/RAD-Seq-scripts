#!/bin/bash

# Initialize variables
REF='' # Reference sequence
JAVA='' # PATH to custom Java binary
JAVAMEM='' # Memory limit for Picard and GATK
JAVAMEMTEST='^[0-9]+m$' # Testing if provided value is a number with k, m, g or t
PICARD='' # Path to directory containing Picard JAR file
PLATFORM='illumina' # Sequencing platform
GATKJ='' # Path to directory containing GATK JAR file
BAMLIST='' # List of BAM files in each directory
BAMFILE='' # Currently processed BAM file
BAMFILEBASE='' # Base name of currently processed sample
BAMFILERUNBASE='' # Base name of currently processed sample, run number
RUNNUMBER='' # Number of currently processed run
RGLB='' # RGLB parameter for Picard
RGPU='' # RGPU parameter for Picard
echo

# Parse initial arguments
while getopts "hrvf:a:j:m:p:g:" INITARGS; do
	case "${INITARGS}" in
		h) # Help and exit
			echo "Usage options:"
			echo -e "\t-h\tPrint this help and exit."
			echo -e "\t-r\tPrint references to used software and exit."
			echo -e "\t-v\tPrint script version and exit."
			echo -e "\t-f\tInput directory with FASTQ files saved as \"*.f*q*\"."
			echo -e "\t-a\tReference FASTA file."
			echo -e "\t-j\tOptional path to custom Java binary (default is output of \`command -v java\`; GATK requires Oracle Java)."
			echo -e "\t-m\tMaximal memory consumption allowed to Picard and GATK in MB per one CPU thread. Input as common for 'jar -Xmx', e.g. 12000m for '-Xmx12000m'. Default is 2000m. Warning! This value will be multiplied by number of CPU threads (-c)."
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
			echo "* Java, https://java.com/"
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
			if [ -d "${OPTARG}" ]; then
				COUNTFASTQ="$(find "${OPTARG}" -name "*.f*q*" 2>/dev/null | wc -l)"
				if [ "${COUNTFASTQ}" != 0 ]; then
					FQDIR="${OPTARG}"
					echo "Input directory: ${FQDIR}"
					echo
					else
						echo "Error! Given input directory does not contain any FASTQ files in *.f*q*!"
						echo
						exit 1
						fi
				else
					echo "Error! You did not provide path to input directory with FASTQ files (-f) \"${OPTARG}\"!"
					echo
					exit 1
					fi
			;;
		a) # Reference FASTA file
			if [ -r "${OPTARG}" ]; then
				REF="$(realpath "${OPTARG}")"
				echo "Reference FASTA file: ${REF}"
				echo
				else
					echo "Error! You did not provide path to reference FASTA file (-a) \"${OPTARG}\"!"
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
			echo "Maximal memory consumption by Java (Picard, GATK) in MB: ${JAVAMEM}"
			echo
			else
				echo "Error! You did not provide correct maximal memory consumption in MB (-m), e.g. 8000m or 6000m, \"${OPTARG}\"!"
				echo
				exit 1
				fi
			;;
		p) # Path to Picard JAR file
			if [ -r "${OPTARG}" ]; then
			PICARD="${OPTARG}"
			echo "Picard JAR file: ${PICARD}"
			echo
			else
				echo "Error! You did not provide path to Picard JAR file (-p) \"${OPTARG}\"!"
				echo
				exit 1
				fi
			;;
		g) # Path to GATK JAR file
			if [ -r "${OPTARG}" ]; then
			GATKJ="${OPTARG}"
			echo "GATK JAR file: ${GATKJ}"
			echo
			else
				echo "Error! You did not provide path to GATK JAR file (-g) \"${OPTARG}\"!"
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

toolcheck bwa
toolcheck samtools

# Checking if all required parameters are provided
if [ -z "${FQDIR}" ]; then
	echo "Error! Input directory with FASTQ files (-f) was not specified!"
	echo "See usage options: \"${0} -h\""
	echo
	exit 1
	fi

if [ -z "${REF}" ]; then
	echo "Error! Reference FASTA file (-a) was not specified!"
	echo "See usage options: \"${0} -h\""
	echo
	exit 1
	fi

if [ -z "${JAVA}" ]; then
	toolcheck java
	echo "Path to custom Java executable (-j) was not specified. Using default 'command -v java'"
	JAVA="$(command -v java)"
	echo
	fi

if [ -z "${JAVAMEM}" ]; then
	echo "Java memory consumption for Picard and GATK (-m) was not set. Using default value of 2000m."
	JAVAMEM='2000m'
	echo
	fi

if [ -z "${PICARD}" ]; then
	echo "Error! Path to Picard JAR file (-p) was not specified!"
	echo "See usage options: \"${0} -h\""
	echo
	exit 1
	fi

if [ -z "${GATKJ}" ]; then
	echo "Error! Path to GATK JAR file (-g) was not specified!"
	echo "See usage options: \"${0} -h\""
	echo
	exit 1
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

# If input files are compressed by BZIP2, uncompress them first
find . -name "*.f*q*.bz2" -print | grep -E '.*' > /dev/null && toolcheck parallel && echo && echo "Input files are compressed by BZIP2. Decompressing at $(date)" && find . -name "*.f*q*.bz2" -print | { parallel -j 2 bunzip2 -v '{}' || operationfailed; }

echo
echo "Mapping of paired and orphaned reads and postprocessing of output BAM files"

# Do the mapping separately paired files and the orphaned reads (reads without a mate)
echo "Starting mapping of paired reads at $(date)"
{ bwa mem "${REF}" ./*.dedup.R1.f*q* ./*.dedup.R2.f*q* | samtools view -bu | samtools sort -l 9 -o '{= s:dedup.+$:paired.bam: =}'; } || operationfailed
echo
echo "Starting mapping of orphaned reads (R1) at $(date)"
{ bwa mem "${REF}" ./*.unp.R1.f*q* | samtools view -bu | samtools sort -l 9 -o '{= s:unp.+$:unpaired.R1.bam: =}'; } || operationfailed
echo
echo "Starting mapping of orphaned reads (R2) at $(date)"
{ bwa mem "${REF}" ./*.unp.R2.f*q* | samtools view -bu | samtools sort -l 9 -o '{= s:unp.+$:unpaired.R2.bam: =}'; } || operationfailed
echo
echo "Finished mapping of paired and orphaned reads at $(date)"
echo

echo "Merging paired and unpaired BAM files, adding RG headers"
# First, merge the PE and SE files for each run and add the read group headers to each file
# Required headers are RGID, RGLB, RGPL, RGSM
BAMLIST=(*paired.bam)
# Process all BAM files
for ((BAMF=0; BAMF<${#BAMLIST[@]}; ++BAMF)); do
	BAMFILE=${BAMLIST[${BAMF}]} # E.g. AA016ac_paired_run02.bam
	BAMFILEBASE="$(basename "${BAMFILE%_run??_*_paired.bam}")" # the file base should be e.g. AA016ac
	BAMFILERUNBASE="${BAMFILE%.bam}" # the run base should be e.g. AA016ac_run02_paired
	# Merge all BAM files (should be 3) with the same file base
	echo "Merging all BAM files"
	echo
	"${JAVA}" -Xmx"${JAVAMEM}" -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar "${PICARD}" MergeSamFiles "$(printf 'INPUT=%s ' "${BAMFILEBASE}"*.bam)" OUTPUT="${BAMFILEBASE}".mergeRun.bam USE_THREADING=true || operationfailed
	echo
	# Extract run number for RGID assignment and add RG info headers
	RUNNUMBER="$(echo "${BAMFILERUNBASE}" | grep -o "run[[:digit:]][[:digit:]]")" # Get the run number from AA016ac_run02_paired, e.g. run02
	RGLB="$("${FQDIR}".lib1 | sed 's/_run[[:digit:]][[:digit:]]_[dipte]\{3\}//')"
	RGPU="${RUNNUMBER}".unit1
	RGSM="$("${FQDIR}" | sed 's/_run[[:digit:]][[:digit:]]_[dipte]\{3\}//')"
	echo "Modifying read groups"
	echo
	"${JAVA}" -Xmx"${JAVAMEM}" -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar "${PICARD}" AddOrReplaceReadGroups INPUT="${BAMFILEBASE}".mergeRun.bam OUTPUT="${BAMFILEBASE}".rg.bam RGID="${RUNNUMBER}" RGLB="${RGLB}" RGPL="${PLATFORM}" RGPU="${RGPU}" RGSM="${RGSM}" || operationfailed
	# Index this BAM
	echo
	echo "Indexing the BAM"
	echo
	"${JAVA}" -Xmx"${JAVAMEM}" -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar "${PICARD}" BuildBamIndex INPUT="${BAMFILEBASE}".rg.bam || operationfailed
	done
# Clean up all those temporary files
rm ./*.mergeRun.bam
echo
echo "Mapping ended at $(date)"
echo

# Produce genomic VCF
echo "Running HaplotypeCaller to produce genomic VCF files"
echo

if ls ./*dip*.rg.bam 1> /dev/null 2>&1; then # Diploids
	echo "Processing diploid at $(date)"
	"${JAVA}" -Xmx"${JAVAMEM}" -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar "${GATKJ}" -R "${REF}" -T HaplotypeCaller -I ./*dip*.rg.bam -ERC GVCF -ploidy 2 -o "${FQDIR}".raw.g.vcf.gz || operationfailed
	elif ls ./*tet*.rg.bam 1> /dev/null 2>&1; then # Tetraploids
		"${JAVA}" -Xmx"${JAVAMEM}" -Djava.io.tmpdir="${SCRATCHDIR}"/tmp -jar "${GATKJ}" -R "${REF}" -T HaplotypeCaller -I ./*tet*.rg.bam -ERC GVCF -ploidy 4 -o "${FQDIR}".raw.g.vcf.gz || operationfailed
		else
			echo "The name of the sample does not allow to find out if it is diploid or tetraploid!"
			echo
			operationfailed
	fi

# Calculating depth of coverage
echo "Calculating statistics of depth of coverage - it will be in file mapping_stats.txt"
{ echo "${FQDIR}"
	echo
	echo "Mapped paired reads"
	echo
	samtools flagstat ./*_paired.bam
	echo
	echo "Read groups"
	echo
	samtools flagstat ./*.rg.bam
	echo
	} > mapping_stats."${FQDIR}".txt || operationfailed
echo

echo "End: $(date)"
echo

exit

