#!/bin/bash

# Initialize variables to ensure they are empty
SUFIX='' # Default sufix of compressed FASTQ sequences
SUFIXTEST='^[a-zA-Z0-9._]+$' # Testing of validity of sufix of FASTQ files
DECOMPRESSER='' # Default archiver decompressing FASTQ files
COLUMNS='' # Testing of list of samples
SAMPLESLIST='' # List of samples (3 columns separated by TAB: names of individuals, barcodes, read IDs)
NCPU='' # Number of CPU threads for parallel processing
NCPUTEST='^[0-9]+$' # Testing if provided value is an integer
COUNTFASTQ='' # Test if input directory contains FASTQ files
OUTDIR='' # Output directory
FASTQINPUTDIR='' # Input directory with compressed FASTQ files to be processed
DIRTY='' # Do not delete temporal files and keep FASTQ files in output directories
echo

# Parse initial arguments
while getopts "hrvs:c:o:f:x:u:d" INITARGS; do
	case "${INITARGS}" in
		h) # Help and exit
			echo "Usage options:"
			echo -e "\t-h\tPrint this help and exit."
			echo -e "\t-r\tPrint references to used software and exit."
			echo -e "\t-v\tPrint script version and exit."
			echo -e "\t-s\tList of samples - TXT file with 3 columns separated by TAB: names of individuals, barcodes and read IDs."
			echo -e "\t-c\tNumber of CPU threads to use for parallel operations (parameter \"-j\" of GNU Parallel). If not provided, default is 2."
			echo -e "\t-o\tOutput directory. It should be empty."
			echo -e "\t-f\tInput directory with FASTQ files saved as \"*.${SUFIX}\"."
			echo -e "\t-x\tCustom sufix of compressed FASTQ sequences. Default is 'txt.gz'. String can contain numbers, letters, dots, or underscores"
			echo -e "\t-u\tCustom uncompressing command for compressed FASTQ files. Use 'zcat' (default) for gunzip (*.gz), 'bzcat' for bunzip2 (*.bz2), or 'lzcat' for LZMA (*.xz) archives."
			echo -e "\t-d\tDo not delete temporal files and keep FASTQ files in output directory sorted according to runs (not in original directory)."
			echo
			exit
			;;
		r) # References to cite and exit
			echo "Software to cite:"
			echo "* FASTX toolkit, http://hannonlab.cshl.edu/fastx_toolkit/"
			echo "* GNU Parallel, https://www.gnu.org/software/parallel/"
			echo
			exit
			;;
		v) # Print script version and exit
			echo "Version: 1.0"
			echo
			exit
			;;
		s) # List of samples (3 columns separated by TAB: names of individuals, barcodes, read IDs)
			if [ -r "${OPTARG}" ]; then
				COLUMNS="$(awk '{print NF}' "${OPTARG}" | sort -nu | tr "\n" " ")"
				if [ "${COLUMNS}" -eq 3 ]; then
					if grep -Pqv '^.+\t[ATCG]{5}\t[a-zA-Z0-9]+$' "${OPTARG}"; then
						echo "Error! List of samples (-s) \"${OPTARG}\" does not contain on all lines sample name, TAB, barcode (5 nucleotides), TAB and run ID!"
						echo
						exit 1
						else
							SAMPLESLIST="${OPTARG}"
							echo "List of samples: ${SAMPLESLIST}"
							echo
							fi
					else
						echo "Error! List of samples (-s) \"${OPTARG}\" does not have exactly three columns separated by TABs!"
						echo
						exit 1
						fi
				else
					echo "Error! List of samples (-s) \"${OPTARG}\" does not exist or is not readable!"
					echo
					exit 1
					fi
			;;
		c) # Number of CPU threads for parallel processing
			if [[ ${OPTARG} =~ ${NCPUTEST} ]]; then
				NCPU="${OPTARG}"
				echo "Number of CPU threads: ${NCPU}"
				echo
				else
					echo "Error! As number of CPU threads (-c) \"${OPTARG}\" you did not provide a number!"
					echo
					exit 1
					fi
			;;
		o) # Output directory
			if [ -d "${OPTARG}" ]; then
				OUTDIR="${OPTARG}"
				echo "Output directory: ${OUTDIR}"
				echo
				else
					echo "Output directory ${OUTDIR} doesn't exist - creating it."
					mkdir "${OUTDIR}" || { echo "Error! Can't create ${OUTDIR}!"; echo; exit 1; }
					fi
			;;
		f) # Input directory with compressed FASTQ files to be processed
			if [ -d "${OPTARG}" ]; then
				COUNTFASTQ="$(ls -1 "${OPTARG}"/*."${SUFIX}" 2>/dev/null | wc -l)"
				if [ "${COUNTFASTQ}" != 0 ]; then
					FASTQINPUTDIR="${OPTARG}"
					echo "Input directory: ${FASTQINPUTDIR}"
					echo
					else
						echo "Error! Given input directory does not contain any FASTQ files in *.${SUFIX}!"
						echo
						exit 1
						fi
				else
					echo "Error! You did not provide path to input directory with FASTQ files (-f) \"${OPTARG}\"!"
					echo
					exit 1
					fi
			;;
		x) # Default sufix of compressed FASTQ sequences
			if [[ ${OPTARG} =~ ${SUFIXTEST} ]]; then
				SUFIX="${OPTARG}"
				echo "Sufix of FASTQ files: ${SUFIX}"
				echo
				else
					echo "Error! As sufix of FASTQ files (-x) \"${OPTARG}\" you did not provide a valid string of number, letters, dots, or underscores!"
					echo
					exit 1
					fi
			;;
		u) # Uncompressing command for compressed FASTQ files
			case "${OPTARG}" in
				zcat) DECOMPRESSER='zcat';; # gunzip (*.gz) archives - default
				bzcat) DECOMPRESSER='bzcat';; # bunzip2 (*.bz2) archives
				lzcat) DECOMPRESSER='lzcat';; # LZMA (*.xz) archives
				*) echo
					echo "Error! Invalid uncompressing command (-u) \"${OPTARG}\"! Use 'zcat' for gunzip (*.gz), 'bzcat' for bunzip2 (*.bz2), or 'lzcat' for LZMA (*.xz) archives!"
					echo
					exit 1
				esac
			;;
		d) # Do not delete temporal files and keep FASTQ files in output directories
			DIRTY='TRUE'
			;;
		*)
			echo "Error! Unknown option!"
			echo "See usage options: \"$0 -h\""
			echo
			exit 1
			;;
		esac
	done

# Checking if all required parameters are provided
if [ -z "${SAMPLESLIST}" ]; then
	echo "Error! List of samples (-s) was not provided!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "${NCPU}" ]; then
	echo "Number of CPU threads (-c) for parallel operations was not set. Using default value of 2."
	echo
	NCPU='2'
	fi

if [ -z "${OUTDIR}" ]; then
	echo "Error! Output directory (-o) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "${FASTQINPUTDIR}" ]; then
	echo "Error! Input directory with FASTQ files (-f) was not specified!"
	echo "See usage options: \"$0 -h\""
	echo
	exit 1
	fi

if [ -z "${SUFIX}" ]; then
	echo "Sufix of FASTQ files (-x) was not provided. Using default 'txt.gz'."
	echo
	SUFIX='txt.gz'
	fi

if [ -z "${DECOMPRESSER}" ]; then
	echo "Uncompressing command for compressed FASTQ files (-u) was not set. Using default 'zcat'."
	echo
	DECOMPRESSER='zcat'
	fi

# Check if all required binaries are available
function toolcheck {
	command -v "${1}" >/dev/null 2>&1 || {
		echo >&2 "Error! $1 is required but not installed. Aborting. Please, install it."
		echo
		exit 1
		}
	}

toolcheck "${DECOMPRESSER}"
toolcheck parallel
toolcheck fastx_barcode_splitter.pl

# Exit on error
function operationfailed {
	echo "Error! Operation failed!"
	echo
	echo "See previous message(s) to be able to trace the problem."
	echo
	exit 1
	}

# # Check list of samples and if necessary, correct it
# # NOTE Modify it according to needs
# cp "${SAMPLESLIST}" "${SAMPLESLIST}".bak
# sed -i 's/^\([0-9]\)/AA\1/' "${SAMPLESLIST}"
# sed -i 's/\(AA[0-9]\{3\}\)\([a-z]\t\)/\1_\2/' "${SAMPLESLIST}"

echo "Start: $(date)"
echo

# Make output directories according to read IDs
# Move archived reads into their respective directories
echo "Making directories according to read IDs and moving respective sequences files into their directories"
cut -f 3 "${SAMPLESLIST}" | sort -u | parallel -j "${NCPU}" "echo '{}' && mkdir -p ${OUTDIR}/{}/demultiplexed && mv ${FASTQINPUTDIR}/*{}*.${SUFIX} ${OUTDIR}/{}/" || operationfailed
echo

# Process list of samples to extract bcfiles (names and barcodes) and save them into respective directories
echo "Creating barcode files"

while read -r FASTQREADGZ; do
	echo "Processing ${FASTQREADGZ}"
	grep "${FASTQREADGZ}" "${SAMPLESLIST}" | cut -f 1,2 > "${OUTDIR}"/"${FASTQREADGZ}"/barcodes.tsv || operationfailed
	done < "$(cut -f 3 "${SAMPLESLIST}" | sort -u)"
echo

# The demultiplexing
echo "Doing the demultiplexing. This may take longer time..."
echo
cut -f 3 "${SAMPLESLIST}" | sort -u | parallel -j $((NCPU-1)) "echo '{}' && date && ${DECOMPRESSER} ${OUTDIR}/{}/*1_sequence.${SUFIX} | fastx_barcode_splitter.pl --bcfile ${OUTDIR}/{}/barcodes.tsv --prefix ${OUTDIR}/{}/demultiplexed/ --suffix '_R1.fq' --bol --mismatches 1 && ${DECOMPRESSER} ${OUTDIR}/{}/*2_sequence.${SUFIX} | fastx_barcode_splitter.pl --bcfile ${OUTDIR}/{}/barcodes.tsv --prefix ${OUTDIR}/{}/demultiplexed/ --suffix '_R2.fq' --bol --mismatches 1 && echo" || operationfailed

# Rename the unmatched files to contain run ID
echo "Renaming the unmatched files to contain run ID"
while read -r DEMULTIPLEXEDDIR; do
	for UNMF in "${OUTDIR}"/"${DEMULTIPLEXEDDIR}"/demultiplexed/unmatched*; do
		echo -ne "Processing \"${DEMULTIPLEXEDDIR}\" : \"${UNMF}\"\r"
		mv "${UNMF}" "$(echo "${UNMF}" | sed "s/unmatched/unmatched_${DEMULTIPLEXEDDIR}/")" || operationfailed
		echo
		done
	done < "$(cut -f 3 "${SAMPLESLIST}" | sort -u)"
echo

# Move files of respective samples into their directories
echo "Moving of processed files into directories specific for each sample"
ls -1 "${OUTDIR}"/*/demultiplexed/*R1.fq | grep -v unmatched | sed 's/.\+\///' | sed 's/_R1.fq//' | parallel -j "${NCPU}" "echo '{}' && mkdir ${OUTDIR}/{} && mv ${OUTDIR}/*/demultiplexed/{}* ${OUTDIR}/{}/" || operationfailed
echo
# Move also unmatched files
echo "Moving unmatched files"
mv "${OUTDIR}"/*/demultiplexed/unmatched* "${OUTDIR}"/ || operationfailed
echo

# Compress all FASTQ files
echo "Compressing FASTQ files"
find "${OUTDIR}"/ -name "*.fq" -print | parallel -j "${NCPU}" "bzip2 -v -9 '{}'" || operationfailed
echo

# Cleanup of temporal files
if [ "${DIRTY}" == "TRUE" ]; then
	echo "Everything will be in ${OUTDIR} (no cleanup)."
	else
		echo "Final cleanup"
		# Remove barcode lists
		rm "${OUTDIR}"/*/barcodes.tsv
		# Move FASTQ files back to input directory
		mv "${OUTDIR}"/*/*."${SUFIX}" "${FASTQINPUTDIR}"/
		# Remove empty directories
		rmdir -p "${OUTDIR}"/*/demultiplexed 2>/dev/null
		fi
echo

echo "End: $(date)"
echo

exit
