#!/bin/bash

# Author: Vojtěch Zeisek, https://trapa.cz/
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# 

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# Checking if exactly one variables is provided
if [ "$#" -ne '1' ]; then
	echo "Error! Exactly 1 parameter (directory with mapped samples) is required! $# parameters received."
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

# Switching to working directory
echo "Going to $1"
cd "$1" || operationfailed
echo

# Mapping statistics
echo "Initializing files with statistics"
printf "Sample name\tNumber of QC passed reads\tNumber of mapped reads\n" > mapping_stats.tsv || operationfailed
for MS in mapping_stats.*.txt; do
	echo "Processing ${MS}"
	{ head -n 1 "${MS}" | xargs printf '%s\t%s' # Sample name
		grep -A 2 "Read groups" "${MS}" | grep -o "^[0-9]\+" | xargs printf '%s\t%s' # Number of QC passed reads
		grep -A 6 "Read groups" "${MS}" | tail -n 1 | grep -o "^[0-9]\+" | xargs printf '%s\t%s' # Number of mapped reads
		printf '\n'
		} >> mapping_stats.tsv || operationfailed
	done
echo

# Sorting into subdirectories according to samples names
# Creating directories
echo "Creating directories according to sample names"
mkdir $(find . -name "mapping_hap_caller.*" | sed 's/^\.\/mapping_hap_caller\.//' | sed 's/\.log$//' | tr "\n" " ") || operationfailed
echo
# Moving files into respective directories
for D in $(find . -name "mapping_hap_caller.*" | sed 's/^\.\/mapping_hap_caller\.//' | sed 's/\.log$//'); do
	echo "Processing ${D}"
	mv "${D}"*.{raw.g.vcf.gz,raw.g.vcf.gz.tbi,bai,bam} mapping_hap_caller."${D}".log mapping_stats."${D}".txt RADSeq_mapping_hapcaller."${D}".* "${D}"/
	done
echo

exit

