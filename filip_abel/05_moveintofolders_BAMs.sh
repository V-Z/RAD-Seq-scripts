#!/bin/bash

# filename must be e.g. AA016ak_paired_run2.bam 
# each "samplename_" file PAIR will be moved into a separate directory named by its name (e.g. AA007A)

for i in *_paired_run?.bam; do
       echo $i
       file_base=${i%_paired_run?.bam}        # get AA016_a from AA016_a_run2_R1.fq.gz
       echo $file_base
       mkdir $file_base
       mv $file_base"_"* $file_base   # move there all BAM files
done

