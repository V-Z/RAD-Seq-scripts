#!/bin/bash

# filename must be SAMPLENAME_RUNNUMBER_R1.fq OR SAMPLENAME_RUNNUMBER_R1.fq (e.g., AA016_a_run2_R1.fq) 
# each sample file PAIR will be moved into a separate directory named by its name (e.g. AA007A)

for i in *R1.fq; do
       echo $i
       file_base=${i%_run2_R1.fq}        # get AA016_a from AA016_a_run2_R1.fq.gz
       gzip $i
       echo $file_base
       mkdir $file_base
       gzip $file_base"_run2_R2.fq"
       mv $file_base"_run"* $file_base   # move there both R1 and R2 files
done

