#!/bin/bash

# Job name:
#SBATCH --job-name=combinegvcfs

# Project:
# dont change
#SBATCH --account=uio       

# Wall clock limit:
# not too high
#SBATCH --time=30:00:00  

# Max memory usage per core (MB) (NB usually 16cores*4G = 64G):
#SBATCH --mem-per-cpu=30G

## Set up job environment      
## never change, it is sourcing an script
source /cluster/bin/jobsetup

## loading modules
module load gatk 

## Copy input files to the work directory:
cp $SUBMITDIR/*.raw.g.vcf* $SCRATCH
cp ~/aa/JIC_reference/*.* $SCRATCH

## Mark outfiles for automatic copying to $SUBMITDIR:
## copy all output files
chkfile "cohort.g.vcf.gz*"

## Do some work:
cd $SCRATCH

#first prepare list of all samples
ls *.raw.g.vcf.gz | sed 's/^/-V /' | sed 's/$/ /' > join.vcf.samplelist.txt
SAMPLELIST=$(<join.vcf.samplelist.txt)
echo "Joint genotyping samples $SAMPLELIST"

#then run CombineGVCFs (needed if > 200 samples)
java -Xmx30g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T CombineGVCFs \
-R alygenomes.fasta \
$SAMPLELIST \
-o cohort.g.vcf.gz

