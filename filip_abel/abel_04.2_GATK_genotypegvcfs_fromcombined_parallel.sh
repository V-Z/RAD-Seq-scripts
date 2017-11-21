#!/bin/bash

# Job name:
#SBATCH --job-name=genotypegvcfs

# Project:
# dont change
#SBATCH --account=uio       

# Wall clock limit:
# not too high
#SBATCH --time=30:00:00  

# Max memory usage per core (MB) (NB usually 16cores*4G = 64G):
#SBATCH --mem-per-cpu=3G

# Number of cores:
#SBATCH --cpus-per-task=16

## Set up job environment      
## never change, it is sourcing an script
source /cluster/bin/jobsetup

## loading modules
module load gatk 

## Copy input files to the work directory:
cp $SUBMITDIR/cohort.g.vcf.gz* $SCRATCH
cp ~/aa/JIC_reference/*.* $SCRATCH

## Mark outfiles for automatic copying to $SUBMITDIR:
## copy all output files
chkfile "*raw.vcf*"

## Do some work:
cd $SCRATCH

# run Genotype GVCFs using one combined file (needed for > 200 *.g.vcf files)
# NB check that (i) -Xmx has cpus-per-task*mem-per-cpu (e.g. 16*3=48) and (ii) -nt = cpus-per-task
java -Xmx48g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-R alygenomes.fasta \
-V cohort.g.vcf.gz \
-nt 16 \
--includeNonVariantSites \
-o joint.raw.vcf.gz

