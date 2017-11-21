#!/bin/bash

# files must be sorted bams named as e.g.

#     AA229_i_run4_dip_paired.bam
#     AA229_i_run4_dip_unpaired_R1.bam
#     AA229_i_run4_dip_unpaired_R2.bam

# files per each sample must be stored in a separate directory named by its name (e.g. AA016ac) =  it is the TaxonDir variable ($1 means just first parameter of the script)
# there could be several BAM files in the individual's folder, but differing by the _runX_ number)
# usage: bash abel_02_GATK_joint_preprocessing_onedir.sh AA016ac


PLATFORM="illumina"
PICARD_HOME="/usit/abel/u1/filipko/programs/picard/1.119"
GATK_HOME="/usit/abel/u1/filipko/programs/GATK/3.5"
REFSEQ="/usit/abel/u1/filipko/aa/JIC_reference/alygenomes.fasta"
# NB .fai and .dict file must be created using picard and samtools
# java -jar /usit/abel/u1/filipko/programs/picard/1.119/CreateSequenceDictionary.jar R= alygenomes.fasta O= alygenomes.dict
# samtools faidx alygenomes.fasta


TaxonDir=$1
echo "$TaxonDir"
cd $TaxonDir
BAMLIST=(*paired.bam)

#First, merge the PE and SE files for each run and add the read group headers to each file
#Required headers are RGID, RGLB, RGPL, RGSM

echo "Merging paired and unpaired BAM files, adding RG headers..."
for ((i=0;i<${#BAMLIST[@]};++i)); do
        echo ${BAMLIST[$i]}
        ThisBam=${BAMLIST[$i]}                                  # AA016ac_paired_run2.bam
        BAM_file_base=${ThisBam%_run?_???_paired.bam}           # the file base should be e.g. AA229_i
	BAM_run_base=${ThisBam%_???_paired.bam}                 # the run base should be e.g. AA229_i_run4
        echo $BAM_file_base
        echo $BAM_run_base                                                              
        #merge ALL bam files (should be 3) with the same file base
        java -Xmx12g -jar $PICARD_HOME/MergeSamFiles.jar $(printf 'INPUT=%s ' $BAM_file_base*.bam) OUTPUT=$BAM_file_base.mergeRun.bam
 
  
        #extract run number for RGID assignment and add RG info headers
        ThisRunNumber=${BAM_run_base##*_}               #get the bit after the last "_" from AA229_i_run4, e.g. run4
        echo $ThisRunNumber
        thisRGLB=$TaxonDir.lib1 
        thisRGPU=$ThisRunNumber.unit1
        java -Xmx12g -jar $PICARD_HOME/AddOrReplaceReadGroups.jar INPUT=$BAM_file_base.mergeRun.bam OUTPUT=$BAM_file_base.rg.bam RGID=$ThisRunNumber RGLB=$thisRGLB RGPL=$PLATFORM RGSM=$TaxonDir RGPU=$thisRGPU
	# index this BAM        
	java -Xmx12g -jar $PICARD_HOME/BuildBamIndex.jar INPUT=$BAM_file_base.rg.bam

        #Next run RealignerTargetCreator and IndelRealigner on these bam files
	echo "Realigning indels..."
        java -Xmx12g -jar $GATK_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REFSEQ -I $BAM_file_base.rg.bam -o $BAM_file_base.intervals

        java -Xmx12g -jar $GATK_HOME/GenomeAnalysisTK.jar -T IndelRealigner -R $REFSEQ -I $BAM_file_base.rg.bam -targetIntervals $BAM_file_base.intervals -o $BAM_file_base.realign.bam

	# index this BAM        
	java -Xmx12g -jar $PICARD_HOME/BuildBamIndex.jar INPUT=$BAM_file_base.realign.bam

done

#Clean up all those temporary files
rm *.mergeRun.bam*
rm *.rg.bam*
rm *.rg.bai
rm *.intervals

#If more runs per this sample, finally merge all pre-processed bam files for each individual
#echo "Merging pre-processed bam files for $TaxonDir..."
#java -Xmx12g -jar $PICARD_HOME/MergeSamFiles.jar $(printf 'INPUT=%s ' *.realign.bam) OUTPUT=$TaxonDir.merged.bam 
#java -Xmx12g -jar $PICARD_HOME/BuildBamIndex.jar INPUT=$TaxonDir.merged.bam

echo "Finished processing files in $TaxonDir"

