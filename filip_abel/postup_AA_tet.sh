# storage on Abel at
/projects/researchers/researchers01/filipko

# calculate number of reads  in each file
for file in *.fastq; do echo $file; cat $file | echo $((`wc -l`/4)); done
for file in *.fq; do wc -l $file ; done        # this must be divided by 4

####################################
### 1) DEMULTIPLEXING
####################################
# samples with different INDEX barcode (that in the fastq header) are already splitted into different fastqs, but we need to split the individuals by the INLINE barcode
# inline barcode is in both R1 and R2

# split into separate folders (e.g raw1 ... raw4 according to teh INDEX barcoded files), files with corresponding pairs into one folder and gunzip

# barcode file - example, first column is newly defined file=sample name
# must not end with emty line
AA016_a	GCATG
AA086_c	AACCA
AA090_c	CGATC
AA096ac	TCGAT
...

# separately for each file
module load fastx-toolkit
cat raw1/C7HW9ANXX_ddRAD_Tatry_16s002650-1-1_Fuxova_lane5Tatry1_1_sequence.fastq | fastx_barcode_splitter.pl --bcfile barcodes_AA_tet_run2_raw1.txt --prefix ./samples/ --suffix "_R1.fq" --bol --mismatches 1 
cat raw1/C7HW9ANXX_ddRAD_Tatry_16s002650-1-1_Fuxova_lane5Tatry1_2_sequence.fastq | fastx_barcode_splitter.pl --bcfile barcodes_AA_tet_run2_raw1.txt --prefix ./samples/ --suffix "_R2.fq" --bol --mismatches 1 

# all via forloop
module load fastx-toolkit
mkdir ./samples
for dir in raw*; do
	echo "processing $dir"
	cd $dir
	gunzip *.gz
        cat *1_sequence.txt | fastx_barcode_splitter.pl --bcfile barcodes_AA_tet_run2_$dir.txt --prefix ../samples/ --suffix "_run2_R1.fq" --bol --mismatches 1 
        cat *2_sequence.txt | fastx_barcode_splitter.pl --bcfile barcodes_AA_tet_run2_$dir.txt --prefix ../samples/ --suffix "_run2_R2.fq" --bol --mismatches 1 
	cd ..
done 


# Then gzip all sample files and move into folders named by the sample name (e.g. AA016_a_R1.fq.gz AND AA016_a_R2.fq.gz into folder AA016_a)
# if necessary also rename
rename s/_R/_run2_R/ AA*  #locally
rename _R _run2_R AA*     #on cluster

# forloop moving into folders
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

####################################
### 2) TRIMMING
####################################

### Trimmomatic
# remove first 8 bases at the beginning of each read, keep reads using default threshold 20 QC score, remove reads shorter than 100 bp, trim some low quality ends
# also save it to a specified folder
# one sample
java -jar ~/programs/trimmomatic/trimmomatic-0.36.jar PE -phred33 AA016_a_run2_R1.fq.gz AA016_a_run2_R2.fq.gz ~/RADseq/data_AA/AA_tet_run2_tatry/trimmed/AA016_a_R1_trm.fq ~/RADseq/data_AA/AA_tet_run2_tatry/trimmed/AA016_a_R1_unpaired.fq ~/RADseq/data_AA/AA_tet_run2_tatry/trimmed/AA016_a_R2_trm.fq ~/RADseq/data_AA/AA_tet_run2_tatry/trimmed/AA016_a_R2_unpaired.fq HEADCROP:8 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:100

# forloop on cluster (specific dirs)
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# I have each fq pair in specific dir named by the sample name
# remove first 8 bases at the beginning of each read, keep reads using default threshold 20 QC score, remove reads shorter than 100 bp, trim some low quality ends
# change the paths if necessary
for dir in AA*; do
	echo "processing $dir"
        mkdir ~/aa/GATK/tetraploids_june_2016/trimmed/$dir   # change this path if necessary
	cp ~/programs/trimmomatic/adapters/TruSeq3-PE.fa ./$dir
	cd $dir
        java -jar ~/programs/trimmomatic/trimmomatic-0.36.jar PE -phred33 $dir"_run2_R1.fq.gz" $dir"_run2_R2.fq.gz" ~/aa/GATK/tetraploids_june_2016/trimmed/$dir/$dir"_trm_run2_R1.fq.gz" ~/aa/GATK/tetraploids_june_2016/trimmed/$dir/$dir"_unpaired_run2_R1.fq.gz" ~/aa/GATK/tetraploids_june_2016/trimmed/$dir/$dir"_trm_run2_R2.fq.gz" ~/aa/GATK/tetraploids_june_2016/trimmed/$dir/$dir"_unpaired_run2_R2.fq.gz" HEADCROP:8 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:100  	
	cd ..
done
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


### calculate stats (N all reads)
nohup bash /usit/abel/u1/filipko/aa/GATK/tetraploids_june_2016/trimmed_backup/04_trimmingstats.sh > nohup_AA_triminfo.out&

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#calculate stats (N all reads)
touch trimming_info.csv
printf "indiv %s\t%s R1_trm_i.e.N_read_PAIRS %s\t%s R1_unpaired_reads %s\t%s R2_unpaired_reads %s\n%s"  >> trimming_info.csv
for dir in AA*
do
echo "processing $dir"
printf '%s\t%s' $dir >> trimming_info.csv  # print sample name
#echo $file |xargs printf '%s\t%s' >> mapping_info.csv   # print file name
less $dir/$dir"_trm_run2_R1.fq.gz" | echo $((`wc -l`/4)) |xargs printf '%s\t%s' >> trimming_info.csv      # No of R1 reads (N of lines /4)
less $dir/$dir"_unpaired_run2_R1.fq.gz" | echo $((`wc -l`/4)) |xargs printf '%s\t%s' >> trimming_info.csv      # No of R1 unpaired reads (N of lines /4)
less $dir/$dir"_unpaired_run2_R2.fq.gz" | echo $((`wc -l`/4)) |xargs printf '%s\t%s' >> trimming_info.csv      # No of R2 unpaired reads (N of lines /4)
printf '%s\n%s' >> trimming_info.csv
done
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


####################################
# 4) MAPPING 
####################################
### map reads of all samples to the refernce using bwa-mem
# ref path
/usit/abel/u1/filipko/aa/JIC_reference/alygenomes.fasta

# first index your reference
module load bwa
bwa index alygenomes.fasta

# rename folders to add sequential numbers at the end (xxx.1, xxx.2 etc)
n=1; for dir in AA*; do mv "$dir" "$(printf "$dir.$n")";   n=$((n+1)); done

# then map, sort and convert to BAM
# we can use gzipped files
# to use the array script, append a rising number to the end of each dir (the.1 should vary from .1 ... .N across the dirs)
# I use following two scripts:
abel_02_bwa_mem_submit.sh
abel_02_bwa_mem_worker.sh

sbatch abel_02_bwa_mem_submit.sh

#  xxxxxxxxxxxxxxxxxxxx PART OF THE SCRIPT WITH THE REAL WORK  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

bwa mem $REF $FASTQDIR/*_trm_run2_R1* $FASTQDIR/*_trm_run2_R2* | samtools view -bS - | samtools sort - $OUTFILE"_paired_run2"
bwa mem $REF $FASTQDIR/*_unpaired_run2_R1*  | samtools view -bS - | samtools sort - $OUTFILE"_unpairedR1_run2"
bwa mem $REF $FASTQDIR/*_unpaired_run2_R2*  | samtools view -bS - | samtools sort - $OUTFILE"_unpairedR2_run2"

#  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


###index all bams
for i in *bam; do samtools index $i; done

###calculate stats (% mapped reads)
nohup bash /usit/abel/u1/filipko/aa/GATK/tetraploids_june_2016/trimmed/calculate_bamstats.sh > nohup_AA_bamstats.out&

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
touch mapping_info.csv
printf "indiv %s\t%s all_reads %s\t%s mapped %s\n%s"  >> mapping_info.csv
for file in *paired_run1.bam
do
echo "processing $file"
printf '%s\t%s' $file >> mapping_info.csv
samtools view $file | wc -l |xargs printf '%s\t%s' >> mapping_info.csv           # No of all reads
samtools view -F 4 $file | wc -l |xargs printf '%s\t%s' >> mapping_info.csv      # No of mapped reads
printf '%s\n%s' >> mapping_info.csv
done
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

## NB now possible to view in IGV
java -Xmx2g -jar /home/filip/programs/IGV/igv.jar

# useful things to set in IGV
#rightclick in the left ... 
#    - expanded/squished
#    - collor alignmnets by ...
#    - view as pairs
#yellow box on the upper bar - tick: show details on click - opens info on the reads only after clicking (not popping automatically)
#save session ... saves it as xml file for opening in next time exactly the same

#check multiple individuals using the VCF file not BAM


################################################
#### 5) BAMS PREPROCESSING AND VARIANT CALLING
#################################################


# 5.1) move BAMs into folders named by the corresponding sample name name using forloop
######################################################################################

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# filename must be e.g. AA016ak_paired_run2.bam 
# each "samplename_" file PAIR will be moved into a separate directory named by its name (e.g. AA007A)
for i in *_paired_run?.bam; do
       echo $i
       file_base=${i%_paired_run?.bam}        # get AA016_a from AA016_a_run2_R1.fq.gz
       echo $file_base
       mkdir $file_base
       mv $file_base"_"* $file_base   # move there all BAM files
done
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


# 5.2) add read group headers, realign around indels and
#######################################################

# SCRIPTS 
abel_01_forloop_specify_samples_preprocessing_dirrange.sh  # for calling the second script that does the work
abel_02_GATK_joint_preprocessing_onedir.sh

# run it via nohup command (or screen)
nohup bash ~/aa/GATK/tetraploids_june_2016/abel_01_forloop_specify_samples_preprocessing_dirrange_bulk1.sh > nohup_AA_1.out&
#  ps -ef | grep 'filipko'


#  xxxxxxxxxxxxxxxxxxxx PART OF THE SCRIPT WITH THE REAL WORK  xxxxxxxxxxxxxxxxxxxx
BAMLIST=(*paired_run?.bam)
#First, merge the PE and SE files for each run and add the read group headers to each file
#Required headers are RGID, RGLB, RGPL, RGSM
echo "Merging paired and unpaired BAM files, adding RG headers..."
for ((i=0;i<${#BAMLIST[@]};++i)); do
        echo ${BAMLIST[$i]}
        ThisBam=${BAMLIST[$i]}                                  # AA016ac_paired_run2.bam
        BAM_file_base=${ThisBam%_paired_run?.bam}               # the file base should be e.g. AA016ac
	BAM_run_base=${ThisBam%.bam}                            # the run base should be e.g. AA016ac_paired_run2
        echo $BAM_file_base
        echo $BAM_run_base                                                              
        #merge ALL bam files (should be 3) with the same file base
        java -Xmx12g -jar $PICARD_HOME/MergeSamFiles.jar $(printf 'INPUT=%s ' $BAM_file_base*.bam) OUTPUT=$BAM_file_base.mergeRun.bam
 
        #extract run number for RGID assignment and add RG info headers
        ThisRunNumber=${BAM_run_base##*_}               #get the bit after the last "_" from AA016ac_paired_run2, e.g. run2
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

#   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



# 5.3) call variants using GATK tool Haplotype caller
######################################################

### IMPORTANT - FIRST SEPARATE DIPLOID AND TETRAPLOID SAMPLES
mv xxxx diploids/
mv xxxx tetraploids/

# rename folders to add sequential numbers at the end (xxx.1, xxx.2 etc) within each
n=1; for dir in AA*; do mv "$dir" "$(printf "$dir.$n")";   n=$((n+1)); done

### variant calling over all DIPLOID/TETRAPLOID samples separately, using GATK tool HaplotypeCaller
# as an ARRAY job
# SCRIPT abel_03_GATK_Hap_caller_submit_DIP.sh
#     submits array jobs of SCRIPT abel_03_GATK_Hap_caller_worker_DIP.sh

#  xxxxxxxxxxxxxxxxxxxx PART OF THE SCRIPT WITH THE REAL WORK  xxxxxxxxxxxxxxxxxxxx
bamfile=(*.realign.bam)
BAM_file_base=${bamfile%.realign.bam}

java -Xmx30g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T HaplotypeCaller \
  -R alygenomes.fasta \
  -I $BAM_file_base".realign.bam" \
  -ERC GVCF \
  -variant_index_type LINEAR \
  -variant_index_parameter 128000 \
  -ploidy 2 \
  -o $BAM_file_base".raw.g.vcf.gz"
#   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



# 5.4) joint genotyping over all samples, using GATK tool Genotype GVCFs
#####################################################################

# SCRIPT abel_04_GATK_genotypegvcfs_parallel.sh

#  xxxxxxxxxxxxxxxxxxxx PART OF THE SCRIPT WITH THE REAL WORK  xxxxxxxxxxxxxxxxxxxx
#first prepare list of all samples
ls *.raw.g.vcf.gz | sed 's/^/-V /' | sed 's/$/ /' > join.vcf.samplelist.txt
SAMPLELIST=$(<join.vcf.samplelist.txt)
echo "Joint genotyping samples $SAMPLELIST"

#then run Genotype GVCFs
# NB check that (i) -Xmx has cpus-per-task*mem-per-cpu (e.g. 16*3=48) and (ii) -nt = cpus-per-task
java -Xmx48g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T GenotypeGVCFs \
-R alygenomes.fasta \
$SAMPLELIST \
-nt 16 \
--includeNonVariantSites \
-o joint.raw.vcf.gz
#   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



################################################
#### 6)   FILTRATION
#################################################

# 6.1.)  first step RAW filtration using GATK (on server or locally)
#####################################################################

# SCRIPT abel_05_GATK_hardfilter_specifydir_covg_percindivs.sh
nohup bash /usit/abel/u1/filipko/aa/GATK/tetraploids_june_2016/abel_05_GATK_hardfilter_specifydir_covg_percindivs.sh jointgenotyping_tet_run12 8 0.5 > nohup_AA_filtering.out&


#  xxxxxxxxxxxxxxxxxxxx PART OF THE SCRIPT WITH THE REAL WORK  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# hardfilter SNPs from this joint file using recommended parameters for SNPS + exclude SNPs from set of putatively paralogou loci (from Harvard)
echo "Extracting non-paralogous SNPs from chr 1-8 joint vcf..."
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T SelectVariants \
   -R $REFSEQ \
   -V joint.raw.vcf.gz \
   -selectType SNP \
   -o joint.raw.SNP.vcf.gz \
   -L scaffold_1 -L scaffold_2 -L scaffold_3 -L scaffold_4 -L scaffold_5 -L scaffold_6 -L scaffold_7 -L scaffold_8 \
   --excludeIntervals $EXCLUDEPARALOGS \
   --excludeNonVariants

#create a filtering criterion
echo "Hard Filtering SNPs from joint vcf..."
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T VariantFiltration \
   -R $REFSEQ \
   -V joint.raw.SNP.vcf.gz \
   -o joint.hardfilt.SNP.vcf.gz \
   --filterExpression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5" \
   --filterName "AA_filter1"

# extract SNPs (both bi- and multiallelic) passing the AA_filter1 criterion
echo "Extracting only passing SNPs from joint vcf...."
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T SelectVariants \
   -R $REFSEQ \
   -V joint.hardfilt.SNP.vcf.gz \
   -o joint.hardfilt.SNP.PASS.vcf.gz \
   --excludeFiltered

# extract BIALLELIC SNPs
echo "Extracting only BIALLELIC SNPs from joint vcf...."
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T SelectVariants \
   -R $REFSEQ \
   -V joint.hardfilt.SNP.PASS.vcf.gz \
   -o joint.hardfilt.SNP.PASS.BIALL.vcf.gz \
   --excludeFiltered \
   --restrictAllelesTo BIALLELIC

# this tool works until 3.4 ... does not work on multiallelic SNPs produced by 3.5, just biallelic
# emit only sites in which x% of samples have sufficient coverage (e.g. 50% at least 8x coverage)
echo "Creating intervals for loci covered min $MinCoverage-times in defined $MinPercIndivs of indivs from joint vcf...."
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.3.0/GenomeAnalysisTK.jar -T CoveredByNSamplesSites \
   -R $REFSEQ \
   -V joint.hardfilt.SNP.PASS.BIALL.vcf.gz \
   -out dp$MinCoverage.perc$MinPercIndivs.BIALL.intervals \
   --minCoverage $MinCoverage \
   --percentageOfSamples $MinPercIndivs

# select variants based on this interval list (NB variants with < defined coverage will be still present in VCF)
echo "Selecting variants based on presence in $MinPercIndivs of indivs"
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T SelectVariants \
   -R $REFSEQ \
   -V joint.hardfilt.SNP.PASS.BIALL.vcf.gz \
   -o joint.hardfilt.SNP.PASS.BIALL.dp$MinCoverage.perc$MinPercIndivs.vcf.gz \
   -L dp$MinCoverage.perc$MinPercIndivs.BIALL.intervals

#  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx




# 6.2) copy *vcf.gz files and then do other filtering and input file conversions locally using GATK/bcftools and PGD Spider (PGDSpider only for diploids)
#########################################################################################################################################################

# for this I need a (i) subfolder called "PGDSpider_out" and inside PGD spider conversion specification files (*spid) and popfile ("popfile_pops_PGDSpider.txt" = first column = individual names, second column = populations/groups)
# SCRIPT GATK_vcftools_filterexport_definegroup.sh


#####################################################
# conversions including tetraploid format - locally
#####################################################

# XXX https://github.com/samtools/bcftools/wiki/HOWTOs

# subsetting using a list of samples (just one sample name per row), excluding nonvariable sites within that subset
bcftools view joint.hardfilt.SNP.PASS.BIALL.dp8.perc0.5.vcf.gz -S seznam_Tatry --min-af 0.0001:nonmajor -Oz > tatry.vcf.gz
bcftools view joint.hardfilt.SNP.PASS.BIALL.dp8.perc0.5.vcf.gz -S tatry_2x --min-af 0.0001:nonmajor -Oz > tatry_2x.vcf.gz
bcftools view joint.hardfilt.SNP.PASS.BIALL.dp8.perc0.5.vcf.gz -S tatry_4x --min-af 0.0001:nonmajor -Oz > tatry_4x.vcf.gz
bcftools view joint.hardfilt.SNP.PASS.BIALL.dp8.perc0.5.vcf.gz -S all_4x --min-af 0.0001:nonmajor -Oz > all_4x.vcf.gz

bcftools view joint.hardfilt.SNP.PASS.BIALL.dp8.perc0.5.vcf.gz -S tatry_2xA --min-af 0.0001:nonmajor -Oz > tatry_2x.vcf.gz
bcftools view joint.hardfilt.SNP.PASS.BIALL.dp8.perc0.5.vcf.gz -S tatry_2xF --min-af 0.0001:nonmajor -Oz > tatry_2x.vcf.gz
bcftools view joint.hardfilt.SNP.PASS.BIALL.dp8.perc0.5.vcf.gz -S tatry_2xLA --min-af 0.0001:nonmajor -Oz > tatry_2x.vcf.gz



# get the same sites from a genome-wide VCF and prepare an intervals list for GATK
bcftools query -f'%CHROM %POS\n' joint.hardfilt.SNP.PASS.BIALL.dp8.perc0.5.vcf.gz | sed 's/ /:/g' - > sites_joint.hardfilt.SNP.PASS.BIALL.dp8.perc0.5.intervals
bcftools query -f'%CHROM %POS\n' tatry.vcf.gz | sed 's/ /:/g' - > sites_tatry.intervals
bcftools query -f'%CHROM %POS\n' all_4x.vcf.gz | sed 's/ /:/g' - > sites_all_4x.intervals

# select variants based on this list
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -T SelectVariants -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -V ~/genomic/300/nomissing_F5/scaf_1_PASS.vcf.gz -o all_fromgenomic.vcf.gz -L sites_joint.hardfilt.SNP.PASS.BIALL.dp8.perc0.5.intervals 




### STRUCTURE
# export a loci object with pegas and then further convert using unix
library(pegas)
aa.loci <- read.vcf("all_4x.vcf.gz", which.loci = 1:74149) # must be exact number otherwise only 10,000 loci read
write.loci(aa.loci, file = "all_4x_onerowperind_raw.stru", loci.sep ="\t", quote = FALSE, allele.sep ="\t", col.names = FALSE) # write structure input for ONEROWPERIND=1 option, takes some time, each indiv. is one row, alleles separated by tab 
# in bash: Insert population number (just "1") into the second column and change letters to integers and missing to -9
gawk '{first = $1; $1 = "1"; print first, $0; }' all_4x_onerowperind_raw.stru | sed 's/ A/ 1/g' - | sed 's/ C/ 2/g' - | sed 's/ G/ 3/g' - | sed 's/ T/ 4/g' - | sed 's/ \./ -9/g' - > all_4x_onerowperind.stru


#### TREEMIX
# first split original COMPLETE VCF by populations/groups
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf tatry_2xA -o tatry_2xA_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf tatry_2xF -o tatry_2xF_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf tatry_2xLA -o tatry_2xLA_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf tatry_4x -o tatry_4x_split.vcf.gz

# move all to treemix folder
mv *split.vcf.gz* ./treemix

#then, get allele counts table using a forloop

for file in *vcf.gz; do java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T VariantsToTable -V $file -F CHROM -F POS -F AC -F AN -F DP -raw --allowMissingData -o $file.table; done

# and create the file
awk '{OFS=","};{print $4-$3, $3}' tatry_2xA_split.vcf.gz.table > file1
awk '{OFS=","};{print $4-$3, $3}' tatry_2xF_split.vcf.gz.table > file2
awk '{OFS=","};{print $4-$3, $3}' tatry_2xLA_split.vcf.gz.table > file3
awk '{OFS=","};{print $4-$3, $3}' tatry_4x_split.vcf.gz.table > file4

printf "tatry_2xA %s\t%s tatry_2xF %s\t%s tatry_2xLA %s\t%s tatry_4x %s\n%s"  > tatry.frq
paste file1 file2 file3 file4 | grep -v 'AC' >> tatry.frq
gzip tatry.frq



# -------------------------------------------------------------

### build ML tree of the populations
# -root = set the outgroup, e.g. named tatry_2xF
# -k = set how many successive SNPs are in LD (here I assume every 100 th, thus setting -k 100)
# -m allow for m migration events
treemix -i tatry.frq.gz -root tatry_2xF -k 100 -m 1 -o outstem

### bootstrap over blocks of -k SNPs ... this generates one replicate
treemix -i tatry.frq.gz -bootstrap -k 100 -o replicate1

### bootstrap over 100 replicates
# creates a file BOOTSTRAP.tre in the bootstrapping folder, could be opened e.g. in Figtree
# package DendroPy needs to be installed
mkdir tatry_bootstraps
cp tatry.frq.gz tatry_bootstraps
cd tatry_bootstraps
for i in {1..100}; do treemix -i tatry.frq.gz -bootstrap -k 100 -o replicate$i; done
gunzip replic*treeout.gz
sumtrees.py --decimals=0 --percentages --output-tree-filepath=BOOTSTRAP.tre --target=replicate1.treeout replicate2.treeout replicate3.treeout replicate4.treeout replicate5.treeout replicate6.treeout replicate7.treeout replicate8.treeout replicate9.treeout replicate10.treeout replicate11.treeout replicate12.treeout replicate13.treeout replicate14.treeout replicate15.treeout replicate16.treeout replicate17.treeout replicate18.treeout replicate19.treeout replicate20.treeout replicate21.treeout replicate22.treeout replicate23.treeout replicate24.treeout replicate25.treeout replicate26.treeout replicate27.treeout replicate28.treeout replicate29.treeout replicate30.treeout replicate31.treeout replicate32.treeout replicate33.treeout replicate34.treeout replicate35.treeout replicate36.treeout replicate37.treeout replicate38.treeout replicate39.treeout replicate40.treeout replicate41.treeout replicate42.treeout replicate43.treeout replicate44.treeout replicate45.treeout replicate46.treeout replicate47.treeout replicate48.treeout replicate49.treeout replicate50.treeout replicate51.treeout replicate52.treeout replicate53.treeout replicate54.treeout replicate55.treeout replicate56.treeout replicate57.treeout replicate58.treeout replicate59.treeout replicate60.treeout replicate61.treeout replicate62.treeout replicate63.treeout replicate64.treeout replicate65.treeout replicate66.treeout replicate67.treeout replicate68.treeout replicate69.treeout replicate70.treeout replicate71.treeout replicate72.treeout replicate73.treeout replicate74.treeout replicate75.treeout replicate76.treeout replicate77.treeout replicate78.treeout replicate79.treeout replicate80.treeout replicate81.treeout replicate82.treeout replicate83.treeout replicate84.treeout replicate85.treeout replicate86.treeout replicate87.treeout replicate88.treeout replicate89.treeout replicate90.treeout replicate91.treeout replicate92.treeout replicate93.treeout replicate94.treeout replicate95.treeout replicate96.treeout replicate97.treeout replicate98.treeout replicate99.treeout replicate100.treeout

### PLOTTING - open R and type following in R
# the whole "/src" dir from the original treemix tarball is necessary to copy to the workingdir
# 1) set working dir to your treemix files location
source("src/plotting_funcs.R")        # source an R script "plotting_funcs.R" 
plot_tree("outstem")                  # draw the plot 
plot_resid("outstem", "pops.poplist.txt") # plot residuals, must contain list of pops in special file



########  for pops ############################
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA016.txt -o AA016_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA023.txt -o AA023_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA096.txt -o AA096_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA162.txt -o AA162_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA165.txt -o AA165_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA167.txt -o AA167_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA178.txt -o AA178_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA183.txt -o AA183_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA208.txt -o AA208_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA225.txt -o AA225_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA227.txt -o AA227_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA229.txt -o AA229_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA234.txt -o AA234_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA235.txt -o AA235_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA236.txt -o AA236_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA243.txt -o AA243_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA248.txt -o AA248_split.vcf.gz
java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T SelectVariants -V tatry.vcf.gz -sf AA321.txt -o AA321_split.vcf.gz


# move all to treemix folder
mkdir ./treemix_pops
mv *split.vcf.gz* ./treemix_pops

#then, get allele counts table using a forloop
cd treemix_pops
for file in *vcf.gz; do java -Xmx2g -jar ~/programs/GATK/3.5/GenomeAnalysisTK.jar -R /home/filip/genomic/arenosa_data/JIC_reference/alygenomes.fasta -T VariantsToTable -V $file -F CHROM -F POS -F AC -F AN -F DP -raw --allowMissingData -o $file.table; done

# and create the file
awk '{OFS=","};{print $4-$3, $3}' AA016_split.vcf.gz.table > file1
awk '{OFS=","};{print $4-$3, $3}' AA023_split.vcf.gz.table > file2
awk '{OFS=","};{print $4-$3, $3}' AA096_split.vcf.gz.table > file3
awk '{OFS=","};{print $4-$3, $3}' AA162_split.vcf.gz.table > file4
awk '{OFS=","};{print $4-$3, $3}' AA165_split.vcf.gz.table > file5
awk '{OFS=","};{print $4-$3, $3}' AA167_split.vcf.gz.table > file6
awk '{OFS=","};{print $4-$3, $3}' AA178_split.vcf.gz.table > file7
awk '{OFS=","};{print $4-$3, $3}' AA183_split.vcf.gz.table > file8
awk '{OFS=","};{print $4-$3, $3}' AA208_split.vcf.gz.table > file9
awk '{OFS=","};{print $4-$3, $3}' AA225_split.vcf.gz.table > file10
awk '{OFS=","};{print $4-$3, $3}' AA227_split.vcf.gz.table > file11
awk '{OFS=","};{print $4-$3, $3}' AA229_split.vcf.gz.table > file12
awk '{OFS=","};{print $4-$3, $3}' AA234_split.vcf.gz.table > file13
awk '{OFS=","};{print $4-$3, $3}' AA235_split.vcf.gz.table > file14
awk '{OFS=","};{print $4-$3, $3}' AA236_split.vcf.gz.table > file15
awk '{OFS=","};{print $4-$3, $3}' AA243_split.vcf.gz.table > file16
awk '{OFS=","};{print $4-$3, $3}' AA248_split.vcf.gz.table > file17
awk '{OFS=","};{print $4-$3, $3}' AA321_split.vcf.gz.table > file18


printf "AA016 %s\t%s AA023 %s\t%s AA096 %s\t%s AA162 %s\t%s AA165 %s\t%s AA167 %s\t%s AA178 %s\t%s AA183 %s\t%s AA208 %s\t%s AA225 %s\t%s AA227 %s\t%s AA229 %s\t%s AA234 %s\t%s AA235 %s\t%s AA236 %s\t%s AA243 %s\t%s AA248 %s\t%s AA321 %s\n%s"  > tatry.frq
paste file1 file2 file3 file4 file5 file6 file7 file8 file9 file10 file11 file12 file13 file14 file15 file16 file17 file18 | grep -v 'AC' >> tatry.frq
gzip tatry.frq
######## for pops ############################

treemix -i tatry.frq.gz -root AA023 -k 100 -m 1 -o outstem







# TODO overlay resulting VCF with genomic ... get Tatrean SNPs from genomic data
# TODO get brian's RAD data














###################################################################################################################
###################################################################################################################
###################################################################################################################

####### leftovers TODO 
# coverage analysis pea all samples
#create list of BAMs, that were used for variant calling (they are in the same dir as the script)
ls ../*.merged.bam | sed 's/^/-I /' | sed 's/$/ /' > merged.bamlist.txt
SAMPLELIST=$(<merged.bamlist.txt)
echo "Joint genotyping samples $SAMPLELIST"
#provide this bamlist to coverage walker
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T FindCoveredIntervals \
   -R $REFSEQ \
   $SAMPLELIST \
   -cov $MinCoverage \
   -o covered.intervals.list

nohup bash ~/aa/GATK/tetraploids_apr_2016/MMM_MODIFY_FOR_3.5abel_05_GATK_hardfilter_specifydir_covg_percindivs.sh jointgenotyping_tet_run1 8 0.5 > nohup_bams.out&

#create list of intervals (one per BAM) with regions ABOVE defined MinCoverage threshold 
# based on BAMs, that were used for variant calling (they are in the same dir as the script)

#cd ../

for bamfile in *.merged.bam; do
java -Xmx12g -jar /usit/abel/u1/filipko/programs/GATK/3.5/GenomeAnalysisTK.jar -T FindCoveredIntervals \
   -R $REFSEQ \
   -I $bamfile \
   -cov $MinCoverage \
   -o $bamfile.covered.intervals.list;
done
