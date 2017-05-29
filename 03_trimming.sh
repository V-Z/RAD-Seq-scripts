# I have each fq pair in specific dir named by the sample name
# remove first 8 bases at the beginning of each read, keep reads using default threshold 20 QC score, remove reads shorter than 100 bp, trim some low quality ends
for dir in AA*; do
	echo "processing $dir"
        mkdir ~/aa/GATK/tetraploids_june_2016/trimmed/$dir   # change this path if necessary
	cp ~/programs/trimmomatic/adapters/TruSeq3-PE.fa ./$dir
	cd $dir
        java -jar ~/programs/trimmomatic/trimmomatic-0.36.jar PE -phred33 $dir"_run2_R1.fq.gz" $dir"_run2_R2.fq.gz" ~/aa/GATK/tetraploids_june_2016/trimmed/$dir/$dir"_trm_run2_R1.fq.gz" ~/aa/GATK/tetraploids_june_2016/trimmed/$dir/$dir"_unpaired_run2_R1.fq.gz" ~/aa/GATK/tetraploids_june_2016/trimmed/$dir/$dir"_trm_run2_R2.fq.gz" ~/aa/GATK/tetraploids_june_2016/trimmed/$dir/$dir"_unpaired_run2_R2.fq.gz" HEADCROP:8 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:100  # change the paths if necessary	
	cd ..
done


