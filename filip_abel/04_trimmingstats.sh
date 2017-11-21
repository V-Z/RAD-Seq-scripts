#calculate stats (N all reads)
touch trimming_info.csv
printf "indiv %s\t%s R1_trm_i.e.N_read_PAIRS %s\t%s R1_unpaired_reads %s\t%s R2_unpaired_reads %s\n%s"  >> trimming_info.csv
for dir in AA*
do
echo "processing $dir"
printf '%s\t%s' $dir >> trimming_info.csv  # print sample name
#echo $file |xargs printf '%s\t%s' >> mapping_info.csv   # print file name
less $dir/$dir"_trm_run2_R1.fq.gz" | echo $((`wc -l`/4)) |xargs printf '%s\t%s' >> trimming_info.csv           # No of R1 reads (N of lines /4)
less $dir/$dir"_unpaired_run2_R1.fq.gz" | echo $((`wc -l`/4)) |xargs printf '%s\t%s' >> trimming_info.csv      # No of R1 unpaired reads (N of lines /4)
less $dir/$dir"_unpaired_run2_R2.fq.gz" | echo $((`wc -l`/4)) |xargs printf '%s\t%s' >> trimming_info.csv      # No of R2 unpaired reads (N of lines /4)
printf '%s\n%s' >> trimming_info.csv
done

