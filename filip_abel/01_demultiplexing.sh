for dir in raw*; do
	echo "processing $dir"
	cd $dir
	gunzip *.gz
        cat *1_sequence.txt | fastx_barcode_splitter.pl --bcfile barcodes_AA_tet_run2_$dir.txt --prefix ../samples/ --suffix "_R1.fq" --bol --mismatches 1 
        cat *2_sequence.txt | fastx_barcode_splitter.pl --bcfile barcodes_AA_tet_run2_$dir.txt --prefix ../samples/ --suffix "_R2.fq" --bol --mismatches 1 
	cd ..
done

