# loop only over directories from all samples in range AA001* to AA400*
# nb works in bash-4 onwards
for dir in AL{001..400}* ; do
	if [[ -d $dir ]]; then
		echo "processing $dict"
		bash abel_02_GATK_joint_preprocessing_onedir.sh $dir    
	fi	
done




