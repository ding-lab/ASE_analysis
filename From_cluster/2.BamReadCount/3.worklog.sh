#Also see here for details on bam-readcount output.

#Tumor
type='Tumor'

for dir in CO OV BRCA GBM HNSCC UCEC LSCC LUAD CCRCC PDA
do
perl Add_Bam_readCount.20230621.pl $dir $type
done

#NO Ref inconsistencies found




#NAT
type='NAT'

for dir in HNSCC UCEC LSCC LUAD CCRCC
do
perl Add_Bam_readCount.20230621.pl $dir $type
done

#NO Ref inconsistencies found
