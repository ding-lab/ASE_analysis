#2023-06-28: also copy PDA: scp -r nvterekhanova@compute1-client-1.ris.wustl.edu:/storage1/fs1/dinglab/Active/Projects/PanCan_Germline_CPTAC/Analysis/Analysis_nadya/Analysis/ASE_analysis/readCount_analysis/20230622_all_detected_vars_upd/Tumor/PDA/ .


for dir in LUAD LSCC HNSCC GBM CCRCC UCEC PDA BRCA CO OV

do 
    for f in Tumor/$dir/*; do perl Parse_results.pl $f; done;

done


for dir in LUAD LSCC HNSCC GBM CCRCC UCEC

do 
    for f in NAT/$dir/*; do perl Parse_results.pl $f; done;

done

