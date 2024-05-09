library(dplyr)
library(reshape)
library(tidyverse)
library(ggrepel)
library(dplyr)


#2023-06-21: lets use >=10 for coverage cutoff
#In the previous paper a more loose cutoff was used: "We retained variants with at least six read counts (reference allele + alternate allele ≥ 6) for ASE analysis."
#https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00964-1



#Need to filter by the SAAV that were detected in either of the assays in a specific cancer type:
all_vars=read.table('Vars_DETECTED_in_either_assay.ALL.10cancers.SAAV_SpliceForms.20230605.tsv',sep='\t',header=T)
all_vars$SAAV=gsub('.*-(.*)','\\1',all_vars$ID)
all_vars_1=all_vars[,c('geneSymbol','Variant','Cancer','SAAV')]

#Kepp a single ident per cancer and across all datasets:
all_vars_1=all_vars_1[!duplicated(all_vars_1),]



#First, do this for Tumor-samples.
samples_tab=read.table('GO_case_list_CPTAC_germline.tsv',sep='\t',header=T)
samples_tab$Cancer_type=gsub('HGSC','OV',samples_tab$Cancer_type)
samples_tab$Cancer_type=gsub('ccRCC','CCRCC',samples_tab$Cancer_type)

cancers=c('LUAD','LSCC','UCEC','CCRCC','HNSCC','GBM','BRCA','CO','OV','PDA')
all_read_c=NULL
for (can in cancers){
	read_c=read.table(paste('DATA/BamReadCount_added/Tumor/',can,'.BamReadCount',sep=''),sep='\t',header=F)
	samples_tab_c=samples_tab[samples_tab$Cancer==can,]
	colnames(read_c)=c('Chromosome','Start','Sample','Reference','Alternate','Variant','Cancer','Ref_count','Alt_count')
	read_c1=read_c[read_c$Sample %in% samples_tab_c$Case_ID,]
	all_read_c=rbind(all_read_c,read_c1)
}

#Filter further to those that were detected in the proteome in a particular cancer type:
all_read_c$Variant_Cancer=paste(all_read_c$Variant,all_read_c$Cancer,sep='__')
all_vars_1$Variant_Cancer=paste(all_vars_1$Variant,all_vars_1$Cancer,sep='__')

all_read_c2=all_read_c[all_read_c$Variant_Cancer %in% all_vars_1$Variant_Cancer,]

write.table(all_read_c2,'out/Read_count_forALL_germline_variants_fromDBs.20230702.tsv',sep='\t',row.names=F)


tab_c=all_read_c2	
tab_c$Ref_count[is.na(tab_c$Ref_count)]<-0
tab_c$Alt_count[is.na(tab_c$Alt_count)]<-0
tab_c$Coverage=tab_c$Ref_count+tab_c$Alt_count
tab_c=tab_c[tab_c$Coverage>=10,]
#N=52,876 with coverage at least 10

#Read GTs with Ref and Alt:
gts=read.table('DATA/GTs_for_SAAV_sites.20230702.tsv',sep='\t',header=T)

#table(gts$GT)
gts_g=gts[gts$GT %in% c('0/1','1/0'),]
gts_g=gts_g[,c('Chromosome','Start','Sample','Reference', 'Alternate', 'GT')]


tab_c2=merge(tab_c,gts_g)

#Perform a two-sided binomial test with a null probability of success 0.5 in a Bernoulli experiment to identify ASE (https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00964-1)
all_str=NULL
for (i in 1:nrow(tab_c2)){
	str=tab_c2[i,]
	st=binom.test(str$Alt_count, str$Coverage, p = 0.5, alternative='two.sided')
	str$p_val=st$p.value
	all_str=rbind(all_str,str)
}

all_str$FDR=p.adjust(all_str$p_val,method='fdr')

all_str$Status='None'
all_str$Status=ifelse(all_str$FDR<0.05,'Significant',all_str$Status)
all_str$VAF=all_str$Alt_count/(all_str$Alt_count+all_str$Ref_count)
all_str$Status_1=all_str$Status
all_str$Status_1=ifelse(all_str$FDR<0.05 & all_str$VAF<0.5,'Pref_REF',all_str$Status_1)
all_str$Status_1=ifelse(all_str$FDR<0.05 & all_str$VAF>0.5,'Pref_ALT',all_str$Status_1)


write.table(all_str,'out/ASE_test_results.10cancers.20230702.tsv',sep='\t',quote=F,row.names=F)

#check that there are no duplicates (we have removed them later and updated FDRs).


###############################
###Also do the same for NAT:###
###############################
samples_tab=read.table('GO_case_list_CPTAC_germline.tsv',sep='\t',header=T)
samples_tab$Cancer_type=gsub('HGSC','OV',samples_tab$Cancer_type)
samples_tab$Cancer_type=gsub('ccRCC','CCRCC',samples_tab$Cancer_type)

cancers=c('LUAD','LSCC','UCEC','CCRCC','HNSCC')
all_read_c=NULL
for (can in cancers){
	read_c=read.table(paste('DATA/BamReadCount_added/NAT/',can,'.BamReadCount',sep=''),sep='\t',header=F)
	samples_tab_c=samples_tab[samples_tab$Cancer==can,]
	colnames(read_c)=c('Chromosome','Start','Sample','Reference','Alternate','Variant','Cancer','Ref_count','Alt_count')
	read_c1=read_c[read_c$Sample %in% samples_tab_c$Case_ID,]
	all_read_c=rbind(all_read_c,read_c1)
}

#Some are NAs, as we didn't have BAM-files for all the samples:
#but we will filter based on coverage later.


#Filter further to those that were detected in the proteome (we had a larger list of sites than needed):
all_read_c$Variant_Cancer=paste(all_read_c$Variant,all_read_c$Cancer,sep='__')
all_vars_1$Variant_Cancer=paste(all_vars_1$Variant,all_vars_1$Cancer,sep='__')

all_read_c2=all_read_c[all_read_c$Variant_Cancer %in% all_vars_1$Variant_Cancer,]

write.table(all_read_c2,'out/Read_count_forALL_germline_variants_fromDBs.NAT.20230702.tsv',sep='\t',row.names=F)

tab_c=all_read_c2	
tab_c$Ref_count[is.na(tab_c$Ref_count)]<-0
tab_c$Alt_count[is.na(tab_c$Alt_count)]<-0
tab_c$Coverage=tab_c$Ref_count+tab_c$Alt_count
tab_c=tab_c[tab_c$Coverage>=10,]
#N=19,061

#Read GTs with Ref and Alt:
gts=read.table('DATA/GTs_for_SAAV_sites.20230702.tsv',sep='\t',header=T)

#table(gts$GT)
gts_g=gts[gts$GT %in% c('0/1','1/0'),]
gts_g=gts_g[,c('Chromosome','Start','Sample','Reference', 'Alternate', 'GT')]


tab_c2=merge(tab_c,gts_g)


#Perform  a two-sided binomial test with a null probability of success 0.5 in a Bernoulli experiment to identify ASE (https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00964-1)
all_str=NULL
for (i in 1:nrow(tab_c2)){
	str=tab_c2[i,]
	st=binom.test(str$Alt_count, str$Coverage, p = 0.5, alternative='two.sided')
	str$p_val=st$p.value
	all_str=rbind(all_str,str)
}

all_str$FDR=p.adjust(all_str$p_val,method='fdr')

all_str$Status='None'
all_str$Status=ifelse(all_str$FDR<0.05,'Significant',all_str$Status)
all_str$VAF=all_str$Alt_count/(all_str$Alt_count+all_str$Ref_count)
all_str$Status_1=all_str$Status
all_str$Status_1=ifelse(all_str$FDR<0.05 & all_str$VAF<0.5,'Pref_REF',all_str$Status_1)
all_str$Status_1=ifelse(all_str$FDR<0.05 & all_str$VAF>0.5,'Pref_ALT',all_str$Status_1)



write.table(all_str,'out/ASE_test_results.5cancers.NAT.20230702.tsv',sep='\t',quote=F,row.names=F)

#check that there are no duplicates (we have removed them later and updated FDRs).
