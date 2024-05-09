library(dplyr)
library(reshape)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(dplyr)

cptac3_cases=read_delim('../../out.bam.paths/Cases_paths_katmai_storage1.20221108.tsv',delim='\t')
cptac3_cases=as.data.frame(cptac3_cases)
cptac2_cases=read_delim('../../out.bam.paths/CPTAC2_Cases_paths_katmai.20221108.tsv',delim='\t')
cptac2_cases=as.data.frame(cptac2_cases)

all_cases=unique(c(cptac2_cases$case,cptac3_cases$case))


samples_tab=read.table('GO_case_list_CPTAC_germline.tsv',sep='\t',header=T)
samples_tab$Cancer_type=gsub('HGSC','OV',samples_tab$Cancer_type)
samples_tab$Cancer_type=gsub('ccRCC','CCRCC',samples_tab$Cancer_type)

#length(intersect(samples_tab$Case_ID,all_cases))
#1,062 -- for 2 samples we don't have BAM-files

all_cases_2=intersect(samples_tab$Case_ID,all_cases)
#samples_tab[!(samples_tab$Case_ID %in% all_cases_2),]
#      	Case_ID 	Cancer_type
#	03BR011         BRCA
#       C3N-02788       GBM


####################################################################
#2023-07-18: get the numbers for T/N samples with RNA-seq BAM-files:
cptac2_bam_f=cptac2_cases[cptac2_cases$case %in% all_cases_2,]
cptac3_bam_f=cptac3_cases[cptac3_cases$case %in% all_cases_2,]

#Check if you have here a single BAM-file listed for each sample ID:
#nrow(cptac2_bam_f[!duplicated(cptac2_bam_f$New_ID),])==nrow(cptac2_bam_f)
#nrow(cptac3_bam_f[!duplicated(cptac3_bam_f$New_ID),])==nrow(cptac3_bam_f)
#TRUE for both
cptac2_bam_f$Sample_type=gsub('.*-([T/N])','\\1',cptac2_bam_f$New_ID)
cptac3_bam_f$Sample_type=gsub('.*-([T/N])','\\1',cptac3_bam_f$New_ID)

#Get numbers for Tumor and NAT samples used in this analysis:
table(cptac2_bam_f$Sample_type)
table(cptac2_bam_f$Sample_type)

#This piece of code is needed for Methods only
####################################################################




#Now read table of variants that were detected in proteome:
tab=read.table('Vars_DETECTED_in_either_assay.ALL.10cancers.SAAV_SpliceForms.20230605.tsv',sep='\t',header=T)

tab_s=tab[tab$Sample %in% all_cases_2,]
tab_s$chr=gsub('(.*)\\:[A-Z]+([0-9]+)[A-Z]+','\\1',tab_s$Variant)
tab_s$new_chr=paste('chr',tab_s$chr,sep='')

tab_s$Coordinate=gsub('.*\\:[A-Z]+([0-9]+)[A-Z]+','\\1',tab_s$Variant)
tab_s$Reference=gsub('.*\\:([A-Z]+)([0-9]+)([A-Z]+)','\\1',tab_s$Variant)
tab_s$Alternate=gsub('.*\\:([A-Z]+)([0-9]+)([A-Z]+)','\\3',tab_s$Variant)

#Remove X/Y chromosomes:
tab_s1=tab_s[!(tab_s$new_chr %in% c('chrX','chrY')),]

#Also remove non-SNPs (only 1 to 1 substitutions should remain after that)
tab_s1$Ref_len=str_length(tab_s1$Reference)
tab_s1$Alt_len=str_length(tab_s1$Alternate)

#Remove splice-variants:
tab_s2=tab_s1[tab_s1$Ref_len==1 & tab_s1$Alt_len==1,]

colnames(tab_s2)<-gsub('new_chr','Chromosome',colnames(tab_s2))
colnames(tab_s2)<-gsub('Coordinate','Start',colnames(tab_s2))

#limit them to cancer genes only:
c_genes=read.table('624_Cancer_genes.txt',sep='\t',header=F)
tab_s3=tab_s2[tab_s2$geneSymbol %in% c_genes$V1,]

write.table(tab_s3, 'out/Filtered_vars_for_ASE_test_input.20230702.tsv',sep='\t',quote=F,row.names=F)

#Keep only unique variants across all datasets:
tab_s4=tab_s3[,c('Chromosome','Start','Sample','Reference','Alternate','Variant','Cancer')]
tab_s4=tab_s4[!duplicated(tab_s4),]

dir.create('Site_lists')
write.table(tab_s4,"Site_lists/All_sites.tsv", sep='\t',quote=F,row.names=F)
