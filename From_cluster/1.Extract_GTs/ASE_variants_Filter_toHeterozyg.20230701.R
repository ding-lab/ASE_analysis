#export TMPDIR=/diskmnt/Projects/Users/$USER/mytmp
library(reshape)
library(tidyverse)
library(dplyr)


tab=read.table(paste('/diskmnt/Projects/Users/nvterekhanova/Germline_analysis/15.ASE_analysis/20230622_all_detected_vars_upd/',
'Site_lists/All_sites.tsv',sep=''),sep='\t',header=T)

#Remove splice-site variants (from INDELs):
#tab_s=tab[!grepl('^ENSP',tab$Start),]
#No splice-sites in the latest list.

gt=read.table('../20230701_Extract_GTs/out/All_germline_var_sites.20230701.tsv',sep='\t',header=F)
colnames(gt)=c('Chromosome', 'Start','Reference', 'Alternate','Other','Sample','GT')
gt$Sample=gsub('(.*)\\.N','\\1',gt$Sample)
gt$ID=paste(gt$Chromosome,gt$Start,gt$Sample,sep='_')


tab$ID=paste(tab$Chromosome,tab$Start,tab$Sample,sep='_')

#Filter to sites that are in the Site_lists file:
#gt_1=gt[gt$ID %in% tab$ID,]


gt_1=gt[,c('ID','GT','Reference','Alternate')]
colnames(gt_1)[3:4]=c('VCF_Reference','VCF_Alternate')

tab_1=merge(tab,gt_1)


tab_2=tab_1[tab_1$Reference==tab_1$VCF_Reference & tab_1$Alternate==tab_1$VCF_Alternate,]
#nrow(tab_2)==nrow(tab) --> TRUE; all sites were found in the VCF file!


#Save the table with only matching Ref and Alt
write.table(tab_2,'out/GTs_for_SAAV_sites.20230702.tsv',sep='\t',row.names=F,quote=F)
