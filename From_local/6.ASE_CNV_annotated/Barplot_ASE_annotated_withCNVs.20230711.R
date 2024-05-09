library(dplyr)
library(reshape)
library(reshape2)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)

theme_nogrid = function(...) theme_bw() + theme(axis.line = element_line(colour = "black"),
                            panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank())



tum=read.table('ASE_test_results.10cancers.20230702.tsv',sep='\t',header=T)

norm=read.table('ASE_test_results.5cancers.NAT.20230702.tsv',sep='\t',header=T)

cancers=c('BRCA','OV','CO','GBM','LUAD','CCRCC','UCEC', 'HNSCC','LSCC','PDA')

#reaad cancer genes
c_genes=read.table('624_Cancer_genes.txt',sep='\t',header=F)

tum_s=tum[,c('Variant','Status_1','Sample','Cancer','VAF')]
norm_s=norm[,c('Variant','Status_1','Sample','Cancer','VAF')]
tum_s$Type='Tumor'
norm_s$Type='NAT'

both_ase=rbind(tum_s,norm_s)

#Read table with selected variants for summary DotPlot:
sel_snps=read.table('For_summary_DotPlot_data.v9.Upd.Annotation.20230725.tsv',sep='\t',header=T)

sel_snps_s=sel_snps[,c('geneSymbol','Variant','New_ID_2')]
sel_snps_s=sel_snps_s[!duplicated(sel_snps_s),]

both_ase_2=both_ase[both_ase$Variant %in% sel_snps_s$Variant,]
both_ase_2=merge(both_ase_2,sel_snps_s,all.x=T)
colnames(both_ase_2)[c(3,7)]=c('Case','Gene')

#Filter to tumor ASEs:
both_ase_2=both_ase_2[both_ase_2$Type=='Tumor',]

cnvs_tab=read.table('GISTIC_results.thresholded.formatted.1064Cases.20230711.tsv',sep='\t',header=T)

both_ase_3=merge(both_ase_2,cnvs_tab,all.x=T)

both_ase_3$CNV='Neutral'
both_ase_3$CNV=ifelse(both_ase_3$CNV_status>0,'AMP',both_ase_3$CNV) 
both_ase_3$CNV=ifelse(both_ase_3$CNV_status<0,'DEL',both_ase_3$CNV) 



both_ase_4=both_ase_3[both_ase_3$Status_1!='None',]
st=as.data.frame(table(both_ase_4[,c('CNV','Gene','New_ID_2')]))
st1=as.data.frame(table(both_ase_4[,c('Gene','New_ID_2')]))
colnames(st)<-gsub('Freq','Count',colnames(st))
colnames(st1)<-gsub('Freq','Total',colnames(st1))

st_all=merge(st,st1)
st_all$Fraction=st_all$Count/st_all$Total
st_all=st_all[st_all$Total!=0,]
st_all=st_all[order(st_all$Total),]

colnames(st_all)<-gsub('New_ID_2','ID',colnames(st_all))
st_all$ID=factor(st_all$ID,levels=unique(st_all$ID))

st_all_v2=st_all[order(-st_all$Total),]
sel_ids=unique(as.character(unlist(st_all_v2$ID)))[1:25]

to_plot=st_all[st_all$ID %in% sel_ids,]

p <- ggplot(to_plot, aes(x = factor(ID), y = Count)) 

p <- p + geom_bar(stat="identity",aes(fill = CNV),color=NA,width=0.9) 
  
p <- p + theme(axis.text.x = element_text(colour="black", size=10, angle=45, vjust = 1,hjust=1), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank()) + labs(x="") 

p <- p + scale_fill_manual(values=c('Neutral'='grey','AMP'='#E41A1C','DEL'='#377eb8'))

p <- p + theme_minimal()

p <- p + theme(axis.title = element_blank(), axis.text.x = element_text(colour="black", size=10,angle=90), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())

p <- p + theme_classic()+coord_flip()
            
pdf(paste("Barplot_CNV_summary.LatestList.top25.20230725.pdf",sep=""), width=7, height=4,useDingbats=FALSE)
p
dev.off()



#The code below is to create CNV-table:
########################
###Look at CNV-levels###
########################
#Use CNV results from Broad pipeline (GISTIC-based).

samples_tab=read.table('GO_case_list_CPTAC_germline.tsv',sep='\t',header=T)
samples_tab$Cancer_type=gsub('HGSC','OV',samples_tab$Cancer_type)
samples_tab$Cancer_type=gsub('ccRCC','CCRCC',samples_tab$Cancer_type)

cnvs=read_delim('all_thresholded.by_genes.txt',delim='\t')
cnvs=as.data.frame(cnvs)
gist_ids=read.table('GISTIC_Matched_Samples_Updated.txt',sep='\t',header=T)
gist_ids=gist_ids[!duplicated(gist_ids),]

rownames(gist_ids)=gist_ids$GISTIC_Label_ID

#Make same row order in Table with IDs as in the CNV-table:
gist_ids=gist_ids[colnames(cnvs)[4:1086],]
colnames(cnvs)[4:1086]=gist_ids$Case_ID
colnames(cnvs)[1]='Gene'
cnvs_1=cnvs[,c(1,4:ncol(cnvs))]
cnvs_1=cnvs_1[,colnames(cnvs_1) %in% c('Gene',samples_tab$Case_ID)]

cnvs_2=melt(cnvs_1,id=c('Gene'))
colnames(cnvs_2)[2:3]=c('Case','CNV_status')

write.table(cnvs_2,'GISTIC_results.thresholded.formatted.1064Cases.20230711.tsv',sep='\t',quote=F,row.names=F)
