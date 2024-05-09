library(dplyr)
library(reshape)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(ggrastr)
library(RColorBrewer)

all_vars=read.table('Vars_DETECTED_in_either_assay.ALL.10cancers.SAAV_SpliceForms.20230605.tsv',sep='\t',header=T)
all_vars_1=all_vars[,c('geneSymbol','Variant')]
all_vars_1=all_vars_1[!duplicated(all_vars_1),]


#Try also plottiing just for the variants of interest:
all_str=read.table('out/ASE_test_results.10cancers.20230702.tsv',sep='\t',header=T)

gene='CHD4'
Var='12:C6601981A'

r_tum=all_str

r_tum=merge(r_tum,all_vars_1)
tab=r_tum[r_tum$geneSymbol==gene,]
tab=tab[tab$Variant==Var,]

tab$DER_count=tab$Ref_count
tab$ANC_count=tab$Alt_count

tab$Status_2=tab$Status_1
tab$Status_2=ifelse(tab$Status_2=='Pref_REF', 'Pref_DER', tab$Status_2)
tab$Status_2=ifelse(tab$Status_2=='Pref_ALT', 'Pref_ANC', tab$Status_2)

#Add more extreme categories:
tab$Status_2=ifelse(tab$Status_2=='Pref_ANC' & tab$DER_count<=10, 'Pref_ANC_e', tab$Status_2)
tab$Status_2=ifelse(tab$Status_2=='Pref_DER' & tab$ANC_count<=10, 'Pref_DER_e', tab$Status_2)


p <- ggplot(tab,aes(x=(DER_count), y=(ANC_count), color = Status_2))

p <- p + geom_point_rast(alpha=0.8, stroke=0, size=2,aes(pch=Status))+geom_abline(linetype="longdash",color='grey') 

p <- p + theme_bw() + theme_classic()

p <- p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10))

p <- p + expand_limits(x = 0, y = 0) 

p <- p + scale_color_manual(values=c('None'='grey','Pref_DER'='#6e016b','Pref_DER_e'='#3f007d','Pref_ANC'='#b10026','Pref_ANC_e'='#67001f'))

p <- p + labs(x = "DER allele read count", y = "ANC allele read count")

pdf(paste("plotting.v.20230706/ASE_results.allTumor.10.cancers.v3.",gene,"_",Var,"20230706.pdf",sep=''), width=6.5, height=5,useDingbats=FALSE)
p
dev.off()

#For TP53 need to switch X and Y axises (Ref corresponds to ANC, and Alt corresponds to DER allele):
all_str=read.table('ASE_test_results.10cancers.20230702.tsv',sep='\t',header=T)

gene='TP53'
Var='17:G7676154C'
r_tum=all_str

r_tum=merge(r_tum,all_vars_1)
tab=r_tum[r_tum$geneSymbol==gene,]
tab=tab[tab$Variant==Var,]

tab$DER_count=tab$Alt_count
tab$ANC_count=tab$Ref_count

tab$Status_2=tab$Status_1
tab$Status_2=ifelse(tab$Status_2=='Pref_REF', 'Pref_ANC', tab$Status_2)
tab$Status_2=ifelse(tab$Status_2=='Pref_ALT', 'Pref_DER', tab$Status_2)

#Decided not to separate these, as no data in personalized proteomes for the samples with extreme/high preference:
#Add more extreme categories:
#tab$Status_2=ifelse(tab$Status_2=='Pref_ANC' & tab$DER_count<=10, 'Pref_ANC_e', tab$Status_2)
#tab$Status_2=ifelse(tab$Status_2=='Pref_DER' & tab$ANC_count<=10, 'Pref_DER_e', tab$Status_2)

p <- ggplot(tab,aes(x=(DER_count), y=(ANC_count), color = Status_2))

p <- p + geom_point_rast(alpha=0.8, stroke=0, size=2,aes(pch=Status))+geom_abline(linetype="longdash",color='grey') 

p <- p + theme_bw() + theme_classic()

p <- p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10))

p <- p + expand_limits(x = 0, y = 0) 

p <- p + scale_color_manual(values=c('None'='grey','Pref_DER'='#6e016b','Pref_ANC'='#b10026'))

p <- p + labs(x = "DER allele read count", y = "ANC allele read count")

pdf(paste("plotting.v.20230706/ASE_results.allTumor.10.cancers.v3.",gene,"_",Var,"20230706.pdf",sep=''), width=6.5, height=5,useDingbats=FALSE)
p
dev.off()




###############################
###Also do the same for NAT:###
###############################


gene='CHD4'
Var='12:C6601981A'

all_str=read.table('ASE_test_results.5cancers.NAT.20230702.tsv',sep='\t',header=T)

#Try also plottiing just for the variants of interest:
r_norm=all_str

r_norm=merge(r_norm,all_vars_1)
tab=r_norm[r_norm$geneSymbol==gene,]
tab=tab[tab$Variant==Var,]

tab$DER_count=tab$Ref_count
tab$ANC_count=tab$Alt_count

tab$Status_2=tab$Status_1
tab$Status_2=ifelse(tab$Status_2=='Pref_REF', 'Pref_DER', tab$Status_2)
tab$Status_2=ifelse(tab$Status_2=='Pref_ALT', 'Pref_ANC', tab$Status_2)

#Add more extreme categories:
tab$Status_2=ifelse(tab$Status_2=='Pref_ANC' & tab$DER_count<=10, 'Pref_ANC_e', tab$Status_2)
tab$Status_2=ifelse(tab$Status_2=='Pref_DER' & tab$ANC_count<=10, 'Pref_DER_e', tab$Status_2)

p <- ggplot(tab,aes(x=(DER_count), y=(ANC_count), color = Status_2))

p <- p + geom_point_rast(alpha=0.8, stroke=0, size=2,aes(pch=Status))+geom_abline(linetype="longdash",color='grey') 

p <- p + theme_bw() + theme_classic()

p <- p + theme(axis.title = element_text(size=16), axis.text.x = element_text(colour="black", size=10), axis.text.y = element_text(colour="black", size=10))

p <- p + expand_limits(x = 0, y = 0) 

p <- p + scale_color_manual(values=c('None'='grey','Pref_DER'='#6e016b','Pref_DER_e'='#3f007d','Pref_ANC'='#b10026','Pref_ANC_e'='#67001f'))

p <- p + labs(x = "DER allele read count", y = "ANC allele read count")

pdf(paste("plotting.v.20230706/NAT__ASE_results.allTumor.10.cancers.v3.",gene,"_",Var,"20230706.pdf",sep=''), width=6.5, height=5,useDingbats=FALSE)
p
dev.off()

##########################
#####ENDS HERE FOR NOW!###
##########################

