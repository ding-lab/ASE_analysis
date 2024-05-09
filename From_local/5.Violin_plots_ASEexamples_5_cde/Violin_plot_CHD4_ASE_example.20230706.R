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


#now do violin plots for examples:
ph_g=read.table('germlineOnly_phosphoSites_expr_quantiles.20230613.tsv',sep='\t',header=T)
ph=read.table('phosphosite_expr_quantiles.20230613.tsv',sep='\t',header=T)
ph$Type='Ref'
ph_g$Type='Alt'

tab=rbind(ph,ph_g)

tab=tab[tab$Gene %in% c('CHD4'),]


all_vars=read.table('Vars_DETECTED_in_either_assay.ALL.10cancers.SAAV_SpliceForms.20230605.tsv',sep='\t',header=T)


#Add GT information:
gts=read.table('DATA/GTs_for_SAAV_sites.20230702.tsv',sep='\t',header=T)
gts_1=gts[,c('Variant','GT','Sample')]


all_vars=all_vars[,c(1:5,7)]
all_vars=all_vars[!duplicated(all_vars),]
all_vars_1=all_vars[all_vars$geneSymbol=='CHD4',]

all_vars_2=merge(all_vars_1,gts_1)

#all_vars_1$Sample[!(all_vars_1$Sample %in% all_vars_2$Sample)]
#Sample C3N-02788 is missing in GTs (GBM sample). we don't have RNA-seq BAM-file for this one.

all_vars_2=all_vars_2[,c('accession_number', 'geneSymbol', 'Sample','GT','Variant')]
colnames(all_vars_2)[c(2,3)]=c('Gene','Case')

#Filter to the ones that are carriers of germline variants
tab_2=merge(tab,all_vars_2)
tab_2$accession_ph=ifelse(tab_2$Type=='Ref',gsub('(.*)\\..*','\\1',tab_2$ID),gsub('(.*)-.*','\\1',tab_2$ID))

#Limit to the ones with matching accessions between the accession from the mapping, and the accession from the phosphoSite data:
tab_3=tab_2[tab_2$accession_number==tab_2$accession_ph,]

#Identify phosphosites that were affected by Variant, and that should be plotted:
table(tab_3$ID[tab_3$Type=='Alt'])

#Filter to the phosphosite of interest (S145), the one which was affected by germline variant:
tab_4=tab_3[tab_3$ID %in% c('ENSP00000496634.1_S145s','ENSP00000496634-E139D_S145s'),]
colnames(tab_4)<-gsub('Type','VM_site_type',colnames(tab_4))

#Also keep only the variant of interest:
tab_4=tab_4[tab_4$Variant=='12:C6601981A',]

#Add sample type annotation:
tab_4$Type=tab_4$Sample
tab_4$Type=gsub('.*\\.T','Tumor',tab_4$Type)
tab_4$Type=gsub('.*\\.N','NAT',tab_4$Type)


#Rename Ref-->DER, and Alt-->ANC
tum$DER_count=tum$Ref_count
tum$ANC_count=tum$Alt_count
norm$DER_count=norm$Ref_count
norm$ANC_count=norm$Alt_count

#Add ASE results (from RNA-seq):
tum_s=tum[,c('Variant','Status_1','Sample','Cancer','VAF','ANC_count','DER_count')]
norm_s=norm[,c('Variant','Status_1','Sample','Cancer','VAF','ANC_count','DER_count')]
tum_s$Type='Tumor'
norm_s$Type='NAT'

both_ase=rbind(tum_s,norm_s)
colnames(both_ase)<-gsub('Sample','Case',colnames(both_ase))

#This will filter to those with GT 0/1, and also the ones with both ASE-results, and proteomic data.
tab_5=merge(tab_4, both_ase)

#Only the ones with GT=='0/1' are here:
to_plot=tab_5
to_plot$Status_2=to_plot$Status_1
to_plot$Status_2=ifelse(to_plot$Status_2=='None','0/1 None',to_plot$Status_2)
to_plot$Status_2=ifelse(to_plot$Status_2=='Pref_ALT','0/1 ASPE, ANC',to_plot$Status_2)
to_plot$Status_2=ifelse(to_plot$Status_2=='Pref_REF','0/1 ASPE, DER',to_plot$Status_2)

to_plot$Status_2=ifelse(to_plot$Status_2=='0/1 ASPE, ANC' & to_plot$DER_count<=10,'0/1 ASPE, ANC e.',to_plot$Status_2)
to_plot$Status_2=ifelse(to_plot$Status_2=='0/1 ASPE, DER' & to_plot$ANC_count<=10,'0/1 ASPE, DER e.',to_plot$Status_2)


to_plot$VM_site_type=ifelse(to_plot$VM_site_type=='Ref','DER','ANC')
to_plot$VM_site_type=factor(to_plot$VM_site_type,levels=c('DER','ANC'))


###Add RNA-seq VAF to the table:
#As we need only RNA-seq VAF, then make rows unique by Sample ID
to_plot_2=to_plot[!duplicated(to_plot$Sample),]

to_plot_vm=to_plot[,c('Status_2','Expr_quant_pancan','Type','Sample','VM_site_type')]
to_plot_rna=to_plot_2[,c('Status_2','VAF','Type','Sample')]
to_plot_rna$Data_type='RNA'

colnames(to_plot_vm) <- gsub('Expr_quant_pancan','value',colnames(to_plot_vm))
colnames(to_plot_vm) <- gsub('VM_site_type','Data_type',colnames(to_plot_vm))
colnames(to_plot_rna) <- gsub('VAF','value',colnames(to_plot_rna))

to_plot_final=rbind(to_plot_vm,to_plot_rna)

to_plot_final$Data_type=factor(to_plot_final$Data_type,levels=c('DER','ANC','RNA'))
to_plot_final$Status_2=factor(to_plot_final$Status_2,levels=c('0/1 ASPE, ANC e.','0/1 ASPE, ANC','0/1 None','0/1 ASPE, DER', '0/1 ASPE, DER e.'))

#Remove NA values:
to_plot_final=to_plot_final[!is.na(to_plot_final$value),]

p <- ggplot(data = to_plot_final, aes(x=Status_2, y=value)) 

p <- p + geom_violin(width=0.6,fill=NA, trim=T)

p <- p + stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width=0.7)

p <- p + geom_jitter(shape=19, position=position_jitter(0.2),size=0.2,alpha=0.5,aes(color=Status_2))

p <- p + facet_grid(Data_type~Type,drop=T,scales = "free", space = "free")

p <- p + theme_bw() + theme_nogrid() + theme_minimal() 

p <- p + labs(title="",x="",y="")

p <- p + theme(axis.text.x = element_text(colour="black", size=8,  angle=90, vjust = 1,hjust=1), axis.text.y = element_text(colour="black", size=10),axis.ticks = element_blank())+theme(legend.position = "none")

P <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))+theme(legend.position='bottom')

p <- p + ggtitle('CHD4-phosphosite (S145)')

p <- p + scale_color_manual(values=c('0/1 None'='#BEBEBE','0/1 ASPE, ANC'='#B11F28', '0/1 ASPE, ANC e.'='#67001f','0/1 ASPE, DER'='#6C2369', '0/1 ASPE DER e.'='#3f007d'))

pdf(paste("CDH4_varD139E.Expr_and_RNA_seq_VAF.20230706.pdf",sep=""), width=3.5, height=4,useDingbats=FALSE)
print(p)
dev.off()