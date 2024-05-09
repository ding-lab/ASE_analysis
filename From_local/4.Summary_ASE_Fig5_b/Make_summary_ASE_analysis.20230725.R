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


samples_tab=read.table('GO_case_list_CPTAC_germline.tsv',sep='\t',header=T)


cancers=c('BRCA','OV','CO','GBM','LUAD','CCRCC','UCEC', 'HNSCC','LSCC','PDA')

#reaad cancer genes
c_genes=read.table('DATA/gene_Lists/624_Cancer_genes.txt',sep='\t',header=F)


#First get all ASE results for Tumor and NAT
tum=read.table('ASE_test_results.10cancers.20230702.tsv',sep='\t',header=T)

norm=read.table('ASE_test_results.5cancers.NAT.20230702.tsv',sep='\t',header=T)

tum_s=tum[,c('Variant','Status_1','Sample','Cancer')]
norm_s=norm[,c('Variant','Status_1','Sample','Cancer')]
tum_s$Type='Tumor'
norm_s$Type='NAT'

both_r=rbind(tum_s,norm_s)
both_r$Count=1

#Get counts of ASE tests (ASE and None combined) for each SNP, for Tumor and NAT separately -- it corresponds to the N of the tested
st=aggregate(both_r$Count, by=list(both_r$Cancer,both_r$Type,both_r$Variant),FUN=sum)

#Filter to variants that have ASE for either  ALT or REF:
both_r1=both_r[both_r$Status_1 %in% c('Pref_ALT','Pref_REF'),]

#Here we calculate counts for ASE, summing UP both REF and ALT alleles counts -- it corresponds to the N with ASE
st_1=aggregate(both_r1$Count, by=list(both_r1$Cancer,both_r1$Type,both_r1$Variant),FUN=sum)

colnames(st)=c('Cancer','Type','Variant','Count_Total')
colnames(st_1)=c('Cancer','Type','Variant','Count')

res=merge(st_1,st)
res$Freq=res$Count/res$Count_Total

res$Count_Total=as.numeric(as.character(unlist(res$Count_Total)))


#Try doing separately for Tum & NAT:
res_1_t=res[res$Type=='Tumor',]
res_1_n=res[res$Type=='NAT',]


#Do the filtering of events to show (by frequency of ASE ecents, and total counts of ASEs); they are different for Tumor and NAT, as the majority are tumor ASEs, and we want to show some NAT ASEs.
res_1_t1=res_1_t$Variant[(res_1_t$Freq>0.5 & res_1_t$Count_Total>30)]
res_1_n1=res_1_n$Variant[(res_1_n$Freq>0.5 & res_1_n$Count_Total>30)]
res_1=res[res$Variant %in% c(res_1_t1, res_1_n1),]



#2023-07-25: Upd, probably would be good to lower this one to 0.2: 
res_1=res_1[res_1$Freq>0.2 & res_1$Count_Total>20,]


cols <- c("#386cb0","#386cb0","#386cb0","white","white","#d73027","#d73027","#d73027")
col_paletteB = colorRampPalette(brewer.pal(9,"Blues"))
col_paletteR = colorRampPalette(brewer.pal(9,"Reds"))
RdBu = brewer.pal(9, "RdBu") 
getPalette = colorRampPalette(rev(RdBu))
RdBu1024 = colorRampPalette(rev(RdBu))(1024)

########################
#Also add info about in what assay those variants were detected:
all_vars=read.table('Vars_DETECTED_in_either_assay.ALL.10cancers.SAAV_SpliceForms.20230605.tsv',sep='\t',header=T)
all_vars$SAAV=gsub('.*-(.*)','\\1',all_vars$ID)
all_vars_1=all_vars[,c('geneSymbol','Variant','Cancer','SAAV')]
all_vars_1=all_vars_1[!duplicated(all_vars_1),]

all_vars_2=all_vars_1
all_vars_2$ID=paste(all_vars_2$geneSymbol,all_vars_2$SAAV,sep='-')

#The N rows don't change, because those were already filtered to those detected in at least one assay.
res_2=merge(res_1,all_vars_2)

#Filter HLA-A/HLA-B variants (we will have 45 variants in total):
res_2=res_2[!(res_2$geneSymbol %in% c('HLA-A','HLA-B')),]


#Now save the file so that we can find the correct matching new IDs in the latest MAF-file.
write.table(res_2,'for_plot/Vars_selected_for_5B_initial_IDs.20230725.tsv',sep='\t',quote=F,row.names=F)




#Update ID using latest AA and DA alleles:
aa_da_annot=read.table('../inputs/ASE_annotated_withAA_DA.20230702.tsv',sep='\t',header=T)
aa_da_annot=aa_da_annot[,c('SNP_ID','Derived_HGVSp_Short','Changed_via_ancestral_determination')]
colnames(aa_da_annot)[1]='Variant'
aa_da_annot$Derived_HGVSp_Short=gsub('^p\\.','',aa_da_annot$Derived_HGVSp_Short)

#All variants were found.
res_3=merge(res_2,aa_da_annot)

#For consistency use the ones from the MAF file.
res_3$New_ID=paste(res_3$geneSymbol,res_3$Derived_HGVSp_Short,sep='-')

#Order ASE based on number of cancers that share the same ASE SNP in either Tumor or NAT:
st=res_3
st$Count=1
st_1=aggregate(st$Count, by=list(st$Variant),FUN='sum')
colnames(st_1)=c('Variant','ASE_count')

st_1$Status=ifelse(st_1$ASE_count>1,'Shared','Specific')
nat_ase=res_3$Variant[res_3$Type=='NAT']

#Add another level of ordering based on average frequency of ASE across cancers:
st_2=aggregate(st$Freq, by=list(st$Variant), FUN='mean')
colnames(st_2)=c('Variant','ASE_freq_mean')

st_both=merge(st_1,st_2)
st_both=st_both[order(-st_both$ASE_freq_mean),]
st_both=st_both[order(-st_both$ASE_count),]

res_3=res_3 %>% mutate(Variant=factor(Variant, levels=c(st_both$Variant))) %>% arrange(Variant) 
res_3$New_ID=factor(res_3$New_ID, levels=unique(res_3$New_ID))
res_3$New_ID_2=gsub('-',' ',res_3$New_ID)
res_3$New_ID_2=factor(res_3$New_ID_2, levels=unique(res_3$New_ID_2))

res_3$Status=ifelse(res_3$Variant %in% nat_ase, 'NAT_ASE','Tumor_ASE')

res_s=res_3[,c('ID','Variant','Status','New_ID_2')]


#Now also add information, if those ASE SNPs were identified in the proteomic data:
t1=all_vars[,c('geneSymbol','Variant','SAAV','Data')]
colnames(t1)[4]='Cancer'
t1=t1[!duplicated(t1),]

res_4=res_3[,c('Variant','geneSymbol','Cancer','Type','Count','Freq','SAAV','ID','Status','New_ID_2')]

###Set same params for all variants detected in the proteome.
t1$Count=20
t1$Freq=0.99
t1$Type='Tumor'

#Retain only variants that were detected in the proteome
t2=merge(t1,res_s) 

t2=t2[,colnames(res_4)]
res_7=rbind(res_4,t2)
res_7$Cancer=factor(res_7$Cancer,levels=c(cancers,c('Proteome','Phosphoproteome','Acetylome')))
res_7$Tech=ifelse(res_7$Cancer %in% cancers, 'WES','MS')

#2023-07-25: Final version of the plotting:

x=(colorRampPalette(append(brewer.pal(9,"YlOrRd"), "#FFFFFF", after=0))(1000))

      p <- ggplot()

      p <- p + geom_point(data = res_7, mapping = aes(x = New_ID_2, y = Type, fill = Freq, size = Count,shape=Tech), alpha =1,color='black')

      p <- p + scale_fill_gradientn(name= "Frequency", na.value=NA, colours=x)+scale_shape_manual(values=c('WES'=21,'MS'=22))

      p <- p + facet_grid(Cancer ~ Status, drop=T, space = "free",scales = "free")

      p <- p + xlab("") + ylab("")

      p <- p + guides(colour = guide_legend(title = "FDR"))

      p <- p + theme_bw()

      p <- p + theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"))

      p <- p + theme(axis.text.y = element_text(size = 10, face = "bold"))

      p <- p + theme(axis.text.x = element_text(size = 10, face = "bold", angle = 90, hjust = 0.95, vjust = 0.2))

      p <- p + theme(panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey80"))

      p <- p + theme(strip.text.x = element_text(angle = 0),strip.text.y= element_text(angle=0))
      
      p <- p + scale_size(range = c(1, 5),breaks=c(10, 30, 60, 90, 120))+theme(legend.position='right')

pdf(paste("plots/Summary_ASE.UpdCutoffs.1064Samples.v10.20230725.pdf",sep=""), width=12, height=4.5,useDingbats=FALSE)
print(p)
dev.off()

#Also save the data used for plottting:
write.table(res_7,'for_plot/For_summary_DotPlot_data.v9.Upd.Annotation.20230725.tsv',sep='\t',row.names=F,quote=F)

