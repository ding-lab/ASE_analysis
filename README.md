# Allele specific expression (ASE) analysis

This is a repo with the scripts for ASE analysis.

  * Allele specific expression:

   + (https://stat.ethz.ch/R-manual/R-devel/library/stats/html/binom.test.html).

   + https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00964-1.

# Bam-readcount installed:
```/diskmnt/Software/bam-readcount-0.7.4/mybuild/bin/bam-readcount```

# Reference:
```/diskmnt/Datasets/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.fa```


# Analysis scripts:

* The order follows the steps that were done for analysis in the paper.

1. ```From_local/1.Make_site_list/Make_site_list.20230702.R``` -- make list of sites to run bam-readcount on cluster.


2. ```From_cluster/1.Extract_GTs``` -- folder with a script to extract genotypes (GTs) for the variants of interest. We would need only heterozygous variants.


3. ```From_cluster/2.BamReadCount``` -- folder with the scripts to run bamreadcount analysis on the RNA-seq data for the variants of interest. We will use the lists of sites for all samples to run bamreadcount on the respective RNA-seq BAM-files. 


2. ```From_local/2.ASE_test/ASE_RNA_seq_1064Samples.20230702.R``` -- use results of bam-readcount to run binomial test on the read-counts of REF and ALT alleles. It uses also the latest genotypes (GTs) from the VCF file.


3. ```From_local/3.ASE_example_Figs5_cde/ASE_RNA_seq.plottingExamples.20230706.R``` -- script to make plot for read counts in REF anf ALT for the two selected variants (colored by the ASE status).


4. ```From_local/4.Summary_ASE_Fig5_b/Make_summary_ASE_analysis.20230725.R``` -- script to make summary dotplot of the most frequent ASE events in tumor and NAT.


5. ```From_local/5.Violin_plots_ASEexamples_5_cde ``` -- this folder contains scripts to make violin plots for the peptide abundances for the two selected variants in CHD4 and TP53.


6. ```From_local/6.ASE_CNV_annotated/Barplot_ASE_annotated_withCNVs.20230711.R``` -- script for annotation of ASE events with CNV status for the respective gene.
