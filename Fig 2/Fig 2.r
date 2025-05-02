### Fig. 2 ###

### load required packages ###
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)


### read in RNA-seq count matrix ###
sample_count <- read.table("./RNA_count_matrix.txt", header=TRUE)
sample_count <- column_to_rownames(sample_count, var="gene_name")

### read in patient metadata ###
patient_metadata <- read.csv("./patient_metadata.csv")
patient_metadata$Patient_ID <- paste0("Acc", patient_metadata$Patient_ID)
sample_metadata <- data.frame(sample_ID = colnames(sample_count))
sample_metadata <- separate(sample_metadata, col="sample_ID", into=c("Timepoint","Patient_ID"), sep="_")
sample_metadata <- left_join(sample_metadata, patient_metadata, by="Patient_ID")
rownames(sample_metadata) <- paste0(sample_metadata$Timepoint,"_", sample_metadata$Patient_ID)
sample_metadata$PPP2R1A_mutation <- factor(sample_metadata$PPP2R1A_mutation)
sample_metadata$ARID1A_mutation <- factor(sample_metadata$ARID1A_mutation)
sample_metadata$AKT_mutation <- factor(sample_metadata$AKT_mutation)

### Fig. 2a ###
sample_ID_pre_treatment <- rownames(sample_metadata)[which(sample_metadata$Timepoint == "Pre")]
sample_count_pre_treatment <- sample_count[,sample_ID_pre_treatment]
sample_metadata_pre_treatment <- sample_metadata[sample_ID_pre_treatment,]

dds_pre_treatment <- DESeqDataSetFromMatrix(countData = sample_count_pre_treatment,
                                            colData = sample_metadata_pre_treatment,
                                            design = ~ PPP2R1A_mutation)
dds_pre_treatment <- DESeq(dds_pre_treatment)

res_pre_treatment <- results(dds_pre_treatment)
res_pre_treatment_table <- as.data.frame(res_pre_treatment)
res_pre_treatment_table$Gene=rownames(res_pre_treatment_table)

res_pre_treatment_table <- filter(res_pre_treatment_table, !is.na(res_pre_treatment_table$stat))
Pre_PPP2R1A_Mut_vs_WT_genes <- res_pre_treatment_table %>% arrange(desc(stat)) %>% dplyr::select(Gene,stat)
Pre_PPP2R1A_Mut_vs_WT_genes <- tibble::deframe(Pre_PPP2R1A_Mut_vs_WT_genes)

hallmark <- msigdbr(species = "Homo sapiens", category = "H") 
hallmark <- hallmark[,c("gs_name", "gene_symbol")]

set.seed(13)
Pre_PPP2R1A_Mut_vs_WT_gsea <- GSEA(Pre_PPP2R1A_Mut_vs_WT_genes,TERM2GENE=hallmark, eps=0, pvalueCutoff = 1)
Pre_PPP2R1A_Mut_vs_WT_gsea <- as.data.frame(Pre_PPP2R1A_Mut_vs_WT_gsea)

Pre_PPP2R1A_Mut_vs_WT_gsea <- Pre_PPP2R1A_Mut_vs_WT_gsea%>%arrange(desc(NES), .by_group = F)
Pre_PPP2R1A_Mut_vs_WT_gsea$ID<-factor(Pre_PPP2R1A_Mut_vs_WT_gsea$ID,levels=rev(Pre_PPP2R1A_Mut_vs_WT_gsea$ID))
Pre_PPP2R1A_Mut_vs_WT_gsea$`-log10(FDR)`<- -log10(Pre_PPP2R1A_Mut_vs_WT_gsea$qvalue)
Pre_PPP2R1A_Mut_vs_WT_gsea$Color <- NA
Pre_PPP2R1A_Mut_vs_WT_gsea$Color[Pre_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.05 & Pre_PPP2R1A_Mut_vs_WT_gsea$NES >0] <- "FDR < 0.05, pos"
Pre_PPP2R1A_Mut_vs_WT_gsea$Color[Pre_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.01 & Pre_PPP2R1A_Mut_vs_WT_gsea$NES >0] <- "FDR < 0.01, pos"
Pre_PPP2R1A_Mut_vs_WT_gsea$Color[Pre_PPP2R1A_Mut_vs_WT_gsea$qvalue< 0.001 & Pre_PPP2R1A_Mut_vs_WT_gsea$NES >0] <- "FDR < 0.001, pos"
Pre_PPP2R1A_Mut_vs_WT_gsea$Color[Pre_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.0001 & Pre_PPP2R1A_Mut_vs_WT_gsea$NES >0] <- "FDR < 0.0001, pos"
Pre_PPP2R1A_Mut_vs_WT_gsea$Color[Pre_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.05 & Pre_PPP2R1A_Mut_vs_WT_gsea$NES < 0] <- "FDR < 0.05, neg"
Pre_PPP2R1A_Mut_vs_WT_gsea$Color[Pre_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.01 & Pre_PPP2R1A_Mut_vs_WT_gsea$NES < 0] <- "FDR < 0.01, neg"
Pre_PPP2R1A_Mut_vs_WT_gsea$Color[Pre_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.001 & Pre_PPP2R1A_Mut_vs_WT_gsea$NES < 0] <- "FDR < 0.001, neg"
Pre_PPP2R1A_Mut_vs_WT_gsea$Color[Pre_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.0001 & Pre_PPP2R1A_Mut_vs_WT_gsea$NES < 0] <- "FDR < 0.0001, neg"
Pre_PPP2R1A_Mut_vs_WT_gsea$Color <- factor(Pre_PPP2R1A_Mut_vs_WT_gsea$Color, levels = c("NS","FDR < 0.05, pos","FDR < 0.01, pos","FDR < 0.001, pos","FDR < 0.0001, pos",
                                                                                        "FDR < 0.05, neg","FDR < 0.01, neg","FDR < 0.001, neg","FDR < 0.0001, neg"))


ggplot(Pre_PPP2R1A_Mut_vs_WT_gsea,aes(x=NES, y=`-log10(FDR)`, color=Color)) +
  geom_point(aes(size=abs(NES))) +
  scale_colour_manual(values = c(
    `FDR < 0.0001, neg` = "#081d58",
    `FDR < 0.001, neg` = "#225ea8",
    `FDR < 0.01, neg` = "#1d91c0",
    `FDR < 0.05, neg` = "#deebf7",
    `FDR < 0.0001, pos` = "#662506",
    `FDR < 0.001, pos` = "#993404",
    `FDR < 0.01, pos` = "#cc4c02",
    `FDR < 0.05, pos` = "#ec7014",
    `NS` = "#878787"),
    guide = guide_legend(override.aes = list(size = 4))) +
  ggtitle("Pre-treatment PPP2R1A Mut vs. WT GSEA") +
  #ylim(0,150) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "right") +
  labs(x = "Enriched in WT <- NES -> Enriched in Mut", y = "-log10(FDR)", color = "Significance", size="|NES|") +
  geom_text_repel(data = Pre_PPP2R1A_Mut_vs_WT_gsea[Pre_PPP2R1A_Mut_vs_WT_gsea$`-log10(FDR)`>-log10(0.05),],
                  aes(label = gsub('HALLMARK_','', rownames(Pre_PPP2R1A_Mut_vs_WT_gsea[Pre_PPP2R1A_Mut_vs_WT_gsea$`-log10(FDR)`>-log10(0.05),]))),size = 3,
                  box.padding = unit(0.35, "lines"),  point.padding = unit(0.3, "lines"), color="black", max.overlaps = 100) 


### Fig. 2b ###
sample_ID_on_treatment <- rownames(sample_metadata)[which(sample_metadata$Timepoint == "On")]
sample_count_on_treatment <- sample_count[,sample_ID_on_treatment]
sample_metadata_on_treatment <- sample_metadata[sample_ID_on_treatment,]

dds_on_treatment <- DESeqDataSetFromMatrix(countData = sample_count_on_treatment,
                                           colData = sample_metadata_on_treatment,
                                           design = ~ PPP2R1A_mutation)
dds_on_treatment <- DESeq(dds_on_treatment)

res_on_treatment <- results(dds_on_treatment)
res_on_treatment_table <- as.data.frame(res_on_treatment)
res_on_treatment_table$Gene=rownames(res_on_treatment_table)

res_on_treatment_table <- filter(res_on_treatment_table, !is.na(res_on_treatment_table$stat))
On_PPP2R1A_Mut_vs_WT_genes <- res_on_treatment_table %>% arrange(desc(stat)) %>% dplyr::select(Gene,stat)
On_PPP2R1A_Mut_vs_WT_genes <- tibble::deframe(On_PPP2R1A_Mut_vs_WT_genes)

set.seed(13)
On_PPP2R1A_Mut_vs_WT_gsea <- GSEA(On_PPP2R1A_Mut_vs_WT_genes,TERM2GENE=hallmark, eps=0, pvalueCutoff = 1)
On_PPP2R1A_Mut_vs_WT_gsea <- as.data.frame(On_PPP2R1A_Mut_vs_WT_gsea)

On_PPP2R1A_Mut_vs_WT_gsea <- On_PPP2R1A_Mut_vs_WT_gsea%>%arrange(desc(NES), .by_group = F)
On_PPP2R1A_Mut_vs_WT_gsea$ID<-factor(On_PPP2R1A_Mut_vs_WT_gsea$ID,levels=rev(On_PPP2R1A_Mut_vs_WT_gsea$ID))
On_PPP2R1A_Mut_vs_WT_gsea$`-log10(FDR)`<- -log10(On_PPP2R1A_Mut_vs_WT_gsea$qvalue)
On_PPP2R1A_Mut_vs_WT_gsea$Color <- NA
On_PPP2R1A_Mut_vs_WT_gsea$Color[On_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.05 & On_PPP2R1A_Mut_vs_WT_gsea$NES >0] <- "FDR < 0.05, pos"
On_PPP2R1A_Mut_vs_WT_gsea$Color[On_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.01 & On_PPP2R1A_Mut_vs_WT_gsea$NES >0] <- "FDR < 0.01, pos"
On_PPP2R1A_Mut_vs_WT_gsea$Color[On_PPP2R1A_Mut_vs_WT_gsea$qvalue< 0.001 & On_PPP2R1A_Mut_vs_WT_gsea$NES >0] <- "FDR < 0.001, pos"
On_PPP2R1A_Mut_vs_WT_gsea$Color[On_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.0001 & On_PPP2R1A_Mut_vs_WT_gsea$NES >0] <- "FDR < 0.0001, pos"
On_PPP2R1A_Mut_vs_WT_gsea$Color[On_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.05 & On_PPP2R1A_Mut_vs_WT_gsea$NES < 0] <- "FDR < 0.05, neg"
On_PPP2R1A_Mut_vs_WT_gsea$Color[On_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.01 & On_PPP2R1A_Mut_vs_WT_gsea$NES < 0] <- "FDR < 0.01, neg"
On_PPP2R1A_Mut_vs_WT_gsea$Color[On_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.001 & On_PPP2R1A_Mut_vs_WT_gsea$NES < 0] <- "FDR < 0.001, neg"
On_PPP2R1A_Mut_vs_WT_gsea$Color[On_PPP2R1A_Mut_vs_WT_gsea$qvalue < 0.0001 & On_PPP2R1A_Mut_vs_WT_gsea$NES < 0] <- "FDR < 0.0001, neg"
On_PPP2R1A_Mut_vs_WT_gsea$Color <- factor(On_PPP2R1A_Mut_vs_WT_gsea$Color, levels = c("NS","FDR < 0.05, pos","FDR < 0.01, pos","FDR < 0.001, pos","FDR < 0.0001, pos",
                                                                                        "FDR < 0.05, neg","FDR < 0.01, neg","FDR < 0.001, neg","FDR < 0.0001, neg"))

ggplot(On_PPP2R1A_Mut_vs_WT_gsea,aes(x=NES, y=`-log10(FDR)`, color=Color)) +
  geom_point(aes(size=abs(NES))) +
  scale_colour_manual(values = c(
    `FDR < 0.0001, neg` = "#081d58",
    `FDR < 0.001, neg` = "#225ea8",
    `FDR < 0.01, neg` = "#1d91c0",
    `FDR < 0.05, neg` = "#deebf7",
    `FDR < 0.0001, pos` = "#662506",
    `FDR < 0.001, pos` = "#993404",
    `FDR < 0.01, pos` = "#cc4c02",
    `FDR < 0.05, pos` = "#ec7014",
    `NS` = "#878787"),
    guide = guide_legend(override.aes = list(size = 4))) +
  ggtitle("On-treatment PPP2R1A Mut vs. WT GSEA") +
  #ylim(0,150) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "right") +
  labs(x = "Enriched in WT <- NES -> Enriched in Mut", y = "-log10(FDR)", color = "Significance", size="|NES|") +
  #geom_text_repel(data = head(subset(result, gn.selected),30),aes(label = rownames(head(subset(result, gn.selected),30))),size = 2,
  geom_text_repel(data = On_PPP2R1A_Mut_vs_WT_gsea[On_PPP2R1A_Mut_vs_WT_gsea$`-log10(FDR)`>-log10(0.05),],
                  aes(label = gsub('HALLMARK_','', rownames(On_PPP2R1A_Mut_vs_WT_gsea[On_PPP2R1A_Mut_vs_WT_gsea$`-log10(FDR)`>-log10(0.05),]))),size = 3,
                  box.padding = unit(0.35, "lines"),  point.padding = unit(0.3, "lines"), color="black", max.overlaps = 100) 


### Fig. 2c-f ###
### Immune cell deconvolution is performed with CIBERSORTx ###
### read in deconvolution results ###
Cibersort_res <- read.csv("./Cibersort_res.csv")
sample_metadata$Sample_ID <- rownames(sample_metadata)
Cibersort_res <- left_join(Cibersort_res, sample_metadata, by="Sample_ID")

paired_patient_ID <- c('Acc160','Acc162','Acc166','Acc119','Acc121','Acc124','Acc137','Acc139','Acc149','Acc150','Acc19','Acc31','Acc36','Acc38')
Cibersort_res_paired <- filter(Cibersort_res, Cibersort_res$Patient_ID %in% paired_patient_ID)
Cibersort_res_paired$Timepoint_PPP2R1A <- paste0(Cibersort_res_paired$Timepoint, "_", Cibersort_res_paired$PPP2R1A_mutation)
Cibersort_res_paired$Timepoint_PPP2R1A <- factor(Cibersort_res_paired$Timepoint_PPP2R1A, levels=c("Pre_1","On_1","Pre_0","On_0"))

ggplot(Cibersort_res_paired, aes(x=Timepoint_PPP2R1A, y=T.cells.CD8))+
    geom_boxplot(outlier.size = 0,width=0.5, aes(color=Timepoint_PPP2R1A))+
    geom_point(aes(color=Timepoint_PPP2R1A), size=3.5, alpha=1)+
    geom_line(aes(group=Patient_ID),color="#969696")+
    scale_color_manual(values=c("#fee391","#ec7014","#fee391","#ec7014"))+
    theme_classic()+
    theme(axis.ticks.length = unit(0.2,'cm'))+
    theme(legend.position = "none", plot.title = element_text(size=14), 
          axis.title = element_text(size=14), 
          axis.text.x = element_text(angle = 45, hjust=1, size=10))+
    guides(fill=guide_legend(title = NULL))+
    xlab("Timepoint_PPP2R1A")+
    ylab("Prop. of cells")+
    ggtitle("CD8 T cells")

### The lognitudinal changes of other cells are plotted in a similar way.


### Fig. 2g-h ###
### TCR and BCR reconstruction was performed using TRUST4 ###
### read in TCR/BCR results ###
TCR_BCR_richness <- read.csv("./TCR_BCR_richness.csv")
TCR_BCR_richness$Sample_ID <- paste0(TCR_BCR_richness$Timepoint, "_", TCR_BCR_richness$Patient_ID)
TCR_BCR_richness$Patient_ID <- NULL; TCR_BCR_richness$Timepoint <- NULL
TCR_BCR_richness$TCR.Volume <- log10(TCR_BCR_richness$TCR.Volume)
TCR_BCR_richness$BCR.Volume <- log10(TCR_BCR_richness$BCR.Volume)
TCR_BCR_richness <- left_join(TCR_BCR_richness, sample_metadata, by="Sample_ID")

TCR_paired_patient_ID <- c('Acc160','Acc162','Acc166','Acc119','Acc121','Acc124','Acc137','Acc139','Acc149','Acc150','Acc19','Acc31','Acc36')
TCR_richness_paired <- filter(TCR_BCR_richness, TCR_BCR_richness$Patient_ID %in% TCR_paired_patient_ID)
TCR_richness_paired$Timepoint_PPP2R1A <- paste0(TCR_richness_paired$Timepoint, "_", TCR_richness_paired$PPP2R1A_mutation)
TCR_richness_paired$Timepoint_PPP2R1A <- factor(TCR_richness_paired$Timepoint_PPP2R1A, levels=c("Pre_1","On_1","Pre_0","On_0"))

ggplot(TCR_richness_paired, aes(x=Timepoint_PPP2R1A, y=TCR.Volume))+
  geom_boxplot(outlier.size = 0,width=0.5, aes(color=Timepoint_PPP2R1A))+
  geom_point(aes(color=Timepoint_PPP2R1A), size=3, alpha=1)+
  geom_line(aes(group=Patient_ID),color="#969696")+
  scale_color_manual(values=c("#fee391","#ec7014","#fee391","#ec7014"))+
  theme_classic()+
  theme(axis.ticks.length = unit(0.2,'cm'))+
  theme(legend.position = "none", plot.title = element_text(size=14), 
        axis.title = element_text(size=14), 
        axis.text.x = element_text(angle = 45, hjust=1, size=10))+
  guides(fill=guide_legend(title = NULL))+
  xlab("Timepoint_PPP2R1A")+
  ylab("TCR.Volume")+
  ggtitle("TCR.Volume")


BCR_paired_patient_ID <- c('Acc160','Acc162','Acc166','Acc119','Acc121','Acc124','Acc137','Acc139','Acc149','Acc150','Acc19','Acc31','Acc36','Acc38')
BCR_richness_paired <- filter(TCR_BCR_richness, TCR_BCR_richness$Patient_ID %in% BCR_paired_patient_ID)
BCR_richness_paired$Timepoint_PPP2R1A <- paste0(BCR_richness_paired$Timepoint, "_", BCR_richness_paired$PPP2R1A_mutation)
BCR_richness_paired$Timepoint_PPP2R1A <- factor(BCR_richness_paired$Timepoint_PPP2R1A, levels=c("Pre_1","On_1","Pre_0","On_0"))

ggplot(BCR_richness_paired, aes(x=Timepoint_PPP2R1A, y=BCR.Volume))+
  geom_boxplot(outlier.size = 0,width=0.5, aes(color=Timepoint_PPP2R1A))+
  geom_point(aes(color=Timepoint_PPP2R1A), size=3, alpha=1)+
  geom_line(aes(group=Patient_ID),color="#969696")+
  scale_color_manual(values=c("#fee391","#ec7014","#fee391","#ec7014"))+
  theme_classic()+
  theme(axis.ticks.length = unit(0.2,'cm'))+
  theme(legend.position = "none", plot.title = element_text(size=14), 
        axis.title = element_text(size=14), 
        axis.text.x = element_text(angle = 45, hjust=1, size=10))+
  guides(fill=guide_legend(title = NULL))+
  xlab("Timepoint_PPP2R1A")+
  ylab("BCR.Volume")+
  ggtitle("BCR.Volume")
