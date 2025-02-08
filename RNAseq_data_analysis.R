library(dplyr)
library(tidyverse)
library(limma)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(fgsea)
library(enrichplot)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(pheatmap)


########################################## data processing ###################################################################
{
file<-read.table('sample_HTSeq.txt',header=F)

result<-NULL;

for(i in 1:nrow(file)){
  data<-read.table(file[i,2],header = F,sep='\t');
  rownames(data)<-data$V1;
  result<-cbind(result,data$V2);
  colnames(result)[i]<-gsub('.*-.*-(.*-.*).txt','\\1',file[i,1])
}

rownames(result)<-rownames(data);
result<-result[-c((nrow(data)-4):nrow(data)),];

if(!require(rtracklayer)) stop('error loading rtracklayer')
gff = readGFF('gencode.v19.nochr.annotation.gtf', tags = c('gene_id','gene_type','gene_name','transcript_id'));
gff<-gff[gff$type=='gene',]

resultall<-data.frame(gff[match(rownames(result),gff$gene_id),],result,check.names = F);
write.table(resultall,'HTSeq_combined_read_count.txt',row.names = F,quote=F,sep='\t')

####coding
resultcoding<-resultall[resultall$gene_type=='protein_coding',];
resultcoding<-resultcoding[,-c(1:10,12)];
colnames(resultcoding)[1]<-'Gene.symbol'
write.table(resultcoding,'HTSeq_combined_read_count_protein-coding.txt',row.names = F,quote=F,sep='\t')

####linc.miRNA
result.linc.miRNA<-resultall[resultall$gene_type=='lincRNA' |resultall$gene_type=='miRNA',];
result.linc.miRNA<-result.linc.miRNA[,-c(1:10,12)];
colnames(result.linc.miRNA)[1]<-'Gene.symbol'
write.table(result.linc.miRNA,'HTSeq_combined_read_count_linc-miRNA.txt',row.names = F,quote=F,sep='\t')

###FPKM, TPM
geneLen<-read.table('gencode.v19.nochr.annotation_exon_geneLen.txt',header = T,sep='\t');

resultcoding<-na.omit(resultcoding);
countMtr = select(resultcoding, -Gene.symbol)
libSz = apply(countMtr, 2, sum) ##RCpc Number of reads mapped to all protein-coding genes
geneLen<-geneLen[match(resultcoding$Gene.symbol,geneLen$gene_name),]
geneLen = unlist(geneLen$GeneLen) ##L Length of the gene in base pairs; Calculated as the sum of all exons in a gene

fpkm = t(apply(countMtr*1e9,1,function(rw,libSz) rw/libSz,libSz = libSz))
fpkm = apply(fpkm, 2, function(cl,geneLen) cl/geneLen,geneLen = geneLen)

fpkm = tbl_df(data.frame(
  Gene.symbol=resultcoding$Gene.symbol,
  round(fpkm,4),
  check.names =F
))

write.table(fpkm,'HTSeq_combined_read_count_protein-coding_fpkm.txt',row.names = F,quote=F,sep='\t');


aa=apply(countMtr,2,function(x){x*10^3/geneLen})
tpm=apply(aa,2,function(x){x*10^6/sum(x)})

#colnames(fpkm)[1]<-'Gene.symbol';
tpm=tbl_df(data.frame(
  Gene.symbol=resultcoding$Gene.symbol,
  round(tpm,4),
  check.names =F
))

write.table(tpm,'HTSeq_combined_read_count_protein-coding_tpm.txt',row.names = F,quote=F,sep='\t')

#tpm=aggregate(tpm[,-1],by=list(Gene.symbol=tpm$Gene.symbol),FUN=sum)
#tpm<-read.table('/Volumes/scratch/AJ-Immune_RNA730/process/data/HTSeq_combined_read_count_protein-coding_tpm.txt',header=T)
log2tpm=data.frame('Gene.symbol'=tpm$Gene.symbol,round(log2(tpm[,-1]+1),4),check.names = F)
write.table(log2tpm,'HTSeq_combined_read_count_protein-coding_tpm_log2.txt',row.names = F,quote=F,sep='\t')
}


########################################## DEG pre-treatment  PPP2R1Amut vs. wt #########################################################

#sample_table_pre_treatment_filter <- filter(sample_table_pre_treatment, sample_table_pre_treatment$AKT.alt=="No")
#resultcoding_count_pre_treatment_filter <- resultcoding_count_pre_treatment[, rownames(sample_table_pre_treatment_filter)]
dds_PPP2R1A <- DESeqDataSetFromMatrix(countData = resultcoding_count_pre_treatment,
                                      colData = sample_table_pre_treatment,
                                      design = ~ PPP2R1A.mut)
dds_PPP2R1A <- DESeq(dds_PPP2R1A)
res_PPP2R1A <- results(dds_PPP2R1A)
resultsNames(dds_PPP2R1A)

res2_PPP2R1A_shrink <- lfcShrink(dds_PPP2R1A, coef="PPP2R1A.mut_Yes_vs_No", type="ashr")
res2_PPP2R1A_shrink_table <- as.data.frame(res2_PPP2R1A_shrink)
#res2_PPP2R1A
#summary(res2_PPP2R1A)

res2_PPP2R1A_table <- as.data.frame(res_PPP2R1A)
res2_PPP2R1A_table$invert_P <- (-log10(res2_PPP2R1A_table$padj)) 
res2_PPP2R1A_table$Gene=rownames(res2_PPP2R1A_table)
res2_PPP2R1A_table$Color <- "NS"
res2_PPP2R1A_table$Color[res2_PPP2R1A_table$padj < 0.1] <- "FDR < 0.1"
res2_PPP2R1A_table$Color[res2_PPP2R1A_table$padj < 0.01] <- "FDR < 0.01"
res2_PPP2R1A_table$Color[res2_PPP2R1A_table$padj< 0.001] <- "FDR < 0.001"
res2_PPP2R1A_table$Color[res2_PPP2R1A_table$padj < 0.0001] <- "FDR < 0.0001"
res2_PPP2R1A_table$Color <- factor(res2_PPP2R1A_table$Color, levels = c("NS","FDR < 0.1","FDR < 0.01","FDR < 0.001","FDR < 0.0001"))

ggplot(res2_PPP2R1A_table, aes(x = log2FoldChange, y = invert_P, color = Color,label = Gene)) +
  geom_vline(xintercept = c(2, -2), lty = "dashed") +
  geom_hline(yintercept = -log10(0.1), lty = "dashed") +
  geom_point(size=0.3) +
  labs(x = "Enriched in WT <- log2(FC) -> Enriched in PPP2R1A Mut",
       y = "-log10(p.adj)",
       color = "Significance") +
  scale_color_manual(values = c(
    `FDR < 0.0001` = "#08306b",
    `FDR < 0.001` = "#2171b5",
    `FDR < 0.01` = "#6baed6",
    `FDR < 0.1` = "#c6dbef",
    `NS` = "#f0f0f0"),
    guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(res2_PPP2R1A_table, padj < 0.01 & abs(log2FoldChange) > 2),
                  size = 3, point.padding = 0.2, color = "black",
                  min.segment.length = .5, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  ggtitle("DEGs in pre-treatment samples")+
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom") 


res2_PPP2R1A_table_top30 <- res2_PPP2R1A_table %>% filter(res2_PPP2R1A_table$padj<0.1 & res2_PPP2R1A_table$log2FoldChange>0) %>% slice_max(n=30, order_by = stat) 
res2_PPP2R1A_table_bottom30 <- res2_PPP2R1A_table %>% filter(res2_PPP2R1A_table$padj<0.1 & res2_PPP2R1A_table$log2FoldChange<0) %>% slice_min(n=30, order_by = stat) 
res2_PPP2R1A_top_bottom_DEGs <- c(rownames(res2_PPP2R1A_table_top30),rownames(res2_PPP2R1A_table_bottom30))



data_PPP2R1A<-log2tpm[immune_related_genes$symbol,
                      c(rownames(sample_table_pre_treatment)[sample_table_pre_treatment$PPP2R1A.mut=="No"],
                        rownames(sample_table_pre_treatment)[sample_table_pre_treatment$PPP2R1A.mut=="Yes"])]

data_PPP2R1A<-t(scale(t(data_PPP2R1A)))
data_PPP2R1A[data_PPP2R1A>1.5]=1.5
data_PPP2R1A[data_PPP2R1A< -1.5]= -1.5

sample_table_pre_treatment_for_anno <- sample_table_pre_treatment
colnames(sample_table_pre_treatment_for_anno)[11]<-"PPP2R1A mut"
sample_table_pre_treatment_for_anno$`PPP2R1A mut` <- as.numeric(sample_table_pre_treatment_for_anno$`PPP2R1A mut`)
sample_table_pre_treatment_for_anno$`PPP2R1A mut`[sample_table_pre_treatment_for_anno$`PPP2R1A mut`== 0]<-"No"
sample_table_pre_treatment_for_anno$`PPP2R1A mut`[sample_table_pre_treatment_for_anno$`PPP2R1A mut`== 1]<-"Yes"

pheatmap::pheatmap(data_PPP2R1A,show_colnames =T,show_rownames =T, fontsize_row = 7,cluster_rows = F, cluster_cols = F,
                   color=colorRampPalette(c('blue','black','yellow'))(100),
                   annotation_col=sample_table_pre_treatment_for_anno[,c('PPP2R1A.mut'),drop=FALSE])

########################################## *Fig 3a GSEA pre-treatment  PPP2R1Amut vs. wt #########################################################

res2_PPP2R1A_table <- filter(res2_PPP2R1A_table, !is.na(res2_PPP2R1A_table$stat))
Pre_PPP2R1A_Mut_vs_WT_genes <- res2_PPP2R1A_table %>% arrange(desc(stat)) %>% dplyr::select(Gene,stat)
Pre_PPP2R1A_Mut_vs_WT_genes <- tibble::deframe(Pre_PPP2R1A_Mut_vs_WT_genes)

TLS_signatue <- data.frame(gs_name="TLS", gene_symbol=c("CD19","CD1D","CD22","CD3D","CD3E","CD52","CD74","CD79A",
                                                        "CD79B","CD8A","CD8B","CR2","CXCL13","CXCR5","FCER2","IL7R",
                                                        "ITGAE","MS4A1","PDCD1","PTGDS","CD4","TRBC2"))
hallmark <- rbind(hallmark, TLS_signatue)
hallmark <- hallmark[-c(8210:8231),]

set.seed(13)
Pre_PPP2R1A_Mut_vs_WT_gsea <- GSEA(Pre_PPP2R1A_Mut_vs_WT_genes,TERM2GENE=hallmark, eps=0, pvalueCutoff = 1)
Pre_PPP2R1A_Mut_vs_WT_gsea <- as.data.frame(Pre_PPP2R1A_Mut_vs_WT_gsea)

Pre_PPP2R1A_Mut_vs_WT_gsea <- Pre_PPP2R1A_Mut_vs_WT_gsea%>%arrange(desc(NES), .by_group = F)
Pre_PPP2R1A_Mut_vs_WT_gsea$ID<-factor(Pre_PPP2R1A_Mut_vs_WT_gsea$ID,levels=rev(Pre_PPP2R1A_Mut_vs_WT_gsea$ID))
Pre_PPP2R1A_Mut_vs_WT_gsea$`-log10(FDR)`<- -log10(Pre_PPP2R1A_Mut_vs_WT_gsea$qvalue)

p<-ggplot(Pre_PPP2R1A_Mut_vs_WT_gsea, aes(ID, NES)) +
  geom_col(aes(fill=`-log10(FDR)`)) +
  #scale_fill_viridis() +
  #scale_fill_gradientn(colors=jet()[c(1:11,15:25)],na.value = "#e6e6e6") +
  coord_flip() +
  labs(x="", y="NES", title="")+
  theme(axis.text.x = element_text(size = 12,angle=0,vjust=0.2,hjust=0.95),
        strip.text.x = element_text(size=12, color="black", face="bold"),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

ggsave('On_vs_Pre.rm.outliers.NES.vs.FDR.barplot.pdf',p,height=10,width=10)

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
  geom_hline(yintercept=-log10(0.05), lty = "dashed")+
  labs(x = "Enriched in WT <- NES -> Enriched in Mut", y = "-log10(FDR)", color = "Significance", size="|NES|") +
  #geom_text_repel(data = head(subset(result, gn.selected),30),aes(label = rownames(head(subset(result, gn.selected),30))),size = 2,
  geom_text_repel(data = Pre_PPP2R1A_Mut_vs_WT_gsea[Pre_PPP2R1A_Mut_vs_WT_gsea$`-log10(FDR)`>-log10(0.05),],
                  aes(label = gsub('HALLMARK_','', rownames(Pre_PPP2R1A_Mut_vs_WT_gsea[Pre_PPP2R1A_Mut_vs_WT_gsea$`-log10(FDR)`>-log10(0.05),]))),size = 3,
                  box.padding = unit(0.35, "lines"),  point.padding = unit(0.3, "lines"), color="black", max.overlaps = 100) 


ggsave(paste0('On_vs_Pre.rm.outliers.NES.vs.FDR.',Sys.Date(),'.pdf'),plotx,width = 6,height = 6)

write.csv(Pre_PPP2R1A_Mut_vs_WT_gsea,"Fig3a_Pre_PPP2R1A_Mut_vs_WT_gsea.csv")








########################################## DEG pre-treatment  PPP2R1Amut AKTwt vs. PPP2R1Awt AKTwt #########################################################

sample_table_pre_treatment_filter <- filter(sample_table_pre_treatment, sample_table_pre_treatment$AKT.alt=="No")
resultcoding_count_pre_treatment_filter <- resultcoding_count_pre_treatment[, rownames(sample_table_pre_treatment_filter)]
dds_PPP2R1A_no_alt <- DESeqDataSetFromMatrix(countData = resultcoding_count_pre_treatment_filter,
                                              colData = sample_table_pre_treatment_filter,
                                              design = ~ PPP2R1A.mut)
dds_PPP2R1A_no_alt <- DESeq(dds_PPP2R1A_no_alt)
res_PPP2R1A_no_alt <- results(dds_PPP2R1A_no_alt)

#res2_PPP2R1A_shrink <- lfcShrink(dds_PPP2R1A, coef="PPP2R1A.mut_Yes_vs_No", type="ashr")
#res2_PPP2R1A_shrink_table <- as.data.frame(res2_PPP2R1A_shrink)
#res2_PPP2R1A
#summary(res2_PPP2R1A)

res2_PPP2R1A_no_alt_table <- as.data.frame(res_PPP2R1A_no_alt)
res2_PPP2R1A_no_alt_table$invert_P <- (-log10(res2_PPP2R1A_no_alt_table$padj)) 
res2_PPP2R1A_no_alt_table$Gene=rownames(res2_PPP2R1A_no_alt_table)
res2_PPP2R1A_no_alt_table$Color <- "NS"
res2_PPP2R1A_no_alt_table$Color[res2_PPP2R1A_no_alt_table$padj < 0.1] <- "FDR < 0.1"
res2_PPP2R1A_no_alt_table$Color[res2_PPP2R1A_no_alt_table$padj < 0.01] <- "FDR < 0.01"
res2_PPP2R1A_no_alt_table$Color[res2_PPP2R1A_no_alt_table$padj< 0.001] <- "FDR < 0.001"
res2_PPP2R1A_no_alt_table$Color[res2_PPP2R1A_no_alt_table$padj < 0.0001] <- "FDR < 0.0001"
res2_PPP2R1A_no_alt_table$Color <- factor(res2_PPP2R1A_no_alt_table$Color, levels = c("NS","FDR < 0.1","FDR < 0.01","FDR < 0.001","FDR < 0.0001"))

ggplot(res2_PPP2R1A_no_alt_table, aes(x = log2FoldChange, y = invert_P, color = Color,label = Gene)) +
  geom_vline(xintercept = c(2, -2), lty = "dashed") +
  geom_hline(yintercept = -log10(0.1), lty = "dashed") +
  geom_point(size=0.3) +
  labs(x = "Enriched in WT <- log2(FC) -> Enriched in PPP2R1A Mut",
       y = "-log10(p.adj)",
       color = "Significance") +
  scale_color_manual(values = c(
    `FDR < 0.0001` = "#08306b",
    `FDR < 0.001` = "#2171b5",
    `FDR < 0.01` = "#6baed6",
    `FDR < 0.1` = "#c6dbef",
    `NS` = "#f0f0f0"),
    guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(res2_PPP2R1A_no_alt_table, padj < 0.01 & abs(log2FoldChange) > 2),
                  size = 3, point.padding = 0.2, color = "black",
                  min.segment.length = .5, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  ggtitle("DEGs in pre-treatment samples")+
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom") 


res2_PPP2R1A_table_top30 <- res2_PPP2R1A_table %>% filter(res2_PPP2R1A_table$padj<0.1 & res2_PPP2R1A_table$log2FoldChange>0) %>% slice_max(n=30, order_by = stat) 
res2_PPP2R1A_table_bottom30 <- res2_PPP2R1A_table %>% filter(res2_PPP2R1A_table$padj<0.1 & res2_PPP2R1A_table$log2FoldChange<0) %>% slice_min(n=30, order_by = stat) 
res2_PPP2R1A_top_bottom_DEGs <- c(rownames(res2_PPP2R1A_table_top30),rownames(res2_PPP2R1A_table_bottom30))



data_PPP2R1A<-log2tpm[immune_related_genes$symbol,
                      c(rownames(sample_table_pre_treatment)[sample_table_pre_treatment$PPP2R1A.mut=="No"],
                        rownames(sample_table_pre_treatment)[sample_table_pre_treatment$PPP2R1A.mut=="Yes"])]

data_PPP2R1A<-t(scale(t(data_PPP2R1A)))
data_PPP2R1A[data_PPP2R1A>1.5]=1.5
data_PPP2R1A[data_PPP2R1A< -1.5]= -1.5

sample_table_pre_treatment_for_anno <- sample_table_pre_treatment
colnames(sample_table_pre_treatment_for_anno)[11]<-"PPP2R1A mut"
sample_table_pre_treatment_for_anno$`PPP2R1A mut` <- as.numeric(sample_table_pre_treatment_for_anno$`PPP2R1A mut`)
sample_table_pre_treatment_for_anno$`PPP2R1A mut`[sample_table_pre_treatment_for_anno$`PPP2R1A mut`== 0]<-"No"
sample_table_pre_treatment_for_anno$`PPP2R1A mut`[sample_table_pre_treatment_for_anno$`PPP2R1A mut`== 1]<-"Yes"

pheatmap::pheatmap(data_PPP2R1A,show_colnames =T,show_rownames =T, fontsize_row = 7,cluster_rows = F, cluster_cols = F,
                   color=colorRampPalette(c('blue','black','yellow'))(100),
                   annotation_col=sample_table_pre_treatment_for_anno[,c('PPP2R1A.mut'),drop=FALSE])

########################################## *ED Fig 4a GSEA pre-treatment  PPP2R1Amut AKTwt vs. PPP2R1Awt AKTwt #########################################################

res2_PPP2R1A_no_alt_table_filter <- filter(res2_PPP2R1A_no_alt_table, !is.na(res2_PPP2R1A_no_alt_table$stat))
Pre_PPP2R1A_Mut_vs_no_alt_genes <- res2_PPP2R1A_no_alt_table_filter %>% arrange(desc(stat)) %>% dplyr::select(Gene,stat)
Pre_PPP2R1A_Mut_vs_no_alt_genes <- tibble::deframe(Pre_PPP2R1A_Mut_vs_no_alt_genes)

TLS_signatue <- data.frame(gs_name="TLS", gene_symbol=c("CD19","CD1D","CD22","CD3D","CD3E","CD52","CD74","CD79A",
                                                        "CD79B","CD8A","CD8B","CR2","CXCL13","CXCR5","FCER2","IL7R",
                                                        "ITGAE","MS4A1","PDCD1","PTGDS","CD4","TRBC2"))
hallmark <- rbind(hallmark, TLS_signatue)
hallmark <- hallmark[-c(8210:8231),]

set.seed(10)
Pre_PPP2R1A_Mut_vs_no_alt_gsea <- GSEA(Pre_PPP2R1A_Mut_vs_no_alt_genes,TERM2GENE=hallmark, eps=0, pvalueCutoff = 1)
Pre_PPP2R1A_Mut_vs_no_alt_gsea <- as.data.frame(Pre_PPP2R1A_Mut_vs_no_alt_gsea)

Pre_PPP2R1A_Mut_vs_no_alt_gsea <- Pre_PPP2R1A_Mut_vs_no_alt_gsea%>%arrange(desc(NES), .by_group = F)
Pre_PPP2R1A_Mut_vs_no_alt_gsea$ID<-factor(Pre_PPP2R1A_Mut_vs_no_alt_gsea$ID,levels=rev(Pre_PPP2R1A_Mut_vs_no_alt_gsea$ID))
Pre_PPP2R1A_Mut_vs_no_alt_gsea$`-log10(FDR)`<- -log10(Pre_PPP2R1A_Mut_vs_no_alt_gsea$qvalue)

p<-ggplot(Pre_PPP2R1A_Mut_vs_no_alt_gsea, aes(ID, NES)) +
  geom_col(aes(fill=`-log10(FDR)`)) +
  #scale_fill_viridis() +
  #scale_fill_gradientn(colors=jet()[c(1:11,15:25)],na.value = "#e6e6e6") +
  coord_flip() +
  labs(x="", y="NES", title="")+
  theme(axis.text.x = element_text(size = 12,angle=0,vjust=0.2,hjust=0.95),
        strip.text.x = element_text(size=12, color="black", face="bold"),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

ggsave('On_vs_Pre.rm.outliers.NES.vs.FDR.barplot.pdf',p,height=10,width=10)

Pre_PPP2R1A_Mut_vs_no_alt_gsea$Color <- NA
Pre_PPP2R1A_Mut_vs_no_alt_gsea$Color[Pre_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.05 & Pre_PPP2R1A_Mut_vs_no_alt_gsea$NES >0] <- "FDR < 0.05, pos"
Pre_PPP2R1A_Mut_vs_no_alt_gsea$Color[Pre_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.01 & Pre_PPP2R1A_Mut_vs_no_alt_gsea$NES >0] <- "FDR < 0.01, pos"
Pre_PPP2R1A_Mut_vs_no_alt_gsea$Color[Pre_PPP2R1A_Mut_vs_no_alt_gsea$qvalue< 0.001 & Pre_PPP2R1A_Mut_vs_no_alt_gsea$NES >0] <- "FDR < 0.001, pos"
Pre_PPP2R1A_Mut_vs_no_alt_gsea$Color[Pre_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.0001 & Pre_PPP2R1A_Mut_vs_no_alt_gsea$NES >0] <- "FDR < 0.0001, pos"
Pre_PPP2R1A_Mut_vs_no_alt_gsea$Color[Pre_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.05 & Pre_PPP2R1A_Mut_vs_no_alt_gsea$NES < 0] <- "FDR < 0.05, neg"
Pre_PPP2R1A_Mut_vs_no_alt_gsea$Color[Pre_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.01 & Pre_PPP2R1A_Mut_vs_no_alt_gsea$NES < 0] <- "FDR < 0.01, neg"
Pre_PPP2R1A_Mut_vs_no_alt_gsea$Color[Pre_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.001 & Pre_PPP2R1A_Mut_vs_no_alt_gsea$NES < 0] <- "FDR < 0.001, neg"
Pre_PPP2R1A_Mut_vs_no_alt_gsea$Color[Pre_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.0001 & Pre_PPP2R1A_Mut_vs_no_alt_gsea$NES < 0] <- "FDR < 0.0001, neg"
Pre_PPP2R1A_Mut_vs_no_alt_gsea$Color <- factor(Pre_PPP2R1A_Mut_vs_no_alt_gsea$Color, levels = c("NS","FDR < 0.05, pos","FDR < 0.01, pos","FDR < 0.001, pos","FDR < 0.0001, pos",
                                                                                        "FDR < 0.05, neg","FDR < 0.01, neg","FDR < 0.001, neg","FDR < 0.0001, neg"))


ggplot(Pre_PPP2R1A_Mut_vs_no_alt_gsea,aes(x=NES, y=`-log10(FDR)`, color=Color)) +
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
  ggtitle("Pre-treatment PPP2R1A Mut vs. No alt GSEA") +
  #ylim(0,150) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "right") +
  geom_hline(yintercept=-log10(0.05), lty = "dashed")+
  labs(x = "Enriched in WT <- NES -> Enriched in Mut", y = "-log10(FDR)", color = "Significance", size="|NES|") +
  #geom_text_repel(data = head(subset(result, gn.selected),30),aes(label = rownames(head(subset(result, gn.selected),30))),size = 2,
  geom_text_repel(data = Pre_PPP2R1A_Mut_vs_no_alt_gsea[Pre_PPP2R1A_Mut_vs_no_alt_gsea$`-log10(FDR)`>-log10(0.05),],
                  aes(label = gsub('HALLMARK_','', rownames(Pre_PPP2R1A_Mut_vs_no_alt_gsea[Pre_PPP2R1A_Mut_vs_no_alt_gsea$`-log10(FDR)`>-log10(0.05),]))),size = 3,
                  box.padding = unit(0.35, "lines"),  point.padding = unit(0.3, "lines"), color="black", max.overlaps = 100) 


ggsave(paste0('On_vs_Pre.rm.outliers.NES.vs.FDR.',Sys.Date(),'.pdf'),plotx,width = 6,height = 6)

write.csv(Pre_PPP2R1A_Mut_vs_no_alt_gsea,"EDFig4a_Pre_PPP2R1A_Mut_vs_no_alt_gsea.csv")








########################################## DEG pre-treatment  AKT alt vs. no alt #########################################################

sample_table_pre_treatment_filter_AKT <- filter(sample_table_pre_treatment, sample_table_pre_treatment$PPP2R1A.mut=="No")
resultcoding_count_pre_treatment_filter_AKT <- resultcoding_count_pre_treatment[, rownames(sample_table_pre_treatment_filter_AKT)]
dds_AKT_no_alt <- DESeqDataSetFromMatrix(countData = resultcoding_count_pre_treatment_filter_AKT,
                                             colData = sample_table_pre_treatment_filter_AKT,
                                             design = ~ AKT.alt)
dds_AKT_no_alt <- DESeq(dds_AKT_no_alt)
res_AKT_no_alt <- results(dds_AKT_no_alt)

#res2_PPP2R1A_shrink <- lfcShrink(dds_PPP2R1A, coef="PPP2R1A.mut_Yes_vs_No", type="ashr")
#res2_PPP2R1A_shrink_table <- as.data.frame(res2_PPP2R1A_shrink)
#res2_PPP2R1A
#summary(res2_PPP2R1A)

res2_AKT_no_alt_table <- as.data.frame(res_AKT_no_alt)
res2_AKT_no_alt_table$invert_P <- (-log10(res2_AKT_no_alt_table$padj)) 
res2_AKT_no_alt_table$Gene=rownames(res2_AKT_no_alt_table)
res2_AKT_no_alt_table$Color <- "NS"
res2_AKT_no_alt_table$Color[res2_AKT_no_alt_table$padj < 0.1] <- "FDR < 0.1"
res2_AKT_no_alt_table$Color[res2_AKT_no_alt_table$padj < 0.01] <- "FDR < 0.01"
res2_AKT_no_alt_table$Color[res2_AKT_no_alt_table$padj< 0.001] <- "FDR < 0.001"
res2_AKT_no_alt_table$Color[res2_AKT_no_alt_table$padj < 0.0001] <- "FDR < 0.0001"
res2_AKT_no_alt_table$Color <- factor(res2_AKT_no_alt_table$Color, levels = c("NS","FDR < 0.1","FDR < 0.01","FDR < 0.001","FDR < 0.0001"))

ggplot(res2_AKT_no_alt_table, aes(x = log2FoldChange, y = invert_P, color = Color,label = Gene)) +
  geom_vline(xintercept = c(2, -2), lty = "dashed") +
  geom_hline(yintercept = -log10(0.1), lty = "dashed") +
  geom_point(size=0.3) +
  labs(x = "Enriched in WT <- log2(FC) -> Enriched in PPP2R1A Mut",
       y = "-log10(p.adj)",
       color = "Significance") +
  scale_color_manual(values = c(
    `FDR < 0.0001` = "#08306b",
    `FDR < 0.001` = "#2171b5",
    `FDR < 0.01` = "#6baed6",
    `FDR < 0.1` = "#c6dbef",
    `NS` = "#f0f0f0"),
    guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(res2_AKT_no_alt_table, padj < 0.05 & abs(log2FoldChange) > 2),
                  size = 3, point.padding = 0.2, color = "black",
                  min.segment.length = .5, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  ggtitle("DEGs in pre-treatment samples")+
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom") 


res2_PPP2R1A_table_top30 <- res2_AKT_no_alt_table %>% filter(res2_AKT_no_alt_table$padj<0.1 & res2_PPP2R1A_table$log2FoldChange>0) %>% slice_max(n=30, order_by = stat) 
res2_PPP2R1A_table_bottom30 <- res2_AKT_no_alt_table %>% filter(res2_AKT_no_alt_table$padj<0.1 & res2_PPP2R1A_table$log2FoldChange<0) %>% slice_min(n=30, order_by = stat) 
res2_PPP2R1A_top_bottom_DEGs <- c(rownames(res2_PPP2R1A_table_top30),rownames(res2_PPP2R1A_table_bottom30))



data_AKT<-log2tpm[c("PPP2R1A"),
                      c(rownames(sample_table_pre_treatment_filter_AKT)[sample_table_pre_treatment_filter_AKT$AKT.alt=="No"],
                        rownames(sample_table_pre_treatment_filter_AKT)[sample_table_pre_treatment_filter_AKT$AKT.alt=="Yes"])]

data_AKT <- as.data.frame(t(data_AKT))
data_AKT$Sample_ID <- rownames(data_AKT)
data_AKT <- left_join(data_AKT, sample_table_pre_treatment_filter_AKT, by="Sample_ID")

ggplot(data_AKT, aes(x=AKT.alt, y=PPP2R1A))+
  geom_violin(scale = "width")+
  geom_boxplot(width=0.3)+
  theme_classic()

data_AKT<-t(scale(t(data_AKT)))
data_AKT[data_AKT>1.5]=1.5
data_AKT[data_AKT< -1.5]= -1.5

sample_table_pre_treatment_for_anno <- sample_table_pre_treatment
colnames(sample_table_pre_treatment_for_anno)[11]<-"PPP2R1A mut"
sample_table_pre_treatment_for_anno$`PPP2R1A mut` <- as.numeric(sample_table_pre_treatment_for_anno$`PPP2R1A mut`)
sample_table_pre_treatment_for_anno$`PPP2R1A mut`[sample_table_pre_treatment_for_anno$`PPP2R1A mut`== 0]<-"No"
sample_table_pre_treatment_for_anno$`PPP2R1A mut`[sample_table_pre_treatment_for_anno$`PPP2R1A mut`== 1]<-"Yes"

pheatmap::pheatmap(data_AKT,show_colnames =T,show_rownames =T, fontsize_row = 7,cluster_rows = F, cluster_cols = F,
                   color=colorRampPalette(c('blue','black','yellow'))(100),
                   annotation_col=sample_table_pre_treatment_filter_AKT[,c('AKT.alt'),drop=FALSE])

########################################## *ED Fig 8 GSEA pre-treatment  AKT alt vs. no alt #########################################################

res2_AKT_no_alt_table_filter <- filter(res2_AKT_no_alt_table, !is.na(res2_AKT_no_alt_table$stat))
Pre_AKT_alt_vs_no_alt_genes <- res2_AKT_no_alt_table_filter %>% arrange(desc(stat)) %>% dplyr::select(Gene,stat)
Pre_AKT_alt_vs_no_alt_genes <- tibble::deframe(Pre_AKT_alt_vs_no_alt_genes)

TLS_signatue <- data.frame(gs_name="TLS", gene_symbol=c("CD19","CD1D","CD22","CD3D","CD3E","CD52","CD74","CD79A",
                                                        "CD79B","CD8A","CD8B","CR2","CXCL13","CXCR5","FCER2","IL7R",
                                                        "ITGAE","MS4A1","PDCD1","PTGDS","CD4","TRBC2"))
hallmark <- rbind(hallmark, TLS_signatue)
hallmark <- hallmark[-c(8210:8231),]

set.seed(10)
Pre_AKT_alt_vs_no_alt_gsea <- GSEA(Pre_AKT_alt_vs_no_alt_genes,TERM2GENE=hallmark, eps=0, pvalueCutoff = 1)
Pre_AKT_alt_vs_no_alt_gsea <- as.data.frame(Pre_AKT_alt_vs_no_alt_gsea)

Pre_AKT_alt_vs_no_alt_gsea <- Pre_AKT_alt_vs_no_alt_gsea%>%arrange(desc(NES), .by_group = F)
Pre_AKT_alt_vs_no_alt_gsea$ID<-factor(Pre_AKT_alt_vs_no_alt_gsea$ID,levels=rev(Pre_AKT_alt_vs_no_alt_gsea$ID))
Pre_AKT_alt_vs_no_alt_gsea$`-log10(FDR)`<- -log10(Pre_AKT_alt_vs_no_alt_gsea$qvalue)

p<-ggplot(Pre_AKT_alt_vs_no_alt_gsea, aes(ID, NES)) +
  geom_col(aes(fill=`-log10(FDR)`)) +
  #scale_fill_viridis() +
  #scale_fill_gradientn(colors=jet()[c(1:11,15:25)],na.value = "#e6e6e6") +
  coord_flip() +
  labs(x="", y="NES", title="")+
  theme(axis.text.x = element_text(size = 12,angle=0,vjust=0.2,hjust=0.95),
        strip.text.x = element_text(size=12, color="black", face="bold"),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

ggsave('On_vs_Pre.rm.outliers.NES.vs.FDR.barplot.pdf',p,height=10,width=10)

Pre_AKT_alt_vs_no_alt_gsea$Color <- NA
Pre_AKT_alt_vs_no_alt_gsea$Color[Pre_AKT_alt_vs_no_alt_gsea$qvalue < 0.05 & Pre_AKT_alt_vs_no_alt_gsea$NES >0] <- "FDR < 0.05, pos"
Pre_AKT_alt_vs_no_alt_gsea$Color[Pre_AKT_alt_vs_no_alt_gsea$qvalue < 0.01 & Pre_AKT_alt_vs_no_alt_gsea$NES >0] <- "FDR < 0.01, pos"
Pre_AKT_alt_vs_no_alt_gsea$Color[Pre_AKT_alt_vs_no_alt_gsea$qvalue< 0.001 & Pre_AKT_alt_vs_no_alt_gsea$NES >0] <- "FDR < 0.001, pos"
Pre_AKT_alt_vs_no_alt_gsea$Color[Pre_AKT_alt_vs_no_alt_gsea$qvalue < 0.0001 & Pre_AKT_alt_vs_no_alt_gsea$NES >0] <- "FDR < 0.0001, pos"
Pre_AKT_alt_vs_no_alt_gsea$Color[Pre_AKT_alt_vs_no_alt_gsea$qvalue < 0.05 & Pre_AKT_alt_vs_no_alt_gsea$NES < 0] <- "FDR < 0.05, neg"
Pre_AKT_alt_vs_no_alt_gsea$Color[Pre_AKT_alt_vs_no_alt_gsea$qvalue < 0.01 & Pre_AKT_alt_vs_no_alt_gsea$NES < 0] <- "FDR < 0.01, neg"
Pre_AKT_alt_vs_no_alt_gsea$Color[Pre_AKT_alt_vs_no_alt_gsea$qvalue < 0.001 & Pre_AKT_alt_vs_no_alt_gsea$NES < 0] <- "FDR < 0.001, neg"
Pre_AKT_alt_vs_no_alt_gsea$Color[Pre_AKT_alt_vs_no_alt_gsea$qvalue < 0.0001 & Pre_AKT_alt_vs_no_alt_gsea$NES < 0] <- "FDR < 0.0001, neg"
Pre_AKT_alt_vs_no_alt_gsea$Color <- factor(Pre_AKT_alt_vs_no_alt_gsea$Color, levels = c("NS","FDR < 0.05, pos","FDR < 0.01, pos","FDR < 0.001, pos","FDR < 0.0001, pos",
                                                                                                "FDR < 0.05, neg","FDR < 0.01, neg","FDR < 0.001, neg","FDR < 0.0001, neg"))


ggplot(Pre_AKT_alt_vs_no_alt_gsea,aes(x=NES, y=`-log10(FDR)`, color=Color)) +
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
  ggtitle("Pre-treatment AKT alt vs. No alt GSEA") +
  #ylim(0,150) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "right") +
  geom_hline(yintercept=-log10(0.05), lty = "dashed")+
  labs(x = "Enriched in WT <- NES -> Enriched in Alt", y = "-log10(FDR)", color = "Significance", size="|NES|") +
  #geom_text_repel(data = head(subset(result, gn.selected),30),aes(label = rownames(head(subset(result, gn.selected),30))),size = 2,
  geom_text_repel(data = Pre_AKT_alt_vs_no_alt_gsea[Pre_AKT_alt_vs_no_alt_gsea$`-log10(FDR)`>-log10(0.05),],
                  aes(label = gsub('HALLMARK_','', rownames(Pre_AKT_alt_vs_no_alt_gsea[Pre_AKT_alt_vs_no_alt_gsea$`-log10(FDR)`>-log10(0.05),]))),size = 3,
                  box.padding = unit(0.35, "lines"),  point.padding = unit(0.3, "lines"), color="black", max.overlaps = 100) 


ggsave(paste0('On_vs_Pre.rm.outliers.NES.vs.FDR.',Sys.Date(),'.pdf'),plotx,width = 6,height = 6)

write.csv(Pre_AKT_alt_vs_no_alt_gsea,"EDFig8_Pre_AKT_alt_vs_no_alt_gsea.csv")








########################################## *ED Fig 6b DEG pre-treatment  PPP2R1A mut, good vs. bad  #########################################################

sample_table_pre_treatment_mut_PPP2R1A <- filter(sample_table_pre_treatment, sample_table_pre_treatment$PPP2R1A.mut=="Yes")
resultcoding_count_pre_treatment_PPP2R1A <- resultcoding_count[,rownames(sample_table_pre_treatment_mut_PPP2R1A)]

dds_pre_PPP2R1A <- DESeqDataSetFromMatrix(countData = resultcoding_count_pre_treatment_PPP2R1A,
                                      colData = sample_table_pre_treatment_mut_PPP2R1A,
                                      design = ~ Response)
dds_pre_PPP2R1A <- DESeq(dds_pre_PPP2R1A)
res_pre_PPP2R1A <- results(dds_pre_PPP2R1A)
#res_pre_PPP2R1A
#resultsNames(dds_pre_PPP2R1A)

res2_pre_PPP2R1A_shrink <- lfcShrink(dds_pre_PPP2R1A, coef="Response_GOOD_vs_BAD", type="ashr")
#res2_pre_PPP2R1A
#summary(res2_pre_PPP2R1A)

res2_pre_PPP2R1A_table <- as.data.frame(res_pre_PPP2R1A)
res2_pre_PPP2R1A_table$invert_P <- (-log10(res2_pre_PPP2R1A_table$padj)) 
res2_pre_PPP2R1A_table$Gene=rownames(res2_pre_PPP2R1A_table)

res2_pre_PPP2R1A_table_shrink <- as.data.frame(res2_pre_PPP2R1A_shrink)
res2_pre_PPP2R1A_table_shrink$invert_P <- (-log10(res2_pre_PPP2R1A_table_shrink$padj)) 
res2_pre_PPP2R1A_table_shrink$Gene=rownames(res2_pre_PPP2R1A_table_shrink)
res2_pre_PPP2R1A_table_shrink$Color <- "NS"
res2_pre_PPP2R1A_table_shrink$Color[res2_pre_PPP2R1A_table_shrink$padj < 0.1 & res2_pre_PPP2R1A_table_shrink$log2FoldChange >0] <- "FDR < 0.1, pos"
res2_pre_PPP2R1A_table_shrink$Color[res2_pre_PPP2R1A_table_shrink$padj < 0.05 & res2_pre_PPP2R1A_table_shrink$log2FoldChange >0] <- "FDR < 0.05, pos"
res2_pre_PPP2R1A_table_shrink$Color[res2_pre_PPP2R1A_table_shrink$padj < 0.01 & res2_pre_PPP2R1A_table_shrink$log2FoldChange >0] <- "FDR < 0.01, pos"
res2_pre_PPP2R1A_table_shrink$Color[res2_pre_PPP2R1A_table_shrink$padj < 0.001 & res2_pre_PPP2R1A_table_shrink$log2FoldChange >0] <- "FDR < 0.001, pos"
res2_pre_PPP2R1A_table_shrink$Color[res2_pre_PPP2R1A_table_shrink$padj < 0.1 & res2_pre_PPP2R1A_table_shrink$log2FoldChange < 0] <- "FDR < 0.1, neg"
res2_pre_PPP2R1A_table_shrink$Color[res2_pre_PPP2R1A_table_shrink$padj < 0.05 & res2_pre_PPP2R1A_table_shrink$log2FoldChange < 0] <- "FDR < 0.05, neg"
res2_pre_PPP2R1A_table_shrink$Color[res2_pre_PPP2R1A_table_shrink$padj < 0.01 & res2_pre_PPP2R1A_table_shrink$log2FoldChange < 0] <- "FDR < 0.01, neg"
res2_pre_PPP2R1A_table_shrink$Color[res2_pre_PPP2R1A_table_shrink$padj < 0.001 & res2_pre_PPP2R1A_table_shrink$log2FoldChange < 0] <- "FDR < 0.001, neg"
res2_pre_PPP2R1A_table_shrink$Color <- factor(res2_pre_PPP2R1A_table_shrink$Color, levels = c("NS","FDR < 0.1, pos","FDR < 0.05, pos","FDR < 0.01, pos","FDR < 0.001, pos",
                                                                                                    "FDR < 0.1, neg","FDR < 0.05, neg","FDR < 0.01, neg","FDR < 0.001, neg"))


ggplot(res2_pre_PPP2R1A_table_shrink, aes(x = log2FoldChange, y = invert_P, color = Color,label = Gene)) +
  #geom_vline(xintercept = c(2, -2), lty = "dashed") +
  #geom_hline(yintercept = -log10(0.1), lty = "dashed") +
  geom_point(size=0.3) +
  labs(x = "Enriched in bad responders <- log2(FC) -> Enriched in good responders",
       y = "-log10(p.adj)",
       color = "Significance") +
  scale_color_manual(values = c(
    `FDR < 0.001, pos` = "#bd0026",
    `FDR < 0.01, pos` = "#fc4e2a",
    `FDR < 0.05, pos` = "#feb24c",
    `FDR < 0.1, pos` = "#ffeda0",
    `FDR < 0.001, neg` = "#08306b",
    `FDR < 0.01, neg` = "#2171b5",
    `FDR < 0.05, neg` = "#6baed6",
    `FDR < 0.1, neg` = "#c6dbef",
    `NS` = "#f0f0f0"),
    guide = guide_legend(override.aes = list(size = 4))) +
  scale_x_continuous(limits = c(-10,11),expand = expansion(mult = c(0,0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(res2_pre_PPP2R1A_table_shrink, Gene %in% immune_gene_highlight_pre),
                  size = 3, point.padding = 0.2, color = "black",
                  min.segment.length = 0, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  ggtitle("DEGs in PPP2R1A mut patients, pre-treatment samples")+
  theme_bw(base_size = 10) +
  theme(legend.position = "right") 

immune_gene_highlight_pre <- c("CD69","CD1C","FCER2","CCR7","BACH2","GPR183","CD70")

write.csv(res2_pre_PPP2R1A_table_shrink,"EDFig6b_res2_pre_PPP2R1A_table_shrink.csv")

res2_pre_PPP2R1A_table_top50 <- res2_pre_PPP2R1A_table %>% filter(res2_pre_PPP2R1A_table$padj<0.1 & res2_pre_PPP2R1A_table$log2FoldChange>0) %>% slice_max(n=50, order_by = log2FoldChange) 
res2_pre_PPP2R1A_table_bottom50 <- res2_pre_PPP2R1A_table %>% filter(res2_pre_PPP2R1A_table$padj<0.1 & res2_pre_PPP2R1A_table$log2FoldChange<0) %>% slice_min(n=50, order_by = log2FoldChange) 
res2_pre_PPP2R1A_top_bottom_DEGs <- c(res2_pre_PPP2R1A_table_top50$Gene,res2_pre_PPP2R1A_table_bottom50$Gene)

genes=intersect(res2_pre_PPP2R1A_top_bottom_DEGs, rownames(log2tpm))
data_PPP2R1A<-log2tpm[genes,
                      c(rownames(sample_table_pre_treatment_mut_PPP2R1A)[sample_table_pre_treatment_mut_PPP2R1A$Response=="BAD"],
                        rownames(sample_table_pre_treatment_mut_PPP2R1A)[sample_table_pre_treatment_mut_PPP2R1A$Response=="GOOD"])]

data_PPP2R1A<-t(scale(t(data_PPP2R1A)))
data_PPP2R1A[data_PPP2R1A>1.5]=1.5
data_PPP2R1A[data_PPP2R1A< -1.5]= -1.5

pheatmap::pheatmap(data_PPP2R1A,show_colnames =T,show_rownames =T, fontsize_row = 6,cluster_rows = F, cluster_cols = F,
                   color=colorRampPalette(c('blue','black','yellow'))(100),
                   annotation_col=sample_table_pre_treatment_mut_PPP2R1A[,c('Response'),drop=FALSE])

########################################## *ED Fig 6a GSEA pre-treatment  PPP2R1A mut, good vs. bad #########################################################

res2_pre_PPP2R1A_table_filter <- filter(res2_pre_PPP2R1A_table, !is.na(res2_pre_PPP2R1A_table$stat))
Pre_PPP2R1A_Mut_response_genes <- res2_pre_PPP2R1A_table_filter %>% arrange(desc(stat)) %>% dplyr::select(Gene,stat)
Pre_PPP2R1A_Mut_response_genes <- tibble::deframe(Pre_PPP2R1A_Mut_response_genes)

set.seed(42)
Pre_PPP2R1A_Mut_response_gsea <- GSEA(Pre_PPP2R1A_Mut_response_genes,TERM2GENE=hallmark, eps=0, pvalueCutoff = 1)
Pre_PPP2R1A_Mut_response_gsea <- as.data.frame(Pre_PPP2R1A_Mut_response_gsea)

Pre_PPP2R1A_Mut_response_gsea <- Pre_PPP2R1A_Mut_response_gsea%>%arrange(desc(NES), .by_group = F)
Pre_PPP2R1A_Mut_response_gsea$ID<-factor(Pre_PPP2R1A_Mut_response_gsea$ID,levels=rev(Pre_PPP2R1A_Mut_response_gsea$ID))
Pre_PPP2R1A_Mut_response_gsea$`-log10(FDR)`<- -log10(Pre_PPP2R1A_Mut_response_gsea$qvalue)

p<-ggplot(Pre_PPP2R1A_Mut_response_gsea, aes(ID, NES)) +
  geom_col(aes(fill=`-log10(FDR)`)) +
  #scale_fill_viridis() +
  #scale_fill_gradientn(colors=jet()[c(1:11,15:25)],na.value = "#e6e6e6") +
  coord_flip() +
  labs(x="", y="NES", title="")+
  theme(axis.text.x = element_text(size = 12,angle=0,vjust=0.2,hjust=0.95),
        strip.text.x = element_text(size=12, color="black", face="bold"),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

ggsave('On_vs_Pre.rm.outliers.NES.vs.FDR.barplot.pdf',p,height=10,width=10)

Pre_PPP2R1A_Mut_response_gsea$Color <- NA
Pre_PPP2R1A_Mut_response_gsea$Color[Pre_PPP2R1A_Mut_response_gsea$qvalue < 0.05 & Pre_PPP2R1A_Mut_response_gsea$NES >0] <- "FDR < 0.05, pos"
Pre_PPP2R1A_Mut_response_gsea$Color[Pre_PPP2R1A_Mut_response_gsea$qvalue < 0.01 & Pre_PPP2R1A_Mut_response_gsea$NES >0] <- "FDR < 0.01, pos"
Pre_PPP2R1A_Mut_response_gsea$Color[Pre_PPP2R1A_Mut_response_gsea$qvalue< 0.001 & Pre_PPP2R1A_Mut_response_gsea$NES >0] <- "FDR < 0.001, pos"
Pre_PPP2R1A_Mut_response_gsea$Color[Pre_PPP2R1A_Mut_response_gsea$qvalue < 0.0001 & Pre_PPP2R1A_Mut_response_gsea$NES >0] <- "FDR < 0.0001, pos"
Pre_PPP2R1A_Mut_response_gsea$Color[Pre_PPP2R1A_Mut_response_gsea$qvalue < 0.05 & Pre_PPP2R1A_Mut_response_gsea$NES < 0] <- "FDR < 0.05, neg"
Pre_PPP2R1A_Mut_response_gsea$Color[Pre_PPP2R1A_Mut_response_gsea$qvalue < 0.01 & Pre_PPP2R1A_Mut_response_gsea$NES < 0] <- "FDR < 0.01, neg"
Pre_PPP2R1A_Mut_response_gsea$Color[Pre_PPP2R1A_Mut_response_gsea$qvalue < 0.001 & Pre_PPP2R1A_Mut_response_gsea$NES < 0] <- "FDR < 0.001, neg"
Pre_PPP2R1A_Mut_response_gsea$Color[Pre_PPP2R1A_Mut_response_gsea$qvalue < 0.0001 & Pre_PPP2R1A_Mut_response_gsea$NES < 0] <- "FDR < 0.0001, neg"
Pre_PPP2R1A_Mut_response_gsea$Color <- factor(Pre_PPP2R1A_Mut_response_gsea$Color, levels = c("NS","FDR < 0.05, pos","FDR < 0.01, pos","FDR < 0.001, pos","FDR < 0.0001, pos",
                                                                                        "FDR < 0.05, neg","FDR < 0.01, neg","FDR < 0.001, neg","FDR < 0.0001, neg"))


ggplot(Pre_PPP2R1A_Mut_response_gsea,aes(x=NES, y=`-log10(FDR)`, color=Color)) +
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
  ggtitle("Pre-treatment PPP2R1A Mut, good vs. bad responders,  GSEA") +
  #ylim(0,150) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "right") +
  geom_hline(yintercept=-log10(0.1), lty = "dashed")+
  labs(x = "Enriched in bad responders <- NES -> Enriched in good responders", y = "-log10(FDR)", color = "Significance", size="|NES|") +
  #geom_text_repel(data = head(subset(result, gn.selected),30),aes(label = rownames(head(subset(result, gn.selected),30))),size = 2,
  geom_text_repel(data = Pre_PPP2R1A_Mut_response_gsea[Pre_PPP2R1A_Mut_response_gsea$`-log10(FDR)`>-log10(0.1),],
                  aes(label = gsub('HALLMARK_','', rownames(Pre_PPP2R1A_Mut_response_gsea[Pre_PPP2R1A_Mut_response_gsea$`-log10(FDR)`>-log10(0.1),]))),size = 3,
                  box.padding = unit(0.35, "lines"),  point.padding = unit(0.3, "lines"), color="black", max.overlaps = 100) 


ggsave(paste0('On_vs_Pre.rm.outliers.NES.vs.FDR.',Sys.Date(),'.pdf'),plotx,width = 6,height = 6)

write.csv(Pre_PPP2R1A_Mut_response_gsea,"EDFig6a_Pre_PPP2R1A_Mut_response_gsea.csv")










########################################## DEG on-treatment  PPP2R1A mut vs. WT #########################################################

#sample_table_on_treatment_filter <- filter(sample_table_on_treatment, sample_table_on_treatment$AKT.alt=="No")
#resultcoding_count_on_treatment_filter <- resultcoding_count_on_treatment[,rownames(sample_table_on_treatment_filter)]
dds_on_PPP2R1A <- DESeqDataSetFromMatrix(countData = resultcoding_count_on_treatment,
                                         colData = sample_table_on_treatment,
                                         design = ~ PPP2R1A.mut)
dds_on_PPP2R1A <- DESeq(dds_on_PPP2R1A)
res_on_PPP2R1A <- results(dds_on_PPP2R1A)

res_on_PPP2R1A
resultsNames(dds_on_PPP2R1A)

#res2_on_PPP2R1A <- lfcShrink(dds_on_PPP2R1A, coef="PPP2R1A.mut_Yes_vs_No", type="ashr")
#res_on_PPP2R1A
#summary(res_on_PPP2R1A)

res2_on_PPP2R1A_table <- as.data.frame(res_on_PPP2R1A)
res2_on_PPP2R1A_table$invert_P <- (-log10(res2_on_PPP2R1A_table$padj)) 
res2_on_PPP2R1A_table$Gene=rownames(res2_on_PPP2R1A_table)
res2_on_PPP2R1A_table$Color <- "NS"
res2_on_PPP2R1A_table$Color[res2_on_PPP2R1A_table$padj < 0.1] <- "FDR < 0.1"
res2_on_PPP2R1A_table$Color[res2_on_PPP2R1A_table$padj < 0.05] <- "FDR < 0.05"
res2_on_PPP2R1A_table$Color[res2_on_PPP2R1A_table$padj< 0.01] <- "FDR < 0.01"
res2_on_PPP2R1A_table$Color[res2_on_PPP2R1A_table$padj < 0.001] <- "FDR < 0.001"
res2_on_PPP2R1A_table$Color <- factor(res2_on_PPP2R1A_table$Color, levels = c("NS","FDR < 0.1","FDR < 0.05","FDR < 0.01","FDR < 0.001"))

ggplot(res2_on_PPP2R1A_table, aes(x = log2FoldChange, y = invert_P, color = Color,label = Gene)) +
  geom_vline(xintercept = c(1, -1), lty = "dashed") +
  geom_hline(yintercept = -log10(0.1), lty = "dashed") +
  geom_point(size=0.3) +
  labs(x = "Enriched in WT <- log2(FC) -> Enriched in PPP2R1A Mut",
       y = "-log10(p.adj)",
       color = "Significance") +
  scale_color_manual(values = c(
    `FDR < 0.001` = "#08306b",
    `FDR < 0.01` = "#2171b5",
    `FDR < 0.05` = "#6baed6",
    `FDR < 0.1` = "#c6dbef",
    `NS` = "#f0f0f0"),
    guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(res2_on_PPP2R1A_table, padj < 0.05 & abs(log2FoldChange) > 2),
                  size = 3, point.padding = 0.2, color = "black",
                  min.segment.length = .5, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  ggtitle("DEGs in on-treatment samples")+
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom") 


res2_on_PPP2R1A_table_top30 <- res2_on_PPP2R1A_table %>% filter(res2_on_PPP2R1A_table$padj<0.2 & res2_on_PPP2R1A_table$log2FoldChange>0)  %>% slice_max(n=30, order_by = stat) 
res2_on_PPP2R1A_table_bottom30 <- res2_on_PPP2R1A_table %>% filter(res2_on_PPP2R1A_table$padj<0.2 & res2_on_PPP2R1A_table$log2FoldChange<0) %>% slice_min(n=30, order_by = stat) 
res2_on_PPP2R1A_on_top_bottom_DEGs <- c(res2_on_PPP2R1A_table_top30$Gene,res2_on_PPP2R1A_table_bottom30$Gene)


data_on_PPP2R1A<-log2tpm[immune_related_genes$symbol,
                         c(rownames(sample_table_on_treatment)[sample_table_on_treatment$PPP2R1A.mut=="No"],
                           rownames(sample_table_on_treatment)[sample_table_on_treatment$PPP2R1A.mut=="Yes"])]

data_on_PPP2R1A<-t(scale(t(data_on_PPP2R1A)))
data_on_PPP2R1A[data_on_PPP2R1A>1.5]=1.5
data_on_PPP2R1A[data_on_PPP2R1A< -1.5]= -1.5

pheatmap::pheatmap(data_on_PPP2R1A,show_colnames =T,show_rownames =T, fontsize_row = 7,cluster_rows = F, cluster_cols = F,
                   color=colorRampPalette(c('blue','black','yellow'))(100),
                   annotation_col = sample_table_on_treatment[,c("PPP2R1A.mut"),drop=FALSE])



########################################## *Fig 3b GSEA on-treatment  PPP2R1A mut vs. WT #########################################################

res2_on_PPP2R1A_table <- filter(res2_on_PPP2R1A_table, !is.na(res2_on_PPP2R1A_table$stat))
On_PPP2R1A_Mut_vs_WT_genes <- res2_on_PPP2R1A_table %>% arrange(desc(stat)) %>% dplyr::select(Gene,stat)
On_PPP2R1A_Mut_vs_WT_genes <- tibble::deframe(On_PPP2R1A_Mut_vs_WT_genes)

set.seed(13)
On_PPP2R1A_Mut_vs_WT_gsea <- GSEA(On_PPP2R1A_Mut_vs_WT_genes,TERM2GENE=hallmark, eps=0, pvalueCutoff = 1)
On_PPP2R1A_Mut_vs_WT_gsea <- as.data.frame(On_PPP2R1A_Mut_vs_WT_gsea)

On_PPP2R1A_Mut_vs_WT_gsea <- On_PPP2R1A_Mut_vs_WT_gsea%>%arrange(desc(NES), .by_group = F)
On_PPP2R1A_Mut_vs_WT_gsea$ID<-factor(On_PPP2R1A_Mut_vs_WT_gsea$ID,levels=rev(On_PPP2R1A_Mut_vs_WT_gsea$ID))
On_PPP2R1A_Mut_vs_WT_gsea$`-log10(FDR)`<- -log10(On_PPP2R1A_Mut_vs_WT_gsea$qvalue)

p<-ggplot(On_PPP2R1A_Mut_vs_WT_gsea, aes(ID, NES)) +
  geom_col(aes(fill=`-log10(FDR)`)) +
  #scale_fill_viridis() +
  #scale_fill_gradientn(colors=jet()[c(1:11,15:25)],na.value = "#e6e6e6") +
  coord_flip() +
  labs(x="", y="NES", title="")+
  theme(axis.text.x = element_text(size = 12,angle=0,vjust=0.2,hjust=0.95),
        strip.text.x = element_text(size=12, color="black", face="bold"),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

ggsave('On_vs_Pre.rm.outliers.NES.vs.FDR.barplot.pdf',p,height=10,width=10)


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
  geom_hline(yintercept=-log10(0.1), lty = "dashed")+
  labs(x = "Enriched in WT <- NES -> Enriched in Mut", y = "-log10(FDR)", color = "Significance", size="|NES|") +
  #geom_text_repel(data = head(subset(result, gn.selected),30),aes(label = rownames(head(subset(result, gn.selected),30))),size = 2,
  geom_text_repel(data = On_PPP2R1A_Mut_vs_WT_gsea[On_PPP2R1A_Mut_vs_WT_gsea$`-log10(FDR)`>-log10(0.05),],
                  aes(label = gsub('HALLMARK_','', rownames(On_PPP2R1A_Mut_vs_WT_gsea[On_PPP2R1A_Mut_vs_WT_gsea$`-log10(FDR)`>-log10(0.05),]))),size = 3,
                  box.padding = unit(0.35, "lines"),  point.padding = unit(0.3, "lines"), color="black", max.overlaps = 100) 


ggsave(paste0('On_vs_Pre.rm.outliers.NES.vs.FDR.',Sys.Date(),'.pdf'),plotx,width = 6,height = 6)

write.csv(On_PPP2R1A_Mut_vs_WT_gsea,"Fig3b_On_PPP2R1A_Mut_vs_WT_gsea.csv")














########################################## DEG on-treatment  PPP2R1Amut AKTwt vs. PPP2R1Awt AKTwt #########################################################

sample_table_on_treatment_filter <- filter(sample_table_on_treatment, sample_table_on_treatment$AKT.alt=="No")
resultcoding_count_on_treatment_filter <- resultcoding_count_on_treatment[,rownames(sample_table_on_treatment_filter)]
dds_on_PPP2R1A_no_alt <- DESeqDataSetFromMatrix(countData = resultcoding_count_on_treatment_filter,
                                                colData = sample_table_on_treatment_filter,
                                                design = ~ PPP2R1A.mut)
dds_on_PPP2R1A_no_alt <- DESeq(dds_on_PPP2R1A_no_alt)
res_on_PPP2R1A_no_alt <- results(dds_on_PPP2R1A_no_alt)


#res2_on_PPP2R1A <- lfcShrink(dds_on_PPP2R1A, coef="PPP2R1A.mut_Yes_vs_No", type="ashr")
#res_on_PPP2R1A
#summary(res_on_PPP2R1A)

res2_on_PPP2R1A_no_alt_table <- as.data.frame(res_on_PPP2R1A_no_alt)
res2_on_PPP2R1A_no_alt_table$invert_P <- (-log10(res2_on_PPP2R1A_no_alt_table$padj)) 
res2_on_PPP2R1A_no_alt_table$Gene=rownames(res2_on_PPP2R1A_no_alt_table)
res2_on_PPP2R1A_no_alt_table$Color <- "NS"
res2_on_PPP2R1A_no_alt_table$Color[res2_on_PPP2R1A_no_alt_table$padj < 0.1] <- "FDR < 0.1"
res2_on_PPP2R1A_no_alt_table$Color[res2_on_PPP2R1A_no_alt_table$padj < 0.05] <- "FDR < 0.05"
res2_on_PPP2R1A_no_alt_table$Color[res2_on_PPP2R1A_no_alt_table$padj< 0.01] <- "FDR < 0.01"
res2_on_PPP2R1A_no_alt_table$Color[res2_on_PPP2R1A_no_alt_table$padj < 0.001] <- "FDR < 0.001"
res2_on_PPP2R1A_no_alt_table$Color <- factor(res2_on_PPP2R1A_no_alt_table$Color, levels = c("NS","FDR < 0.1","FDR < 0.05","FDR < 0.01","FDR < 0.001"))

ggplot(res2_on_PPP2R1A_no_alt_table, aes(x = log2FoldChange, y = invert_P, color = Color,label = Gene)) +
  geom_vline(xintercept = c(1, -1), lty = "dashed") +
  geom_hline(yintercept = -log10(0.1), lty = "dashed") +
  geom_point(size=0.3) +
  labs(x = "Enriched in WT <- log2(FC) -> Enriched in PPP2R1A Mut",
       y = "-log10(p.adj)",
       color = "Significance") +
  scale_color_manual(values = c(
    `FDR < 0.001` = "#08306b",
    `FDR < 0.01` = "#2171b5",
    `FDR < 0.05` = "#6baed6",
    `FDR < 0.1` = "#c6dbef",
    `NS` = "#f0f0f0"),
    guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(res2_on_PPP2R1A_no_alt_table, padj < 0.05 & abs(log2FoldChange) > 2),
                  size = 3, point.padding = 0.2, color = "black",
                  min.segment.length = .5, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  ggtitle("DEGs in on-treatment samples")+
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom") 


res2_on_PPP2R1A_table_top30 <- res2_on_PPP2R1A_table %>% filter(res2_on_PPP2R1A_table$padj<0.2 & res2_on_PPP2R1A_table$log2FoldChange>0)  %>% slice_max(n=30, order_by = stat) 
res2_on_PPP2R1A_table_bottom30 <- res2_on_PPP2R1A_table %>% filter(res2_on_PPP2R1A_table$padj<0.2 & res2_on_PPP2R1A_table$log2FoldChange<0) %>% slice_min(n=30, order_by = stat) 
res2_on_PPP2R1A_on_top_bottom_DEGs <- c(res2_on_PPP2R1A_table_top30$Gene,res2_on_PPP2R1A_table_bottom30$Gene)


data_on_PPP2R1A<-log2tpm[immune_related_genes$symbol,
                         c(rownames(sample_table_on_treatment)[sample_table_on_treatment$PPP2R1A.mut=="No"],
                           rownames(sample_table_on_treatment)[sample_table_on_treatment$PPP2R1A.mut=="Yes"])]

data_on_PPP2R1A<-t(scale(t(data_on_PPP2R1A)))
data_on_PPP2R1A[data_on_PPP2R1A>1.5]=1.5
data_on_PPP2R1A[data_on_PPP2R1A< -1.5]= -1.5

pheatmap::pheatmap(data_on_PPP2R1A,show_colnames =T,show_rownames =T, fontsize_row = 7,cluster_rows = F, cluster_cols = F,
                   color=colorRampPalette(c('blue','black','yellow'))(100),
                   annotation_col = sample_table_on_treatment[,c("PPP2R1A.mut"),drop=FALSE])



########################################## *ED Fig 4b GSEA on-treatment  PPP2R1Amut AKTwt vs. PPP2R1Awt AKTwt #########################################################

res2_on_PPP2R1A_no_alt_table_filter <- filter(res2_on_PPP2R1A_no_alt_table, !is.na(res2_on_PPP2R1A_no_alt_table$stat))
On_PPP2R1A_Mut_vs_no_alt_genes <- res2_on_PPP2R1A_no_alt_table_filter %>% arrange(desc(stat)) %>% dplyr::select(Gene,stat)
On_PPP2R1A_Mut_vs_no_alt_genes <- tibble::deframe(On_PPP2R1A_Mut_vs_no_alt_genes)

set.seed(13)
On_PPP2R1A_Mut_vs_no_alt_gsea <- GSEA(On_PPP2R1A_Mut_vs_no_alt_genes,TERM2GENE=hallmark, eps=0, pvalueCutoff = 1)
On_PPP2R1A_Mut_vs_no_alt_gsea <- as.data.frame(On_PPP2R1A_Mut_vs_no_alt_gsea)

On_PPP2R1A_Mut_vs_no_alt_gsea <- On_PPP2R1A_Mut_vs_no_alt_gsea%>%arrange(desc(NES), .by_group = F)
On_PPP2R1A_Mut_vs_no_alt_gsea$ID<-factor(On_PPP2R1A_Mut_vs_no_alt_gsea$ID,levels=rev(On_PPP2R1A_Mut_vs_no_alt_gsea$ID))
On_PPP2R1A_Mut_vs_no_alt_gsea$`-log10(FDR)`<- -log10(On_PPP2R1A_Mut_vs_no_alt_gsea$qvalue)

p<-ggplot(On_PPP2R1A_Mut_vs_no_alt_gsea, aes(ID, NES)) +
  geom_col(aes(fill=`-log10(FDR)`)) +
  #scale_fill_viridis() +
  #scale_fill_gradientn(colors=jet()[c(1:11,15:25)],na.value = "#e6e6e6") +
  coord_flip() +
  labs(x="", y="NES", title="")+
  theme(axis.text.x = element_text(size = 12,angle=0,vjust=0.2,hjust=0.95),
        strip.text.x = element_text(size=12, color="black", face="bold"),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

ggsave('On_vs_Pre.rm.outliers.NES.vs.FDR.barplot.pdf',p,height=10,width=10)


On_PPP2R1A_Mut_vs_no_alt_gsea$Color <- NA
On_PPP2R1A_Mut_vs_no_alt_gsea$Color[On_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.05 & On_PPP2R1A_Mut_vs_no_alt_gsea$NES >0] <- "FDR < 0.05, pos"
On_PPP2R1A_Mut_vs_no_alt_gsea$Color[On_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.01 & On_PPP2R1A_Mut_vs_no_alt_gsea$NES >0] <- "FDR < 0.01, pos"
On_PPP2R1A_Mut_vs_no_alt_gsea$Color[On_PPP2R1A_Mut_vs_no_alt_gsea$qvalue< 0.001 & On_PPP2R1A_Mut_vs_no_alt_gsea$NES >0] <- "FDR < 0.001, pos"
On_PPP2R1A_Mut_vs_no_alt_gsea$Color[On_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.0001 & On_PPP2R1A_Mut_vs_no_alt_gsea$NES >0] <- "FDR < 0.0001, pos"
On_PPP2R1A_Mut_vs_no_alt_gsea$Color[On_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.05 & On_PPP2R1A_Mut_vs_no_alt_gsea$NES < 0] <- "FDR < 0.05, neg"
On_PPP2R1A_Mut_vs_no_alt_gsea$Color[On_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.01 & On_PPP2R1A_Mut_vs_no_alt_gsea$NES < 0] <- "FDR < 0.01, neg"
On_PPP2R1A_Mut_vs_no_alt_gsea$Color[On_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.001 & On_PPP2R1A_Mut_vs_no_alt_gsea$NES < 0] <- "FDR < 0.001, neg"
On_PPP2R1A_Mut_vs_no_alt_gsea$Color[On_PPP2R1A_Mut_vs_no_alt_gsea$qvalue < 0.0001 & On_PPP2R1A_Mut_vs_no_alt_gsea$NES < 0] <- "FDR < 0.0001, neg"
On_PPP2R1A_Mut_vs_no_alt_gsea$Color <- factor(On_PPP2R1A_Mut_vs_no_alt_gsea$Color, levels = c("NS","FDR < 0.05, pos","FDR < 0.01, pos","FDR < 0.001, pos","FDR < 0.0001, pos",
                                                                                      "FDR < 0.05, neg","FDR < 0.01, neg","FDR < 0.001, neg","FDR < 0.0001, neg"))

ggplot(On_PPP2R1A_Mut_vs_no_alt_gsea,aes(x=NES, y=`-log10(FDR)`, color=Color)) +
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
  ggtitle("On-treatment PPP2R1A Mut vs. No alt GSEA") +
  #ylim(0,150) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "right") +
  geom_hline(yintercept=-log10(0.1), lty = "dashed")+
  labs(x = "Enriched in No alt <- NES -> Enriched in Mut", y = "-log10(FDR)", color = "Significance", size="|NES|") +
  #geom_text_repel(data = head(subset(result, gn.selected),30),aes(label = rownames(head(subset(result, gn.selected),30))),size = 2,
  geom_text_repel(data = On_PPP2R1A_Mut_vs_no_alt_gsea[On_PPP2R1A_Mut_vs_no_alt_gsea$`-log10(FDR)`>-log10(0.05),],
                  aes(label = gsub('HALLMARK_','', rownames(On_PPP2R1A_Mut_vs_no_alt_gsea[On_PPP2R1A_Mut_vs_no_alt_gsea$`-log10(FDR)`>-log10(0.05),]))),size = 3,
                  box.padding = unit(0.35, "lines"),  point.padding = unit(0.3, "lines"), color="black", max.overlaps = 100) 


ggsave(paste0('On_vs_Pre.rm.outliers.NES.vs.FDR.',Sys.Date(),'.pdf'),plotx,width = 6,height = 6)

write.csv(On_PPP2R1A_Mut_vs_no_alt_gsea,"EDFig4b_On_PPP2R1A_Mut_vs_no_alt_gsea.csv")













########################################## *ED Fig 6d DEG on-treatment  PPP2R1A mut, good vs. poor responder #########################################################

sample_table_on_treatment_PPP2R1A <- filter(sample_table_on_treatment, sample_table_on_treatment$PPP2R1A.mut=="Yes")
resultcoding_count_on_treatment_PPP2R1A <- resultcoding_count_on_treatment[,rownames(sample_table_on_treatment_PPP2R1A)]
dds_on_PPP2R1A_mut <- DESeqDataSetFromMatrix(countData = resultcoding_count_on_treatment_PPP2R1A,
                                             colData = sample_table_on_treatment_PPP2R1A,
                                             design = ~ Response)
dds_on_PPP2R1A_mut <- DESeq(dds_on_PPP2R1A_mut)
res_on_PPP2R1A_mut <- results(dds_on_PPP2R1A_mut)
res_on_PPP2R1A_mut
resultsNames(dds_on_PPP2R1A_mut)

res2_on_PPP2R1A_mut_shrink <- lfcShrink(dds_on_PPP2R1A_mut, coef="Response_GOOD_vs_BAD", type="ashr")
#res2_on_PPP2R1A_mut
#summary(res2_on_PPP2R1A_mut)

res2_on_PPP2R1A_mut_table <- as.data.frame(res_on_PPP2R1A_mut)
res2_on_PPP2R1A_mut_table$invert_P <- (-log10(res2_on_PPP2R1A_mut_table$padj)) 
res2_on_PPP2R1A_mut_table$Gene=rownames(res2_on_PPP2R1A_mut_table)

res2_on_PPP2R1A_mut_table_shrink <- as.data.frame(res2_on_PPP2R1A_mut_shrink)
res2_on_PPP2R1A_mut_table_shrink$invert_P <- (-log10(res2_on_PPP2R1A_mut_table_shrink$padj)) 
res2_on_PPP2R1A_mut_table_shrink$Gene=rownames(res2_on_PPP2R1A_mut_table_shrink)
res2_on_PPP2R1A_mut_table_shrink$Color <- "NS"
res2_on_PPP2R1A_mut_table_shrink$Color[res2_on_PPP2R1A_mut_table_shrink$padj < 0.1 & res2_on_PPP2R1A_mut_table_shrink$log2FoldChange >0] <- "FDR < 0.1, pos"
res2_on_PPP2R1A_mut_table_shrink$Color[res2_on_PPP2R1A_mut_table_shrink$padj < 0.05 & res2_on_PPP2R1A_mut_table_shrink$log2FoldChange >0] <- "FDR < 0.05, pos"
res2_on_PPP2R1A_mut_table_shrink$Color[res2_on_PPP2R1A_mut_table_shrink$padj < 0.01 & res2_on_PPP2R1A_mut_table_shrink$log2FoldChange >0] <- "FDR < 0.01, pos"
res2_on_PPP2R1A_mut_table_shrink$Color[res2_on_PPP2R1A_mut_table_shrink$padj < 0.001 & res2_on_PPP2R1A_mut_table_shrink$log2FoldChange >0] <- "FDR < 0.001, pos"
res2_on_PPP2R1A_mut_table_shrink$Color[res2_on_PPP2R1A_mut_table_shrink$padj < 0.1 & res2_on_PPP2R1A_mut_table_shrink$log2FoldChange < 0] <- "FDR < 0.1, neg"
res2_on_PPP2R1A_mut_table_shrink$Color[res2_on_PPP2R1A_mut_table_shrink$padj < 0.05 & res2_on_PPP2R1A_mut_table_shrink$log2FoldChange < 0] <- "FDR < 0.05, neg"
res2_on_PPP2R1A_mut_table_shrink$Color[res2_on_PPP2R1A_mut_table_shrink$padj < 0.01 & res2_on_PPP2R1A_mut_table_shrink$log2FoldChange < 0] <- "FDR < 0.01, neg"
res2_on_PPP2R1A_mut_table_shrink$Color[res2_on_PPP2R1A_mut_table_shrink$padj < 0.001 & res2_on_PPP2R1A_mut_table_shrink$log2FoldChange < 0] <- "FDR < 0.001, neg"
res2_on_PPP2R1A_mut_table_shrink$Color <- factor(res2_on_PPP2R1A_mut_table_shrink$Color, levels = c("NS","FDR < 0.1, pos","FDR < 0.05, pos","FDR < 0.01, pos","FDR < 0.001, pos",
                                                                                      "FDR < 0.1, neg","FDR < 0.05, neg","FDR < 0.01, neg","FDR < 0.001, neg"))


ggplot(res2_on_PPP2R1A_mut_table_shrink, aes(x = log2FoldChange, y = invert_P, color = Color,label = Gene)) +
  geom_vline(xintercept = c(1, -1), lty = "dashed") +
  geom_hline(yintercept = -log10(0.1), lty = "dashed") +
  geom_point(size=0.1) +
  labs(x = "Enriched in bad responders <- log2(FC) -> Enriched in good responders",
       y = "-log10(p.adj)",
       color = "Significance") +
  scale_color_manual(values = c(
    `FDR < 0.001, pos` = "#bd0026",
    `FDR < 0.01, pos` = "#fc4e2a",
    `FDR < 0.05, pos` = "#feb24c",
    `FDR < 0.1, pos` = "#ffeda0",
    `FDR < 0.001, neg` = "#08306b",
    `FDR < 0.01, neg` = "#2171b5",
    `FDR < 0.05, neg` = "#6baed6",
    `FDR < 0.1, neg` = "#c6dbef",
    `NS` = "#f0f0f0"),
    guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  scale_x_continuous(limits = c(-8,8),expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(res2_on_PPP2R1A_mut_table_shrink, Gene %in% imune_genes_highlight),
                  size = 3, point.padding = 0, color = "black",
                  min.segment.length = 0, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  ggtitle("DEGs in PPP2R1A mut patients, on-treatment samples")+
  theme_bw(base_size = 10) +
  theme(legend.position = "right") 

write.csv(res2_on_PPP2R1A_mut_table_shrink,"EDFig6d_res2_on_PPP2R1A_mut_table_shrink.csv")

res2_on_PPP2R1A_mut_table_top100 <- res2_on_PPP2R1A_mut_table %>% filter(res2_on_PPP2R1A_mut_table$padj<0.1 & res2_on_PPP2R1A_mut_table$log2FoldChange>0)  %>% slice_max(n=100, order_by = log2FoldChange) 
res2_on_PPP2R1A_mut_table_bottom30 <- res2_on_PPP2R1A_mut_table %>% filter(res2_on_PPP2R1A_mut_table$padj<0.1 & res2_on_PPP2R1A_mut_table$log2FoldChange<0) %>% slice_min(n=30, order_by = log2FoldChange) 
res2_on_PPP2R1A_mut_top_bottom_DEGs <- c(res2_on_PPP2R1A_mut_table_top30$Gene,res2_on_PPP2R1A_mut_table_bottom30$Gene)

immune_genes = intersect(res2_on_PPP2R1A_mut_table_top100$Gene, immune_related_genes$symbol)
imune_genes_highlight <- c("HLA-DQB2","HLA-DQA2","HLA-DPA1","HLA-DRA","HLA-DRB1","CD74","CD3E","CD4","IL2RA","IL2RG",
                           "CD79A","CD79B","GZMB","CCL2","CCL3","CCL4","CD38","BTK","TNFSF13B","CD86","CXCL13","PDCD1","BATF")

data_on_PPP2R1A_mut<-log2tpm[genes,
                             c(rownames(sample_table_on_treatment_PPP2R1A)[sample_table_on_treatment_PPP2R1A$Response=="BAD"],
                               rownames(sample_table_on_treatment_PPP2R1A)[sample_table_on_treatment_PPP2R1A$Response=="GOOD"])]

data_on_PPP2R1A_mut<-t(scale(t(data_on_PPP2R1A_mut)))
data_on_PPP2R1A_mut[data_on_PPP2R1A_mut>1.5]=1.5
data_on_PPP2R1A_mut[data_on_PPP2R1A_mut< -1.5]= -1.5

pheatmap::pheatmap(data_on_PPP2R1A_mut,show_colnames =T,show_rownames =T, fontsize_row = 6,cluster_rows = F, cluster_cols = F,
                   color=colorRampPalette(c('blue','black','yellow'))(100),
                   annotation_col = sample_table_on_treatment_PPP2R1A[,c("Response"),drop=FALSE])



########################################## *ED Fig 6c GSEA on-treatment  PPP2R1A mut, good vs. poor responder #########################################################

res2_on_PPP2R1A_mut_table_filter <- filter(res2_on_PPP2R1A_mut_table, !is.na(res2_on_PPP2R1A_mut_table$stat))
On_PPP2R1A_Mut_response_genes <- res2_on_PPP2R1A_mut_table_filter %>% arrange(desc(stat)) %>% dplyr::select(Gene,stat)
On_PPP2R1A_Mut_response_genes <- tibble::deframe(On_PPP2R1A_Mut_response_genes)

set.seed(42)
On_PPP2R1A_Mut_response_gsea <- GSEA(On_PPP2R1A_Mut_response_genes,TERM2GENE=hallmark, eps=0, pvalueCutoff = 1)
On_PPP2R1A_Mut_response_gsea <- as.data.frame(On_PPP2R1A_Mut_response_gsea)

On_PPP2R1A_Mut_response_gsea <- On_PPP2R1A_Mut_response_gsea%>%arrange(desc(NES), .by_group = F)
On_PPP2R1A_Mut_response_gsea$ID<-factor(On_PPP2R1A_Mut_response_gsea$ID,levels=rev(On_PPP2R1A_Mut_response_gsea$ID))
On_PPP2R1A_Mut_response_gsea$`-log10(FDR)`<- -log10(On_PPP2R1A_Mut_response_gsea$qvalue)

p<-ggplot(On_PPP2R1A_Mut_response_gsea, aes(ID, NES)) +
  geom_col(aes(fill=`-log10(FDR)`)) +
  #scale_fill_viridis() +
  #scale_fill_gradientn(colors=jet()[c(1:11,15:25)],na.value = "#e6e6e6") +
  coord_flip() +
  labs(x="", y="NES", title="")+
  theme(axis.text.x = element_text(size = 12,angle=0,vjust=0.2,hjust=0.95),
        strip.text.x = element_text(size=12, color="black", face="bold"),
        text = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

ggsave('On_vs_Pre.rm.outliers.NES.vs.FDR.barplot.pdf',p,height=10,width=10)

On_PPP2R1A_Mut_response_gsea$Color <- "NS"
On_PPP2R1A_Mut_response_gsea$Color[On_PPP2R1A_Mut_response_gsea$qvalue < 0.1] <- "FDR < 0.1"
On_PPP2R1A_Mut_response_gsea$Color[On_PPP2R1A_Mut_response_gsea$qvalue < 0.05] <- "FDR < 0.05"
On_PPP2R1A_Mut_response_gsea$Color[On_PPP2R1A_Mut_response_gsea$qvalue< 0.01] <- "FDR < 0.01"
On_PPP2R1A_Mut_response_gsea$Color[On_PPP2R1A_Mut_response_gsea$qvalue < 0.001] <- "FDR < 0.001"
On_PPP2R1A_Mut_response_gsea$Color <- factor(On_PPP2R1A_Mut_response_gsea$Color, levels = c("NS","FDR < 0.1","FDR < 0.05","FDR < 0.01","FDR < 0.001"))

On_PPP2R1A_Mut_response_gsea$Color <- NA
On_PPP2R1A_Mut_response_gsea$Color[On_PPP2R1A_Mut_response_gsea$qvalue < 0.05 & On_PPP2R1A_Mut_response_gsea$NES >0] <- "FDR < 0.05, pos"
On_PPP2R1A_Mut_response_gsea$Color[On_PPP2R1A_Mut_response_gsea$qvalue < 0.01 & On_PPP2R1A_Mut_response_gsea$NES >0] <- "FDR < 0.01, pos"
On_PPP2R1A_Mut_response_gsea$Color[On_PPP2R1A_Mut_response_gsea$qvalue< 0.001 & On_PPP2R1A_Mut_response_gsea$NES >0] <- "FDR < 0.001, pos"
On_PPP2R1A_Mut_response_gsea$Color[On_PPP2R1A_Mut_response_gsea$qvalue < 0.0001 & On_PPP2R1A_Mut_response_gsea$NES >0] <- "FDR < 0.0001, pos"
On_PPP2R1A_Mut_response_gsea$Color[On_PPP2R1A_Mut_response_gsea$qvalue < 0.05 & On_PPP2R1A_Mut_response_gsea$NES < 0] <- "FDR < 0.05, neg"
On_PPP2R1A_Mut_response_gsea$Color[On_PPP2R1A_Mut_response_gsea$qvalue < 0.01 & On_PPP2R1A_Mut_response_gsea$NES < 0] <- "FDR < 0.01, neg"
On_PPP2R1A_Mut_response_gsea$Color[On_PPP2R1A_Mut_response_gsea$qvalue < 0.001 & On_PPP2R1A_Mut_response_gsea$NES < 0] <- "FDR < 0.001, neg"
On_PPP2R1A_Mut_response_gsea$Color[On_PPP2R1A_Mut_response_gsea$qvalue < 0.0001 & On_PPP2R1A_Mut_response_gsea$NES < 0] <- "FDR < 0.0001, neg"
On_PPP2R1A_Mut_response_gsea$Color <- factor(On_PPP2R1A_Mut_response_gsea$Color, levels = c("NS","FDR < 0.05, pos","FDR < 0.01, pos","FDR < 0.001, pos","FDR < 0.0001, pos",
                                                                                              "FDR < 0.05, neg","FDR < 0.01, neg","FDR < 0.001, neg","FDR < 0.0001, neg"))


ggplot(On_PPP2R1A_Mut_response_gsea,aes(x=NES, y=`-log10(FDR)`, color=Color)) +
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
  ggtitle("On-treatment PPP2R1A mut, Good vs. bad responders, GSEA") +
  #ylim(0,150) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "right") +
  #geom_hline(yintercept=-log10(0.1), lty = "dashed")+
  labs(x = "Enriched in bad responders <- NES -> Enriched in good responders", y = "-log10(FDR)", color = "Significance", size="|NES|") +
  #geom_text_repel(data = head(subset(result, gn.selected),30),aes(label = rownames(head(subset(result, gn.selected),30))),size = 2,
  geom_text_repel(data = On_PPP2R1A_Mut_response_gsea[On_PPP2R1A_Mut_response_gsea$`-log10(FDR)`>-log10(0.05),],
                  aes(label = gsub('HALLMARK_','', rownames(On_PPP2R1A_Mut_response_gsea[On_PPP2R1A_Mut_response_gsea$`-log10(FDR)`>-log10(0.05),]))),size = 3,
                  box.padding = unit(0.35, "lines"),  point.padding = unit(0.3, "lines"), color="black", max.overlaps = 100) 


ggsave(paste0('On_vs_Pre.rm.outliers.NES.vs.FDR.',Sys.Date(),'.pdf'),plotx,width = 6,height = 6)

write.csv(On_PPP2R1A_Mut_response_gsea,"EDFig6c_On_PPP2R1A_Mut_response_gsea.csv")




########################################## DEG  PPP2R1A mut on vs. pre #########################################################
resultcoding_count_mut_PPP2R1A <- resultcoding_count[, rownames(sample_table)[sample_table$PPP2R1Amutationpresent==1]]
resultcoding_count_wt_PPP2R1A <- resultcoding_count[, rownames(sample_table)[sample_table$PPP2R1Amutationpresent==0]]
sample_table_mut_PPP2R1A <- filter(sample_table, sample_table$PPP2R1Amutationpresent==1)
sample_table_wt_PPP2R1A <- filter(sample_table, sample_table$PPP2R1Amutationpresent==0)
sample_table_mut_PPP2R1A$Timepoint <- factor(sample_table_mut_PPP2R1A$Timepoint, levels = c("Pre","On"))

dds_mut_PPP2R1A <- DESeqDataSetFromMatrix(countData = resultcoding_count_mut_PPP2R1A,
                                          colData = sample_table_mut_PPP2R1A,
                                          design = ~ Timepoint)
dds_mut_PPP2R1A <- DESeq(dds_mut_PPP2R1A)
res_mut_PPP2R1A <- results(dds_mut_PPP2R1A)
  
resultsNames(dds_mut_PPP2R1A)
res_mut_PPP2R1A <- lfcShrink(dds_mut_PPP2R1A, coef="Timepoint_On_vs_Pre", type="ashr")
summary(res_mut_PPP2R1A)

res_mut_PPP2R1A_table <- as.data.frame(res_mut_PPP2R1A)
res_mut_PPP2R1A_table$invert_P <- (-log10(res_mut_PPP2R1A_table$padj)) 
res_mut_PPP2R1A_table$Gene=rownames(res_mut_PPP2R1A_table)
res_mut_PPP2R1A_table$Color <- "NS"
res_mut_PPP2R1A_table$Color[res_mut_PPP2R1A_table$padj < 0.1] <- "FDR < 0.1"
res_mut_PPP2R1A_table$Color[res_mut_PPP2R1A_table$padj < 0.01] <- "FDR < 0.01"
res_mut_PPP2R1A_table$Color[res_mut_PPP2R1A_table$padj< 0.001] <- "FDR < 0.001"
res_mut_PPP2R1A_table$Color[res_mut_PPP2R1A_table$padj < 0.0001] <- "FDR < 0.0001"
res_mut_PPP2R1A_table$Color <- factor(res_mut_PPP2R1A_table$Color, levels = c("NS","FDR < 0.1","FDR < 0.01","FDR < 0.001","FDR < 0.0001"))

ggplot(res_mut_PPP2R1A_table, aes(x = log2FoldChange, y = invert_P, color = Color,label = Gene)) +
  geom_vline(xintercept = c(2, -2), lty = "dashed") +
  geom_hline(yintercept = -log10(0.1), lty = "dashed") +
  geom_point(size=0.3) +
  labs(x = "Enriched in pre-treatment <- log2(FC) -> Enriched in on-treatment",
       y = "-log10(p.adj)",
       color = "Significance") +
  scale_color_manual(values = c(
    `FDR < 0.0001` = "#08306b",
    `FDR < 0.001` = "#2171b5",
    `FDR < 0.01` = "#6baed6",
    `FDR < 0.1` = "#c6dbef",
    `NS` = "#f0f0f0"),
    guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(res_mut_PPP2R1A_table, padj < 0.05 & abs(log2FoldChange) > 2),
                  size = 3, point.padding = 0.2, color = "black",
                  min.segment.length = .5, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  ggtitle("DEGs in PPP2R1A mut samples")+
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom") 


res_mut_PPP2R1A_table_top30 <- res_mut_PPP2R1A_table  %>% filter(res_mut_PPP2R1A_table$padj<0.1 & res_mut_PPP2R1A_table$log2FoldChange>0) %>% slice_max(n=30, order_by = stat) 
res_mut_PPP2R1A_table_bottom30 <- res_mut_PPP2R1A_table %>% filter(res_mut_PPP2R1A_table$padj<0.1 & res_mut_PPP2R1A_table$log2FoldChange<0) %>% slice_min(n=30, order_by = stat) 
res_mut_PPP2R1A_top_bottom_DEGs <- c(res_mut_PPP2R1A_table_top30$Gene,res_mut_PPP2R1A_table_bottom30$Gene)

data_mut_PPP2R1A<-log2tpm[immune_related_genes$symbol,c(rownames(sample_table_mut_PPP2R1A)[sample_table_mut_PPP2R1A$Timepoint=="Pre"],
                                                            rownames(sample_table_mut_PPP2R1A)[sample_table_mut_PPP2R1A$Timepoint=="On"])]

data_mut_PPP2R1A<-t(scale(t(data_mut_PPP2R1A)))
data_mut_PPP2R1A[data_mut_PPP2R1A>1.5]=1.5
data_mut_PPP2R1A[data_mut_PPP2R1A< -1.5]= -1.5

pheatmap::pheatmap(t(data_mut_PPP2R1A),show_colnames =T,show_rownames =T, fontsize_row = 7,cluster_rows = F, cluster_cols = F,
                   color=colorRampPalette(c('blue','black','yellow'))(100),
                   annotation_row=sample_table_mut[,c('Timepoint'),drop=FALSE])

########################################## *ED Fig 4c GSEA PPP2R1A mut on vs. Pre #########################################################
res_mut_PPP2R1A_table <- filter(res_mut_PPP2R1A_table, !is.na(res_mut_PPP2R1A_table$stat))
PPP2R1A_Mut_On_vs_Pre_genes <- res_mut_PPP2R1A_table %>% arrange(desc(stat)) %>% dplyr::select(Gene,stat)
PPP2R1A_Mut_On_vs_Pre_genes <- tibble::deframe(PPP2R1A_Mut_On_vs_Pre_genes)

set.seed(15)
PPP2R1A_Mut_On_vs_Pre_gsea <- GSEA(PPP2R1A_Mut_On_vs_Pre_genes,TERM2GENE=hallmark, eps=0, pvalueCutoff = 1)
PPP2R1A_Mut_On_vs_Pre_gsea <- as.data.frame(PPP2R1A_Mut_On_vs_Pre_gsea)

PPP2R1A_Mut_On_vs_Pre_gsea <- PPP2R1A_Mut_On_vs_Pre_gsea%>%arrange(desc(NES), .by_group = F)
PPP2R1A_Mut_On_vs_Pre_gsea$ID<-factor(PPP2R1A_Mut_On_vs_Pre_gsea$ID,levels=rev(PPP2R1A_Mut_On_vs_Pre_gsea$ID))
PPP2R1A_Mut_On_vs_Pre_gsea$`-log10(FDR)`<- -log10(PPP2R1A_Mut_On_vs_Pre_gsea$qvalue)

PPP2R1A_Mut_On_vs_Pre_gsea$Color <- NA
PPP2R1A_Mut_On_vs_Pre_gsea$Color[PPP2R1A_Mut_On_vs_Pre_gsea$qvalue < 0.05 & PPP2R1A_Mut_On_vs_Pre_gsea$NES >0] <- "FDR < 0.05, pos"
PPP2R1A_Mut_On_vs_Pre_gsea$Color[PPP2R1A_Mut_On_vs_Pre_gsea$qvalue < 0.01 & PPP2R1A_Mut_On_vs_Pre_gsea$NES >0] <- "FDR < 0.01, pos"
PPP2R1A_Mut_On_vs_Pre_gsea$Color[PPP2R1A_Mut_On_vs_Pre_gsea$qvalue< 0.001 & PPP2R1A_Mut_On_vs_Pre_gsea$NES >0] <- "FDR < 0.001, pos"
PPP2R1A_Mut_On_vs_Pre_gsea$Color[PPP2R1A_Mut_On_vs_Pre_gsea$qvalue < 0.0001 & PPP2R1A_Mut_On_vs_Pre_gsea$NES >0] <- "FDR < 0.0001, pos"
PPP2R1A_Mut_On_vs_Pre_gsea$Color[PPP2R1A_Mut_On_vs_Pre_gsea$qvalue < 0.05 & PPP2R1A_Mut_On_vs_Pre_gsea$NES < 0] <- "FDR < 0.05, neg"
PPP2R1A_Mut_On_vs_Pre_gsea$Color[PPP2R1A_Mut_On_vs_Pre_gsea$qvalue < 0.01 & PPP2R1A_Mut_On_vs_Pre_gsea$NES < 0] <- "FDR < 0.01, neg"
PPP2R1A_Mut_On_vs_Pre_gsea$Color[PPP2R1A_Mut_On_vs_Pre_gsea$qvalue < 0.001 & PPP2R1A_Mut_On_vs_Pre_gsea$NES < 0] <- "FDR < 0.001, neg"
PPP2R1A_Mut_On_vs_Pre_gsea$Color[PPP2R1A_Mut_On_vs_Pre_gsea$qvalue < 0.0001 & PPP2R1A_Mut_On_vs_Pre_gsea$NES < 0] <- "FDR < 0.0001, neg"
PPP2R1A_Mut_On_vs_Pre_gsea$Color <- factor(PPP2R1A_Mut_On_vs_Pre_gsea$Color, levels = c("NS","FDR < 0.05, pos","FDR < 0.01, pos","FDR < 0.001, pos","FDR < 0.0001, pos",
                                                                                        "FDR < 0.05, neg","FDR < 0.01, neg","FDR < 0.001, neg","FDR < 0.0001, neg"))


ggplot(PPP2R1A_Mut_On_vs_Pre_gsea,aes(x=NES, y=`-log10(FDR)`, color=Color)) +
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
  ggtitle("On-treatment vs. pre-treatment in PPP2R1A mut GSEA") +
  #ylim(0,150) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "right") +
  geom_hline(yintercept=-log10(0.05), lty = "dashed")+
  labs(x = "Enriched in pre-treatment <- NES -> Enriched in on-treatment", y = "-log10(FDR)", color = "Significance", size="|NES|") +
  #geom_text_repel(data = head(subset(result, gn.selected),30),aes(label = rownames(head(subset(result, gn.selected),30))),size = 2,
  geom_text_repel(data = PPP2R1A_Mut_On_vs_Pre_gsea[PPP2R1A_Mut_On_vs_Pre_gsea$`-log10(FDR)`>-log10(0.05),],
                  aes(label = gsub('HALLMARK_','', rownames(PPP2R1A_Mut_On_vs_Pre_gsea[PPP2R1A_Mut_On_vs_Pre_gsea$`-log10(FDR)`>-log10(0.05),]))),size = 3,
                  box.padding = unit(0.35, "lines"),  point.padding = unit(0.3, "lines"), color="black", max.overlaps = 100) 


ggsave(paste0('On_vs_Pre.rm.outliers.NES.vs.FDR.',Sys.Date(),'.pdf'),plotx,width = 6,height = 6)

write.csv(PPP2R1A_Mut_On_vs_Pre_gsea,"EDFig4c_PPP2R1A_Mut_On_vs_Pre_gsea.csv")


########################################## Deconvolution #######################################

Cibersort_res <- read.csv("CIBERSORTx_Job6_Results.csv")
colnames(Cibersort_res)[1]<-"Sample_ID"
sample_table$Sample_ID <- rownames(sample_table)
Cibersort_res <- left_join(Cibersort_res, sample_table, by="Sample_ID")
Cibersort_res_filter <- filter(Cibersort_res, Cibersort_res$Sample_ID %notin% c("Sample_On-Acc69", "Sample_On-Acc122")) 
Cibersort_res_filter$Timepoint <- factor(Cibersort_res_filter$Timepoint, levels=c("Pre","On"))

Cibersort_res_filter$Genetic_alterations<-"No alt"
Cibersort_res_filter$Genetic_alterations[Cibersort_res_filter$PPP2R1A.mut=="Yes"]<-"PPP2R1A mut"
Cibersort_res_filter$Genetic_alterations[Cibersort_res_filter$AKT.alt=="Yes"]<-"AKT alt"
Cibersort_res_filter$Genetic_alterations <- factor(Cibersort_res_filter$Genetic_alterations, levels = c("PPP2R1A mut", "AKT alt", "No alt"))

Cibersort_res_mut <- filter(Cibersort_res_filter, Cibersort_res_filter$PPP2R1A.mut.or.AKT.alt=="Yes")
Cibersort_res_wt <- filter(Cibersort_res_filter, Cibersort_res_filter$PPP2R1A.mut.or.AKT.alt=="No")

Cibersort_res_mut_good <- filter(Cibersort_res_mut, Cibersort_res_mut$Response=="GOOD")
Cibersort_res_mut_bad <- filter(Cibersort_res_mut, Cibersort_res_mut$Response=="BAD")

Cibersort_res_filter_paired <- filter(Cibersort_res_filter, 
                                      Cibersort_res_filter$Patient_ID %in% c('Acc160','Acc162','Acc166','Acc119','Acc121','Acc124',
                                                                              'Acc137','Acc139','Acc149','Acc150','Acc19','Acc31',
                                                                             'Acc36','Acc38'))

Cibersort_res_mut_PPP2R1A <- filter(Cibersort_res_filter, Cibersort_res_filter$PPP2R1A.mut == "Yes")
Cibersort_res_wt_PPP2R1A <- filter(Cibersort_res_filter, Cibersort_res_filter$PPP2R1A.mut == "No")

Cibersort_res_mut_PPP2R1A_paired <- filter(Cibersort_res_filter_paired, Cibersort_res_filter_paired$PPP2R1A.mut == "Yes")
Cibersort_res_wt_PPP2R1A_paired <- filter(Cibersort_res_filter_paired, Cibersort_res_filter_paired$PPP2R1A.mut == "No")



Cibersort_res_mut_PPP2R1A_good <- filter(Cibersort_res_mut_PPP2R1A, Cibersort_res_mut_PPP2R1A$Response == "GOOD")
Cibersort_res_mut_PPP2R1A_bad <- filter(Cibersort_res_mut_PPP2R1A, Cibersort_res_mut_PPP2R1A$Response == "BAD")

Cibersort_res_pre <- filter(Cibersort_res, Cibersort_res$Timepoint=="Pre")
Cibersort_res_on <- filter(Cibersort_res, Cibersort_res$Timepoint=="On")

Cibersort_res_pre_mut <- filter(Cibersort_res_pre, Cibersort_res_pre$PPP2R1A.mut.or.AKT.alt=="Yes")
Cibersort_res_pre_mut_PPP2R1A <- filter(Cibersort_res_pre, Cibersort_res_pre$PPP2R1A.mut=="Yes")
Cibersort_res_on_mut <- filter(Cibersort_res_on, Cibersort_res_on$PPP2R1A.mut.or.AKT.alt=="Yes")
Cibersort_res_on_mut_PPP2R1A <- filter(Cibersort_res_on, Cibersort_res_on$PPP2R1A.mut=="Yes")


Cibersort_res_pre$PPP2R1A_AKT_mut <- "No"
Cibersort_res_pre$PPP2R1A_AKT_mut[Cibersort_res_pre$dichotomPPP2R1AorAKTstatus==1]<-"Yes"
Cibersort_res_pre$PPP2R1A_mut <- "No"
Cibersort_res_pre$PPP2R1A_mut[Cibersort_res_pre$PPP2R1Amutationpresent==1]<-"Yes"

Cibersort_res_on$PPP2R1A_AKT_mut <- "No"
Cibersort_res_on$PPP2R1A_AKT_mut[Cibersort_res_on$dichotomPPP2R1AorAKTstatus==1]<-"Yes"
Cibersort_res_on$PPP2R1A_mut <- "No"
Cibersort_res_on$PPP2R1A_mut[Cibersort_res_on$PPP2R1Amutationpresent==1]<-"Yes"


ggplot(Cibersort_res_mut_PPP2R1A_paired, aes(x=Timepoint, y=NK.cells.activated))+
    geom_boxplot(outlier.size = 0,width=0.5, aes(color=Timepoint))+
    geom_point(aes(color=Timepoint), size=5, alpha=1)+
    geom_line(aes(group=Patient_ID))+
    theme_classic()+
    #scale_fill_manual(values=c("#fee391","#ec7014"))+
    scale_color_manual(values=c("#fee391","#ec7014"))+
    theme(axis.ticks.length = unit(0.5,'cm'))+
    theme(legend.position = "none", plot.title = element_text(size=14), axis.title = element_text(size=14), axis.text.x = element_text(angle = 45, hjust=1, size=10))+
    guides(fill=guide_legend(title = NULL))+
    xlab("Timepoint")+
    ylab("Prop. of cells")


