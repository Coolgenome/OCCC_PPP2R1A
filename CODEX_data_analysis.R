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

######################## Comparison of cellular density ###############################################################################
CODEX_data <- read.csv("/Users/ydai6/Library/CloudStorage/OneDrive-InsideMDAnderson/Mac data/OCCC/CODEX/CODEX_Data.csv")

CODEX_data_pre <- filter(CODEX_data, CODEX_data$Baseline.or.Treatment=="Pre-tx")
CODEX_data_pre_tumor <- filter(CODEX_data_pre, CODEX_data_pre$Intratumoral.Compartments=="Tumor")
CODEX_data_pre_stroma <- filter(CODEX_data_pre, CODEX_data_pre$Intratumoral.Compartments=="Stroma")
CODEX_data_pre_total <- filter(CODEX_data_pre, CODEX_data_pre$Intratumoral.Compartments=="Total")

CODEX_data_on <- filter(CODEX_data, CODEX_data$Baseline.or.Treatment=="On-tx")
CODEX_data_on_tumor <- filter(CODEX_data_on, CODEX_data_on$Intratumoral.Compartments=="Tumor")
CODEX_data_on_stroma <- filter(CODEX_data_on, CODEX_data_on$Intratumoral.Compartments=="Stroma")
CODEX_data_on_total <- filter(CODEX_data_on, CODEX_data_on$Intratumoral.Compartments=="Total")


ggplot(CODEX_data_on_stroma, aes(x=PPP2R1A_mut, y=CD45pos.CD56pos..n.mm2.))+
  geom_boxplot(outlier.size = 0,width=0.5, aes(color=PPP2R1A_mut))+
  geom_point(aes(color=PPP2R1A_mut), size=5, alpha=1)+
  #geom_line(aes(group=Acc..))+
  theme_classic()+
  #scale_fill_manual(values=c("#fee391","#ec7014"))+
  scale_color_manual(values=c("#80b1d3","#FDB462"))+
  theme(axis.ticks.length = unit(0.2,'cm'))+
  theme(legend.position = "none", plot.title = element_text(size=12), axis.title = element_text(size=10), axis.text.x = element_text(angle = 45, hjust=1, size=10), 
        axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))+
  guides(fill=guide_legend(title = NULL))+
  xlab("PPP2R1A_mut")+
  ylab("Density of cells")


######################## CODEX neighborhood #############################

cell_position_merge$cell_type <- "Other"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Yes" & cell_position_merge$CD8 == "Yes" & cell_position_merge$CD4 == "Other"] <-"T_CD8_other"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Yes" & cell_position_merge$CD8 == "Yes" & cell_position_merge$CD4 == "Other" & cell_position_merge$CD45RO =="Yes"] <-"T_CD8_CD45RO"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Yes" & cell_position_merge$CD8 == "Yes" & cell_position_merge$CD4 == "Other" & cell_position_merge$`G&B`=="Yes"] <-"T_CD8_GZMB"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Yes" & cell_position_merge$CD8 == "Yes" & cell_position_merge$CD4 == "Other" & cell_position_merge$PD1=="Yes"] <-"T_CD8_PD1"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Yes" & cell_position_merge$CD8 == "Yes" & cell_position_merge$CD4 == "Other" & cell_position_merge$CD39=="Yes"] <-"T_CD8_CD39"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Yes" & cell_position_merge$CD4 == "Yes" & cell_position_merge$CD8 == "Other"] <-"T_CD4_other"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Yes" & cell_position_merge$CD4 == "Yes" & cell_position_merge$CD8 == "Other" & cell_position_merge$OX40=="Yes"] <-"T_CD4_OX40"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Yes" & cell_position_merge$CD4 == "Yes" & cell_position_merge$CD8 == "Other" & cell_position_merge$CTLA4=="Yes"] <-"T_CD4_fh"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Yes" & cell_position_merge$CD4 == "Yes" & cell_position_merge$CD8 == "Other" & cell_position_merge$Foxp3=="Yes"] <-"T_CD4_reg"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD56 == "Yes" & cell_position_merge$CD4 == "Other" & cell_position_merge$CD8 == "Other"] <-"NK"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Other" & cell_position_merge$CD20 == "Yes" ] <-"B"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Other" & cell_position_merge$CD68 == "Yes" & cell_position_merge$CD206 == "Other" & cell_position_merge$HLADR == "Other" ] <-"Mac_other"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Other" & cell_position_merge$CD68 == "Yes" & cell_position_merge$CD206 == "Other" & cell_position_merge$HLADR == "Yes" ] <-"Mac_HLADR"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Other" & cell_position_merge$CD68 == "Yes" & cell_position_merge$CD206 == "Yes" & cell_position_merge$HLADR == "Other" ] <-"Mac_CD206"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Yes" & cell_position_merge$CD3e == "Other" & cell_position_merge$CD68 == "Yes" & cell_position_merge$CD206 == "Yes"& cell_position_merge$HLADR == "Yes"  ] <-"Mac_CD206_HLADR"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Other" & cell_position_merge$CK == "Yes" & cell_position_merge$Ki67 == "Yes" & cell_position_merge$MHCI == "Yes" & cell_position_merge$CD73=="Yes"] <-"Epithelial_prolif_MHCI_CD73"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Other" & cell_position_merge$CK == "Yes" & cell_position_merge$Ki67 == "Yes" & cell_position_merge$MHCI == "Yes" & cell_position_merge$CD73=="Other"] <-"Epithelial_prolif_MHCI"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Other" & cell_position_merge$CK == "Yes" & cell_position_merge$Ki67 == "Yes" & cell_position_merge$MHCI == "Other" & cell_position_merge$CD73=="Yes"] <-"Epithelial_prolif_CD73"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Other" & cell_position_merge$CK == "Yes" & cell_position_merge$Ki67 == "Yes" & cell_position_merge$MHCI == "Other" & cell_position_merge$CD73=="Other"] <-"Epithelial_prolif"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Other" & cell_position_merge$CK == "Yes" & cell_position_merge$Ki67 == "Other" & cell_position_merge$MHCI == "Yes" & cell_position_merge$CD73=="Yes"] <-"Epithelial_nonprolif_MHCI_CD73"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Other" & cell_position_merge$CK == "Yes" & cell_position_merge$Ki67 == "Other" & cell_position_merge$MHCI == "Yes" & cell_position_merge$CD73=="Other"] <-"Epithelial_nonprolif_MHCI"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Other" & cell_position_merge$CK == "Yes" & cell_position_merge$Ki67 == "Other" & cell_position_merge$MHCI == "Other" & cell_position_merge$CD73=="Yes"] <-"Epithelial_nonprolif_CD73"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Other" & cell_position_merge$CK == "Yes" & cell_position_merge$Ki67 == "Other" & cell_position_merge$MHCI == "Other" & cell_position_merge$CD73=="Other"] <-"Epithelial_nonprolif"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Other" & cell_position_merge$CK == "Other" & cell_position_merge$CD31 == "Other"] <-"Stromal_other"
cell_position_merge$cell_type[cell_position_merge$CD45 == "Other" & cell_position_merge$CK == "Other"& cell_position_merge$CD31 == "Yes"] <-"Endothelial"

rownames(cell_position_merge) <- paste0(cell_position_merge$sample_ID,"_",cell_position_merge$Object.ID)

cell_type_order <- c("T_CD8_CD45RO","T_CD8_GZMB","T_CD8_PD1","T_CD8_CD39","T_CD8_other","T_CD4_reg","T_CD4_fh","T_CD4_OX40","T_CD4_other","NK",
                     "B","Mac_CD206","Mac_HLADR","Mac_CD206_HLADR","Mac_other","Endothelial","Stromal_other","Epithelial_prolif","Epithelial_prolif_CD73",
                     "Epithelial_prolif_MHCI","Epithelial_prolif_MHCI_CD73","Epithelial_nonprolif","Epithelial_nonprolif_CD73","Epithelial_nonprolif_MHCI",
                     "Epithelial_nonprolif_MHCI_CD73","Other")

cell_table <- as.data.frame(table(cell_position_merge$cell_type))
cell_table$Var1 <- factor(cell_table$Var1, levels=cell_type_order)
ggplot(cell_table, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  theme_classic()+
  scale_fill_manual(values=c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4", "#46F0F0", "#F032E6", 
                             "#BCF60C", "#FABEBE", "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#AAFFC3", "#FA8072",
                             "#FFD8B1", "#4B0082", "#00FF00", "#FF4500", "#4682B4", "#2E8B57", "#BDB76B", "#B0C4DE",
                             "#80cdc1", "#bebada"))


cell_type_detail = unique(cell_position_merge$cell_type)
my_neighbor_list <- list()
for (clstype in cell_type_detail){
  my.neighbor = c()
  print (clstype)
  
  for(sam in unique(cell_position_merge$sample_ID)) {
    print (sam)
    tmp.meta = cell_position_merge[cell_position_merge$cell_type == clstype & cell_position_merge$sample_ID == sam, ]
    tmp.meta1 = cell_position_merge[cell_position_merge$sample_ID == sam, ]
    
    if(clstype %in% unique(tmp.meta1$cell_type)){
      
      tmp.mtx = matrix(0,nrow = length(cell_type_detail),ncol = nrow(tmp.meta))
      rownames(tmp.mtx) = cell_type_detail
      colnames(tmp.mtx) = rownames(tmp.meta)
      
      for(i in 1:nrow(tmp.meta)) {
        xloc = tmp.meta[i, 'Centroid.X.µm']
        yloc = tmp.meta[i, 'Centroid.Y.µm']
        cutoff = 80
        idx1 = tmp.meta1[, 'Centroid.X.µm'] <=  (xloc + cutoff) & tmp.meta1[,'Centroid.X.µm'] >=  (xloc - cutoff)
        idx2 = tmp.meta1[, 'Centroid.Y.µm'] <=  (yloc + cutoff) & tmp.meta1[,'Centroid.Y.µm'] >=  (yloc - cutoff)
        
        square = tmp.meta1[idx1 & idx2,]
        dis = apply(square, 1, function(x) {sqrt((as.numeric(x['Centroid.X.µm']) - as.numeric(xloc))^2 + (as.numeric(x['Centroid.Y.µm']) - as.numeric(yloc))^2)} )
        
        tmp.freq = table(tmp.meta1[match(names(dis[dis <= cutoff]), rownames(tmp.meta1)), 'cell_type'])
        tmp.mtx[,i] = tmp.freq[match(cell_type_detail,names(tmp.freq))]
      }
      
      tmp.mtx[which(is.na(tmp.mtx))] = 0
      
      my.neighbor = cbind(my.neighbor, tmp.mtx)
    }
  }
  
  my_neighbor_list[[clstype]] <- my.neighbor
  
  write.csv(t(my.neighbor), paste0(clstype,'.neighbor.sub.csv'),quote = F)
}


### statistics ###
neighborhood_df <- c()
for(i in 1: length(my_neighbor_list)){
  name = names(my_neighbor_list)[i]
  df = my_neighbor_list[[i]]
  df <- as.data.frame(t(df))
  df$center_cell_type <- name
  neighborhood_df <- rbind(neighborhood_df, df)
}
neighborhood_df$neighborhood_count_sum <- rowSums(neighborhood_df[,c(1:26)])

neighborhood_df$cell_ID <- rownames(neighborhood_df)
neighborhood_df <- separate(neighborhood_df, col="cell_ID", into=c("C1","C2","C3"),sep="_")
colnames(neighborhood_df)[29] <- "patient_ID"
neighborhood_df$patient_ID <- paste0("Acc",neighborhood_df$patient_ID)
colnames(neighborhood_df)[30] <- "Timepoint"
neighborhood_df$Genetic <- "No"
neighborhood_df$Genetic[neighborhood_df$patient_ID %in% c("Acc31","Acc38","Acc124","Acc149","Acc155","Acc160","Acc162")] <- "PPP2R1A"
neighborhood_df$Genetic[neighborhood_df$patient_ID %in% c("Acc121")] <- "AKT"

neighborhood_df$class <- paste0(neighborhood_df$Genetic,"_", neighborhood_df$Timepoint)
neighborhood_df$class <- factor(neighborhood_df$class, levels = c("PPP2R1A_pre","PPP2R1A_on","AKT_pre","AKT_on","No_pre","No_on"))
neighborhood_df$PPP2R1A <- "No"
neighborhood_df$PPP2R1A[neighborhood_df$Genetic == "PPP2R1A"] <- "Yes"
neighborhood_df$class_PPP2R1A <- paste0(neighborhood_df$PPP2R1A,"_", neighborhood_df$Timepoint)
neighborhood_df$class_PPP2R1A <- factor(neighborhood_df$class_PPP2R1A, levels = c("Yes_pre","Yes_on","No_pre","No_on"))

neighborhood_df_log <- neighborhood_df
neighborhood_df_log[,1:26] <- log10(neighborhood_df_log[,1:26] + 1)

ggplot(filter(neighborhood_df_log, neighborhood_df_log$center_cell_type %in% c("Epithelial_prolif","Epithelial_prolif_CD73","Epithelial_nonprolif","Epithelial_nonprolif_CD73")), 
       aes(x=class_PPP2R1A, y=NK))+
  #geom_violin(scale="width")+
  geom_boxplot(outlier.size = 0,width=0.5, aes(color=class_PPP2R1A))+
  geom_point(aes(color=class_PPP2R1A),size=1)+
  scale_color_manual(values=c("#fee391","#ec7014","#fee391","#ec7014"))+
  theme_classic()+
  theme(axis.ticks.length = unit(0.2,'cm'))+
  theme(legend.position = "none", plot.title = element_text(size=6), axis.title = element_text(size=14), axis.text.x = element_text(angle = 45, hjust=1, size=10),
        axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))+
  guides(fill=guide_legend(title = NULL))
