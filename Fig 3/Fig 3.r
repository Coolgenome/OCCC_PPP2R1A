### Fig. 3 ###

### load required packages ###
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)


### read in CODEX cell density data ###
CODEX_cell_density <- read.csv("./CODEX_cell_density.csv", check.names=F)

### read in patient metadata ###
patient_metadata <- read.csv("./patient_metadata.csv")
patient_metadata$Patient_ID <- paste0("Acc", patient_metadata$Patient_ID)

CODEX_cell_density$Patient_ID <- paste0("Acc", CODEX_cell_density$`Patient ID`)
CODEX_cell_density <- left_join(CODEX_cell_density, patient_metadata, by="Patient_ID")
CODEX_cell_density$PPP2R1A_mutation <- factor(CODEX_cell_density$PPP2R1A_mutation, levels=c("0","1"))


### Fig. 3a was created with BioRender ###


### Fig. 3b-c ###
CODEX_cell_density_pre <- filter(CODEX_cell_density, CODEX_cell_density$Timepoint=="Pre")
CODEX_cell_density_pre_tumor <- filter(CODEX_cell_density_pre, CODEX_cell_density_pre$`Intratumoral Compartments`=="Tumor")
CODEX_cell_density_pre_stroma <- filter(CODEX_cell_density_pre, CODEX_cell_density_pre$`Intratumoral Compartments`=="Stroma")
CODEX_cell_density_pre_total <- filter(CODEX_cell_density_pre, CODEX_cell_density_pre$`Intratumoral Compartments`=="Total")

CODEX_cell_density_on <- filter(CODEX_cell_density, CODEX_cell_density$Timepoint=="On")
CODEX_cell_density_on_tumor <- filter(CODEX_cell_density_on, CODEX_cell_density_on$`Intratumoral Compartments`=="Tumor")
CODEX_cell_density_on_stroma <- filter(CODEX_cell_density_on, CODEX_cell_density_on$`Intratumoral Compartments`=="Stroma")
CODEX_cell_density_on_total <- filter(CODEX_cell_density_on, CODEX_cell_density_on$`Intratumoral Compartments`=="Total")

ggplot(CODEX_cell_density_pre_total, aes(x=PPP2R1A_mutation, y=`CD20+/Ki67+ (n/mm2)`))+
  geom_boxplot(outlier.size = 0,width=0.5, aes(color=PPP2R1A_mutation))+
  geom_point(aes(color=PPP2R1A_mutation), size=5, alpha=1)+
  scale_color_manual(values=c("#80b1d3","#FDB462"))+
  theme_classic()+
  theme(axis.ticks.length = unit(0.2,'cm'))+
  theme(legend.position = "none", plot.title = element_text(size=12), 
        axis.title = element_text(size=10), axis.text.x = element_text(angle = 45, hjust=1, size=10), 
        axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))+
  guides(fill=guide_legend(title = NULL))+
  xlab("PPP2R1A mutation")+
  ylab("Density of CD20+Ki67+ cells")+
  ggtitle("Overall")

### The densities of other cell states or in other tissue compartments are plotted in a similar way. 


### Fig. 3d were from CODEX images ###


### Fig. 3e was created with BioRender ###


### Fig. 3f-g ###
CODEX_cell_expression <- readRDS("./CODEX_cell_expression.rds")
rownames(CODEX_cell_expression) <- CODEX_cell_expression$Cell_ID

CODEX_cell_expression$cell_type <- "Other"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 1 & CODEX_cell_expression$CD3e == 1 & CODEX_cell_expression$CD8 == 1 & CODEX_cell_expression$CD4 == 0] <-"T_CD8_other"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 1 & CODEX_cell_expression$CD3e == 1 & CODEX_cell_expression$CD8 == 1 & CODEX_cell_expression$CD4 == 0 & CODEX_cell_expression$CD45RO ==1] <-"T_CD8_CD45RO"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 1 & CODEX_cell_expression$CD3e == 1 & CODEX_cell_expression$CD8 == 1 & CODEX_cell_expression$CD4 == 0 & CODEX_cell_expression$`G&B`==1] <-"T_CD8_GZMB"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 1 & CODEX_cell_expression$CD3e == 1 & CODEX_cell_expression$CD8 == 1 & CODEX_cell_expression$CD4 == 0 & CODEX_cell_expression$PD1==1] <-"T_CD8_PD1"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 1 & CODEX_cell_expression$CD3e == 1 & CODEX_cell_expression$CD8 == 1 & CODEX_cell_expression$CD4 == 0 & CODEX_cell_expression$CD39==1] <-"T_CD8_CD39"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 1 & CODEX_cell_expression$CD3e == 1 & CODEX_cell_expression$CD4 == 1 & CODEX_cell_expression$CD8 == 0] <-"T_CD4"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 1 & CODEX_cell_expression$CD56 == 1 & CODEX_cell_expression$CD4 == 0 & CODEX_cell_expression$CD8 == 0] <-"NK"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 1 & CODEX_cell_expression$CD3e == 0 & CODEX_cell_expression$CD20 == 1] <-"B"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 1 & CODEX_cell_expression$CD3e == 0 & CODEX_cell_expression$CD68 == 1] <-"Mac"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 0 & CODEX_cell_expression$CK == 1 & CODEX_cell_expression$Ki67 == 1 & CODEX_cell_expression$MHCI == 1 & CODEX_cell_expression$CD73==1] <-"Epithelial_prolif_MHCI_CD73"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 0 & CODEX_cell_expression$CK == 1 & CODEX_cell_expression$Ki67 == 1 & CODEX_cell_expression$MHCI == 1 & CODEX_cell_expression$CD73==0] <-"Epithelial_prolif_MHCI"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 0 & CODEX_cell_expression$CK == 1 & CODEX_cell_expression$Ki67 == 1 & CODEX_cell_expression$MHCI == 0 & CODEX_cell_expression$CD73==1] <-"Epithelial_prolif_CD73"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 0 & CODEX_cell_expression$CK == 1 & CODEX_cell_expression$Ki67 == 1 & CODEX_cell_expression$MHCI == 0 & CODEX_cell_expression$CD73==0] <-"Epithelial_prolif"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 0 & CODEX_cell_expression$CK == 1 & CODEX_cell_expression$Ki67 == 0 & CODEX_cell_expression$MHCI == 1 & CODEX_cell_expression$CD73==1] <-"Epithelial_nonprolif_MHCI_CD73"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 0 & CODEX_cell_expression$CK == 1 & CODEX_cell_expression$Ki67 == 0 & CODEX_cell_expression$MHCI == 1 & CODEX_cell_expression$CD73==0] <-"Epithelial_nonprolif_MHCI"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 0 & CODEX_cell_expression$CK == 1 & CODEX_cell_expression$Ki67 == 0 & CODEX_cell_expression$MHCI == 0 & CODEX_cell_expression$CD73==1] <-"Epithelial_nonprolif_CD73"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 0 & CODEX_cell_expression$CK == 1 & CODEX_cell_expression$Ki67 == 0 & CODEX_cell_expression$MHCI == 0 & CODEX_cell_expression$CD73==0] <-"Epithelial_nonprolif"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 0 & CODEX_cell_expression$CK == 0 & CODEX_cell_expression$CD31 == 1] <-"Endothelial"
CODEX_cell_expression$cell_type[CODEX_cell_expression$CD45 == 0 & CODEX_cell_expression$CK == 0 & CODEX_cell_expression$CD31 == 0] <-"Stromal_other"

### identify cellular neighborhood ###
cell_type_detail = unique(CODEX_cell_expression$cell_type)
my_neighbor_list <- list()
for (clstype in cell_type_detail){
  my.neighbor = c()
  print (clstype)
  
  for(sam in unique(CODEX_cell_expression$Sample_ID)) {
    print (sam)
    tmp.meta = CODEX_cell_expression[CODEX_cell_expression$cell_type == clstype & CODEX_cell_expression$Sample_ID == sam, ]
    tmp.meta1 = CODEX_cell_expression[CODEX_cell_expression$Sample_ID == sam, ]
    
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

### obtian and process neighborhood matrix ###
neighborhood_df <- c()
for(i in 1: length(my_neighbor_list)){
  name = names(my_neighbor_list)[i]
  df = my_neighbor_list[[i]]
  df <- as.data.frame(t(df))
  df$center_cell_type <- name
  neighborhood_df <- rbind(neighborhood_df, df)
}

neighborhood_df$cell_ID <- rownames(neighborhood_df)
neighborhood_df <- separate(neighborhood_df, col="cell_ID", into=c("C1","C2","C3"),sep="_")
colnames(neighborhood_df)[22] <- "patient_ID"
neighborhood_df$patient_ID <- paste0("Acc",neighborhood_df$patient_ID)
colnames(neighborhood_df)[23] <- "Timepoint"
neighborhood_df$Genetic <- "No"
neighborhood_df$Genetic[neighborhood_df$patient_ID %in% c("Acc31","Acc38","Acc124","Acc149","Acc157","Acc160","Acc162")] <- "PPP2R1A"
neighborhood_df$Genetic[neighborhood_df$patient_ID %in% c("Acc121")] <- "AKT"

neighborhood_df$PPP2R1A_mutation <- "No"
neighborhood_df$PPP2R1A_mutation[neighborhood_df$Genetic == "PPP2R1A"] <- "Yes"
neighborhood_df$mutation_timepoint <- paste0(neighborhood_df$PPP2R1A_mutation,"_", neighborhood_df$Timepoint)
neighborhood_df$mutation_timepoint <- factor(neighborhood_df$mutation_timepoint, levels = c("Yes_pre","Yes_on","No_pre","No_on"))

neighborhood_df_log <- neighborhood_df
neighborhood_df_log[,1:20] <- log10(neighborhood_df_log[,1:20] + 1)

ggplot(filter(neighborhood_df_log, neighborhood_df_log$center_cell_type %in% c("Epithelial_prolif","Epithelial_prolif_CD73","Epithelial_nonprolif","Epithelial_nonprolif_CD73",
                                                                               "Epithelial_prolif_MHCI","Epithelial_prolif_MHCI_CD73","Epithelial_nonprolif_MHCI","Epithelial_nonprolif_MHCI_CD73")), 
       aes(x=mutation_timepoint, y=T_CD8_CD45RO))+
  geom_boxplot(outlier.size = 0,width=0.5, aes(color=mutation_timepoint))+
  geom_point(aes(color=mutation_timepoint),size=1)+
  scale_color_manual(values=c("#fee391","#ec7014","#fee391","#ec7014"))+
  theme_classic()+
  theme(axis.ticks.length = unit(0.2,'cm'))+
  theme(legend.position = "none", plot.title = element_text(size=6), 
        axis.title = element_text(size=14), axis.text.x = element_text(angle = 45, hjust=1, size=10),
        axis.title.x = element_text(size=10), axis.title.y = element_text(size=10))+
  guides(fill=guide_legend(title = NULL))

### The abundance of other cell states in tumor neighborhood is plotted in a similar way.



### Fig. 3h-i were from CODEX images ###