### Fig. 1 ###

### load required packages ###
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(swimplot)
library(survival)
library(survminer)


### Fig. 1b ###
### read in clinical data ###
OCCC_response_endpoint <- read.csv("./Response_outcome_endpoint.csv")
OCCC_response_endpoint$Genetic <- "No"
OCCC_response_endpoint$Genetic[OCCC_response_endpoint$PPP2R1A_mutation==1] <- "PPP2R1Amut"
OCCC_response_endpoint$Genetic[OCCC_response_endpoint$AKT_alteration==1] <- "AKTalt"

OCCC_response_timeline <- read.csv("./Response_outcome.csv")

# Swimmer plot #
arm_plot <- swimmer_plot(df=OCCC_response_endpoint,id='Acc',end='Time_mo',name_fill='Genetic',
                         id_order='PPP2R1A_mutation', increasing=F, col="black",alpha=1,width=.8)
AE_plot <- arm_plot + 
           swimmer_points(df_points= OCCC_response_timeline,id='Acc',time='Time_mo',
                          name_shape = 'overallResponse',size=2,fill='white',col='black')+
           scale_y_continuous(breaks = c(0,6,12,18,24,30,36,42,48,54,60,66,72))+
           scale_fill_manual(name="Genetic",values=c("AKTalt"="#8fd2c4", "No"="#9bc1dd", "PPP2R1Amut"='#fdc583'))+
           scale_shape_manual(name="overallResponse",values=c(PD=17,SD=16,PR=15,CR=23,Death=3,Ongoing=9),breaks=c('PD','SD','PR','CR',"Death","Ongoing"))
AE_plot


### Fig. 1c ###
patient_metadata <- read.csv("./patient_metadata.csv")

clear_cell_data_OS <- Surv(time = patient_metadata$OS_mo, 
                           event = patient_metadata$OS_stauts)

clear_cell_data_OS_fit <- survfit(clear_cell_data_OS ~ PPP2R1A_mutation, 
                                  data = patient_metadata)

ggsurvplot(clear_cell_data_OS_fit, 
           data = patient_metadata,
           risk.table = TRUE,
           legend.labs = levels(patient_metadata$PPP2R1A_mutation), 
           palette=c("#83b2d3","#fbb463"),
           break.x.by = 12,
           ggtheme = theme_classic())


### Fig. 1d ###
clear_cell_data_PPP2R1A_wt <- filter(patient_metadata, patient_metadata$PPP2R1A_mutation == 0)

clear_cell_data_PPP2R1A_wt_OS <- Surv(time = clear_cell_data_PPP2R1A_wt$OS_mo, 
                                      event = clear_cell_data_PPP2R1A_wt$OS_stauts)

clear_cell_data_PPP2R1A_wt_OS_fit <- survfit(clear_cell_data_PPP2R1A_wt_OS ~ ARID1A_mutation, 
                                             data = clear_cell_data_PPP2R1A_wt)

ggsurvplot(clear_cell_data_PPP2R1A_wt_OS_fit, 
           data = clear_cell_data_PPP2R1A_wt,
           risk.table = TRUE,
           legend.labs = levels(clear_cell_data_PPP2R1A_wt$ARID1A_mutation), 
           palette=c("#83b2d3","#fbb463"),
           break.x.by = 12,
           ggtheme = theme_classic())


### Fig. 1e ###
clear_cell_data_ARID1A_mut <- filter(patient_metadata, patient_metadata$ARID1A_mutation == 1)

clear_cell_data_ARID1A_mut_OS <- Surv(time = clear_cell_data_ARID1A_mut$OS_mo, 
                                      event = clear_cell_data_ARID1A_mut$OS_stauts)

clear_cell_data_ARID1A_mut_OS_fit <- survfit(clear_cell_data_ARID1A_mut_OS ~ PPP2R1A_mutation, 
                                             data = clear_cell_data_ARID1A_mut)

ggsurvplot(clear_cell_data_ARID1A_mut_OS_fit, 
           data = clear_cell_data_ARID1A_mut,
           risk.table = TRUE,
           legend.labs = levels(clear_cell_data_ARID1A_mut$PPP2R1A_mutation), 
           palette=c("#83b2d3","#fbb463"),
           break.x.by = 12,
           ggtheme = theme_classic())


### Fig. 1g ###
tumor_measurement <- read.csv("tumor_measurement.csv")
OCCC_response_endpoint_meta <- OCCC_response_endpoint[,c(1,4:6)]
tumor_measurement <- left_join(tumor_measurement, OCCC_response_endpoint_meta, by="Acc")

### spider plot ###
ggplot(filter(tumor_measurement, tumor_measurement$Genetic %in% c("AKTalt","No")), 
       aes(x=Time_mo, y=pct_change_from_BL, color=Genetic))+
  geom_point(size=0,alpha=0)+
  geom_line(aes(group=Acc))+
  geom_hline(yintercept = -30)+
  scale_color_manual(values=c("AKTalt"="#8fd2c4", "No"="#9bc1dd","PPP2R1Amut"='#fdc583'))+
  scale_x_continuous(limits=c(0,38), breaks = c(0,6,12,18,24,30,36))+
  scale_y_continuous(limits=c(-100,125), breaks = c(-100, -50, 0, 50, 100))+
  theme_classic()

ggplot(filter(tumor_measurement, tumor_measurement$Genetic %in% c("PPP2R1Amut")), 
       aes(x=Time_mo, y=pct_change_from_BL, color=Genetic))+
  geom_point(size=0,alpha=0)+
  geom_line(aes(group=Acc))+
  geom_hline(yintercept = -30)+
  scale_color_manual(values=c("AKTalt"="#8fd2c4", "No"="#9bc1dd","PPP2R1Amut"='#fdc583'))+
  scale_x_continuous(limits=c(0,38), breaks = c(0,6,12,18,24,30,36))+
  scale_y_continuous(limits=c(-100,125), breaks = c(-100, -50, 0, 50, 100))+
  theme_classic()


### Fig. 1h ###
tumor_measurement_for_waterfall_plot <- tumor_measurement[, c("Acc","Genetic","Best_response")]
tumor_measurement_for_waterfall_plot <- unique(tumor_measurement_for_waterfall_plot)
tumor_measurement_for_waterfall_plot <- tumor_measurement_for_waterfall_plot[order(tumor_measurement_for_waterfall_plot$Best_response, decreasing = T),]
# since Pt 152 and Pt 157 both showed best response of -100% change from baseline, here for visualization , we flip the order of the two patients
tumor_measurement_for_waterfall_plot <- tumor_measurement_for_waterfall_plot[c(1:29, 31, 30),]
tumor_measurement_for_waterfall_plot$plot_id <- seq(1:31)

ggplot(tumor_measurement_for_waterfall_plot, aes(x=plot_id, y=Best_response, fill=Genetic))+
  geom_bar(stat = "identity")+
  theme_classic()+
  geom_hline(yintercept = -30)+
  scale_fill_manual(values=c("AKTalt"="#8fd2c4", "No"="#9bc1dd","PPP2R1Amut"='#fdc583'))+
  theme_classic()