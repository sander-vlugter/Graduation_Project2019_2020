### Title: Mutational Signature Analysis, known and unknown aetiologies
## Author: Sander Vlugter
## Date (final version): 18/May/2020
## Description:
## Here we select certain samples that display a shift in contribution of several mutational signatures between the clonal
## and sub-clonal stage of the sample. The mutational signatures of intrest are based on wether they have a known aetiology
## or unknown aetiology. This is done in order to provide more information on the origin of the mutational signatures with
## unknown aetiology. Finding if the mutational signatures of unknown aetiology are most likely caused due to endogenous 
## processes or exogonous agents.
## using data generated in: "Low_cosine_clonal_subclonal_stages.R"
#####
library(MutationalPatterns)
###### Mutational Signatures of intrest:#####
## Mutational signatures with unknwon aetiology:
aetio_unknown <- c("SBS8","SBS12","SBS16","SBS17a","SBS17b","SBS23","SBS28","SBS33","SBS34","SBS37","SBS38","SBS39","SBS41")
# Additional info on these signatures:
# SBS12: contributes <20% of mutations in liver cancer samples
# SBS28: is found in most samples with SBS10 a/b
# SBS38: found ONLY in ultraviolet light associated melanomas suggesting potential indirect damage from UV-light
## Mutational signatures found in cancers - prior treatment with platinum drug
aetio_plati <- c("SBS31","SBS35")
# Additional info on these signatures:
# exhibit similar components of their signature
## Mutational signatures known to be caused by tabbaco smoke
aetio_smoke <- c("SBS4")
# Additional info:
# SBS 5 also has a link to smoking but spread over many cancer types, therefore measuring a decline or increase between -
# clonal and sub-clonal would make it to difficult as the subclonal location might be one of these cancer types.
# SBS 29 relation to chewing tabacco therefore mostly found in oral and not taken with this comparison
## Mutational signatures caused by UV light
aetio_uv <- c("SBS7a","SBS7b","SBS7c","SBS7d")
###### Load in the required data ######
hmf_meta_pure<- readRDS("~/path/to/hmfdata/hmf_meta_pure.RDS")
fitted_clonal<-readRDS("~/path/to/fitted_sig/clonal_fitted.RDS")
fitted_combined<-readRDS("~/path/to/fitted_sig/combined_fitted.RDS")
fitted_subclonal<-readRDS("~/path/to/fitted_sig/subclonal_fitted.RDS")
cosine_sig_profiles<-readRDS("~/path/to/fitted_sig/cosine_signature_profiles.RDS")

## load the COSMIC mutational signatures from MutationalPatterns
cosmic_signatures <- readRDS(system.file("states/COSMIC_signatures.rds",
                                         package = "MutationalPatterns"))
## remove signatures known to be caused by sequencing artefacts
cosmic_signatures$snv =cosmic_signatures$snv[,-c(32,48,50:65)]

######### Calculate percentile weight of the fitted data ###############
# Use the fitted data to calculate the percentile contribution to the mutational signature profile of the stages
## based on the fitted matrix
relativ_column <- function(df){
  df <- as.data.frame(df,drop=F)
  rel_df <- NULL
  for(x in colnames(df)){
    rel_df[[x]]<- as.data.frame(df[x]/sum(df[x]))
  }
  rel_df <- as.data.frame(rel_df, drop =F)
  return(rel_df)
}
## subset the fitted data contribution to their own matrices
# these are used in calculating the percentile contribution of the mutational signatures in the stages
fitted_clo_contri <- fitted_clonal$contribution
fitted_sub_contri <- fitted_subclonal$contribution
fitted_com_contri <- fitted_combined$contribution
# calculate the percentage of contribution for each mutational signature to the mutational signature of the profile,
# of each stage per sample
relativ_clo_contri <- relativ_column(fitted_clo_contri)
relativ_sub_contri <- relativ_column(fitted_sub_contri)
relativ_com_contri <- relativ_column(fitted_com_contri)
#### combine cosine similarity data with metadata entries
cosine_sig_profiles$sampleId <- rownames(cosine_sig_profiles)
## select only columns of importance
hmf_meta_subset<-hmf_meta_pure[match(colnames(relativ_com_contri),hmf_meta_pure$sampleId),c(3,4,5,8,18)]
library(tidyverse)
cosine_meta<- left_join(cosine_sig_profiles,hmf_meta_subset,by=c("sampleId"))
rownames(cosine_meta)<- cosine_meta$sampleId
############ Known aetiology Mutational signature caused by smoking ############################
#### select samples based on primary tumour location == "Lung" & Biopsy location being the same for the samples
# in the thesis we used lung primary lymph biopsy
lung_lymph_samples <- rownames(cosine_meta)[which(cosine_meta$primaryTumorLocation =="Lung" &
                                                    cosine_meta$biopsySite == "Lymph node")]
# clonal stage relative contribution of the smoke mutational signature
relativ_clo_contri[aetio_smoke,lung_lymph_samples]
# sub-clonal stage relative contribution of the smoke mutational signature
relativ_sub_contri[aetio_smoke,lung_lymph_samples]
# as splitting the mutations into clonal and sub-clonal isnt exact select samples that provide a good example
which(relativ_sub_contri[aetio_smoke,lung_lymph_samples]== 0)
which(relativ_clo_contri[aetio_smoke,lung_lymph_samples]== 0)
# select sample names
lung_lymph_names<- colnames(relativ_clo_contri[,c(lung_lymph_samples[c(5,12,15,17,18,21)])])
############ Known aetiology Mutational signatures caused by UV light ############################
skin_lymph_samples <- rownames(cosine_meta)[which(cosine_meta$primaryTumorLocation == "Skin"&
                                                    cosine_meta$biopsySite == "Lymph node")]
# clonal stage relative contribution of the UV mutational signature
relativ_clo_contri[aetio_uv,skin_lymph_samples]
# sub-clonal stage relative contribution of the UV mutational signature
relativ_sub_contri[aetio_uv,skin_lymph_samples]
# select sample names that show the point clearly
skin_lymph_names <- colnames(relativ_clo_contri[,c(skin_lymph_samples[c(7,12,14,17,18,23,26,28)])])
############ Known aetiology Mutational signatures caused by prior treatment with platinum based drugs ############################
platinum_samples <- cosine_meta %>%
  filter(str_detect(cosine_meta$treatment, "platin"))
platinum_names <- platinum_samples$sampleId

platt_rel_sub_fit <- relativ_sub_contri[,platinum_names]
plat_hi_diff <- colnames(platt_rel_sub_fit)[which(platt_rel_sub_fit["SBS31",platinum_names]>0.25|platt_rel_sub_fit["SBS35",platinum_names]>0.25)]

platinum_hi_meta<- hmf_meta_pure[match(plat_hi_diff,hmf_meta_pure$sampleId),]
ovary_plat_names <- as.vector(platinum_hi_meta$sampleId[which(platinum_hi_meta$primaryTumorLocation =="Ovary")])

############## mutational signature analyses of unknown aetiology ##################################
#### one by one relative contri in both clonal and subclonal for unknown aetio
sig_rela_comp <- function(sbs){
  for(sbs in sbs){
    print(sbs)
    print("Higher relative contribution in clonal-stage for this number samples:")
    print(length(colnames(relativ_clo_contri)[which(relativ_clo_contri[sbs,] > relativ_sub_contri[sbs,] &
                                                      (relativ_clo_contri[sbs,] -relativ_sub_contri[sbs,])> 0.05 )]))
    print("Higher relative contribution in sub-clonal stage for this number of samples:")
    print(length(colnames(relativ_clo_contri)[which(relativ_clo_contri[sbs,] < relativ_sub_contri[sbs,] &
                                                (relativ_clo_contri[sbs,] -relativ_sub_contri[sbs,])< -0.05 )]))
    print("where the difference between contribution is atleast 5%")}
}
sig_rela_comp(aetio_unknown)
### sbs37 and sbs 39 seem higher preff for subclonal
### sbs8 seems like a tie
### sbs 41 seems like higher in clonal

#### function that gives the primary and biopsy location ########
# of samples that show the relative signature contribution to differ by atleast 5%
sig_rela_locations <- function(sbs){
  for(sbs in sbs){
    hi_clo_name <- colnames(relativ_clo_contri)[which(relativ_clo_contri[sbs,] > relativ_sub_contri[sbs,] & (relativ_clo_contri[sbs,] -relativ_sub_contri[sbs,])> 0.05 )]
    hi_sub_name <- colnames(relativ_clo_contri)[which(relativ_clo_contri[sbs,] < relativ_sub_contri[sbs,] & (relativ_clo_contri[sbs,] -relativ_sub_contri[sbs,])< -0.05 )]
    
    print(sbs)
    print("Table of primary site, for higher relative signature in clonal")
    print(table(cosine_meta$primaryTumorLocation[match(hi_clo_name, cosine_meta$sampleId)]))
    print("Table of Biopsy site, for higher relative signatue in sub clonal")
    biop<- as.vector(cosine_meta$biopsySite[match(hi_sub_name, cosine_meta$sampleId)])
    print(length(biop))
    print(table(biop))
  }
  
}
sig_rela_locations("SBS41")
## high colon/rectum primary
rel_hi_sbs41_samples<-cosine_meta$sampleId[match(colnames(relativ_clo_contri)[which(relativ_clo_contri["SBS41",]>0.10)],cosine_meta$sampleId)]
cosine_meta$primaryTumorLocation[match(rel_hi_sbs41_samples, cosine_meta$sampleId)]
colon_sbs41_samples<- rel_hi_sbs41_samples[which(relativ_sub_contri["SBS41",rel_hi_sbs41_samples]==0)]


########################################## Create the combined plot ##############################################
library(MutationalPatterns)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
### get a set of colors
colorvector <- default_colors_ggplot(ncol(cosmic_signatures$snv))
colorvector[45]<- "#000033"
colorvector[c(4,7,8,9,10,35,39)] <- c("#666666","#FFFF33","#FFFF33","#FFFF33","#FFFF33","#CCFFFF","#CCFFFF")
colorvector
### set order of signatures we are intrested in
order_ofint <-c(4,7,8,9,10,35,39,45,1,2,3,5,6,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,36,37,38,40,41,42,43,44,46,47,48,49)
### alter the colors accordingly
colorvector <-colorvector[order_ofint]
colorvector[9:49] <- "#009900"
colorvector[6:7] <- "#3399FF"
colorvector[8]<- "#0000CC"
colorvector[1]<- "#CC0000"
colorvector
###### Add pseudo numbers to the fitted contribution plots so the legends stay matched
fit_aetio_clo_all <- fitted_clonal$contribution + 0.0000001
fit_aetio_sub_all <- fitted_subclonal$contribution +0.0000001

##### UV - Skin to lymph
UV_clonal_skin_lymph <- plot_contribution(fit_aetio_clo_all[order_ofint,skin_lymph_names],cosmic_signatures$snv[,order_ofint], type = "snv",mode = "relative", palette = colorvector)
UV_subclo_skin_lymph <- plot_contribution(fit_aetio_sub_all[order_ofint,skin_lymph_names],cosmic_signatures$snv[,order_ofint],type = "snv",mode = "relative", palette = colorvector)
#### smoke - lung - lymph
Smoke_clonal_lung_lymph<- plot_contribution(fit_aetio_clo_all[order_ofint,lung_lymph_names],cosmic_signatures$snv[,order_ofint], type = "snv",mode = "relative", palette = colorvector)
Smoke_subclo_lung_lymph<- plot_contribution(fit_aetio_sub_all[order_ofint,lung_lymph_names],cosmic_signatures$snv[,order_ofint],type = "snv",mode = "relative", palette = colorvector)
### platinum drugs - ovary
plat_clonal_ovary<- plot_contribution(fit_aetio_clo_all[order_ofint,ovary_plat_names],cosmic_signatures$snv[,order_ofint], type = "snv",mode = "relative", palette = colorvector)
plat_subclo_ovary<- plot_contribution(fit_aetio_sub_all[order_ofint,ovary_plat_names],cosmic_signatures$snv[,order_ofint],type = "snv",mode = "relative", palette = colorvector)
### unknown aetio sig41  rect/colon 
sig41_clonal_colon_liver<- plot_contribution(fit_aetio_clo_all[order_ofint,colon_sbs41_samples],cosmic_signatures$snv[,order_ofint],type = "snv",method = "relative", palette = colorvector)
sig41_subclo_colon_liver<- plot_contribution(fit_aetio_sub_all[order_ofint,colon_sbs41_samples],cosmic_signatures$snv[,order_ofint],type = "snv",method = "relative", palette = colorvector)
legend <- sig41_clonal_colon_liver+theme(legend.direction = "vertical", legend.title = element_text(size = 12,face = "bold"), legend.text = element_text(size = 10))+
  guides(fill=guide_legend(ncol =2))

UV_leftside<- plot_grid(UV_clonal_skin_lymph+theme(legend.position = 'none', axis.text.x.bottom = element_blank(),axis.ticks = element_blank(),
                                                   axis.title.y = element_blank(),
                                                   panel.border = element_blank(), plot.margin = unit(c(0.2,0.2,-0.8,0.1),"cm")),
                        UV_subclo_skin_lymph+theme(legend.position = 'none', axis.text.x.bottom = element_blank(),axis.ticks = element_blank(),
                                                   axis.title.y = element_blank(),
                                                   panel.border = element_blank(), plot.margin = unit(c(0.1,0.2,-0.5,0.1),"cm")),nrow = 2)+
  theme(panel.border = element_rect(fill=NA, colour = "black",size = 1))


Smoke_middle<- plot_grid(Smoke_clonal_lung_lymph+theme(legend.position = 'none', axis.text.x.bottom = element_blank(),axis.ticks = element_blank(),
                                                       axis.text.y = element_blank(),axis.title.y = element_blank(),
                                                       panel.border = element_blank(), plot.margin = unit(c(0.2,0.2,-0.8,0.1),"cm")),
                         Smoke_subclo_lung_lymph+theme(legend.position = 'none', axis.text.x.bottom = element_blank(),axis.ticks = element_blank(),
                                                       axis.text.y = element_blank(),axis.title.y = element_blank(),
                                                       panel.border = element_blank(), plot.margin = unit(c(0.1,0.2,-0.5,0.1),"cm")),nrow=2)+
  theme(panel.border = element_rect(fill=NA, colour = "black",size = 1))

Plat_midright<-plot_grid(plat_clonal_ovary+theme(legend.position = 'none', axis.text.x.bottom = element_blank(),axis.ticks = element_blank(),
                                                 axis.text.y = element_blank(),axis.title.y = element_blank(),
                                                 panel.border = element_blank(), plot.margin = unit(c(0.2,0.2,-0.8,0.1),"cm")),
                         plat_subclo_ovary+theme(legend.position = 'none', axis.text.x.bottom = element_blank(),axis.ticks = element_blank(),
                                                 axis.text.y = element_blank(),axis.title.y = element_blank(),
                                                 panel.border = element_blank(), plot.margin = unit(c(0.1,0.2,-0.5,0.1),"cm")),nrow = 2)+
  theme(panel.border = element_rect(fill=NA, colour = "black",size = 1))

sig41_right<- plot_grid(sig41_clonal_colon_liver+theme(legend.position = 'none', axis.text.x.bottom = element_blank(),axis.ticks = element_blank(),
                                                       axis.text.y = element_blank(),axis.title.y = element_blank(),
                                                       panel.border = element_blank(), plot.margin = unit(c(0.2,0.2,-0.8,0.1),"cm")),
                        sig41_subclo_colon_liver+theme(legend.position = 'none', axis.text.x.bottom = element_blank(),axis.ticks = element_blank(),
                                                       axis.text.y = element_blank(),axis.title.y = element_blank(),
                                                       panel.border = element_blank(), plot.margin = unit(c(0.1,0.2,-0.5,0.1),"cm")),nrow = 2)+
  theme(panel.border = element_rect(fill=NA, colour = "black",size = 1))
### create legend
ofint_legend<- plot_contribution(fit_aetio_clo_all[order_ofint,skin_lymph_names],cosmic_signatures$snv[,order_ofint], type = "snv",mode = "relative", palette = colorvector)
## add the order of importance to the legend
legsig_important <- c("SBS4","SBS7a","SBS7b","SBS7c","SBS7d","SBS31","SBS35","SBS41","SBS1")
## create the plot on which the legend finds its origin
getlegend_alt_2<- ofint_legend+
  scale_fill_manual(breaks = legsig_important,values = colorvector, labels= c("SBS4","SBS7a","SBS7b","SBS7c","SBS7d","SBS31","SBS35","SBS41","Other"))
plot_grid(getlegend_alt_2)
#### grab the legend
overall_legend<-get_legend(getlegend_alt_2)
### add the legend
combo_plot_2<-plot_grid(UV_leftside,Smoke_middle,Plat_midright,sig41_right,overall_legend,nrow = 1, rel_widths = c(2,1,1,2,1),labels = c("A","B","C","D"),
                        hjust = c(-3,-1,-1,-1),vjust = 1.02)
####### add y axis text
y.grob <- textGrob("Relative contribution",
                   gp=gpar(col="black", fontsize=12), rot=90)
################ the plot ########
grid.arrange(arrangeGrob(combo_plot_2,left = y.grob))
### save the plot
pdf("~/path/to/images/unknown_aetio_sbs41.pdf",paper = "a4")
grid.arrange(arrangeGrob(combo_plot_2,left = y.grob))
dev.off()
