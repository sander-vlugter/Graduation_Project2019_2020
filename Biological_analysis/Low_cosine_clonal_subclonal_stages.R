### Title: Analyses Low cosine similarity between clonal and sub-clonal stage mutational signature profiles
## Author: Sander Vlugter
## Date (final version): 17/May/2020
## Description:
## fitting the mutation matrices of the sample's stages (clonal, sub-clonal and combined) to cosmic signatures 
## next we determine the cosine similarity between the stages and the stages to the combined mutational signature profiles
## followed by selecting samples that have a low cosine similarity between the clonal and subclonal stages (lower than 0.2)
## for these samples we see if they have had any prior treatment of the primary cancer tumour.
## here we used the mutational signatures as defined by COSMIC
############
library(MutationalPatterns)

## load in the filtered data created with "Filtering_process.R"
subclonal<- readRDS("~/path/to/hmfdata/subclonal.RDS")
clonal <- readRDS("~/path/to/hmfdata/clonal.RDS")
combined <- readRDS("~/path/to/hmfdata/combined.RDS")
hmf_meta_pure<- readRDS("~/path/to/hmfdata/hmf_meta_pure.RDS")
## load the COSMIC mutational signatures from MutationalPatterns
cosmic_signatures <- readRDS(system.file("states/COSMIC_signatures.rds",
                                         package = "MutationalPatterns"))
## remove signatures known to be caused by sequencing artefacts
cosmic_signatures$snv =cosmic_signatures$snv[,-c(32,48,50:65)]
## fit the data to COSMIC signatures
fitted_combined <- fit_to_signatures(combined, type = 'snv', signatures = cosmic_signatures)
fitted_subclonal <- fit_to_signatures(subclonal, type = 'snv',signatures = cosmic_signatures)
fitted_clonal <- fit_to_signatures(clonal, type = 'snv',signatures = cosmic_signatures)
## save the fitted contribution objects for later use
saveRDS(fiited_clonal, "~/path/to/fitted_sig/clonal_fitted.RDS")
saveRDS(fitted_combined, "~/path/to/fitted_sig/combined_fitted.RDS")
saveRDS(fitted_subclonal, "~/path/to/fitted_sig/subclonal_fitted.RDS")

## calculate the cosine similarity of the fitted stages 
cosine_sig_clonal_subclonal <- diag(cos_sim_matrix(fitted_clonal$contribution, fitted_subclonal$contribution))
cosine_sig_combined_clonal <- diag(cos_sim_matrix(fitted_combined$contribution, fitted_clonal$contribution))
cosine_sig_combined_subclonal <- diag(cos_sim_matrix(fitted_combined$contribution, fitted_subclonal$contribution))
## put them together in a dataframe
cosine_sig_profiles<- cbind(cosine_sig_combined_clonal, cosine_sig_combined_subclonal,cosine_sig_clonal_subclonal)
cosine_sig_profiles <- as.data.frame(cosine_sig_profiles)
## save this object for later use
saveRDS(cosine_sig_profiles, "~/path/to/fitted_sig/cosine_signature_profiles.RDS")

## select samples that have a low cosine similarity between the mutational signature profiles of the clonal and sub-clonal stage
Id_Low_cosine <- rownames(cosine_sig_profiles)[which(cosine_sig_profiles$cosine_sig_clonal_subclonal <=0.2)]
## create the mutational signature profile plots for these samples
pl_subclonal_ex_diff <- plot_contribution(fitted_subclonal$contribution[,Id_Low_cosine],
                                          cosmic_signatures, type = "snv",mode = "relative")
pl_clonal_ex_diff <- plot_contribution(fitted_clonal$contribution[,Id_Low_cosine],
                                       cosmic_signatures, type = "snv",mode = "relative")
### match these samples with low cosine to the meta data, to find out about the pretreatment
df_2<-melt(hmf_meta_pure[match(Id_Low_cosine, hmf_meta_pure$sampleId),c("PreTreatment","sampleId")], id.vars = "sampleId")
## select those who had treatment
df_2$chemo <- grepl("Yes",df_2$value)
## reassigning factors
df_2$sampleId <- as.character(df_2$sampleId)
df_2$sampleId <- as.factor(df_2$sampleId)
## create treatment logical plot with matched sample order
treatment_pl<- ggplot(df_2,aes(sampleId, y=variable))+
  scale_fill_manual(values = c("white","darkgreen"))+
  geom_tile(aes(fill=chemo), color = "black")+
  theme(legend.position = 'none')+
  aes(x=fct_inorder(sampleId))
#### check to make sure the order is correct 
plot(treatment_pl)+theme(axis.text.x = element_text(angle = 90))
plot(pl_clonal_ex_diff)+theme(axis.text.x = element_text(angle = 90))
plot(pl_subclonal_ex_diff)+theme(axis.text.x = element_text(angle = 90))
## Combine the plots together Without any legend (a combined legend will be made)
library(gridExtra)
library(cowplot)
lowcosine_72_plot<- plot_grid(pl_clonal_ex_diff+ theme(legend.position = 'none', axis.text.x.bottom = element_blank(),axis.ticks = element_blank(),
                                                panel.border = element_blank(), plot.margin = unit(c(0.1,0.2,-0.5,0.1),"cm")),
                       pl_subclonal_ex_diff+ theme(legend.position = 'none',axis.text.x = element_blank(),axis.ticks = element_blank(),
                                                   panel.border = element_blank(),plot.margin = unit(c(0.1,0.2,-0.5,0.1),"cm")),
                       treatment_pl+theme(legend.position = 'none',axis.text.x = element_blank(),
                                              axis.text.y.left = element_blank(),
                                              axis.title = element_blank(),panel.background = element_blank(),axis.ticks = element_blank()
                                              ,plot.margin = unit(c(0.1,0.2,0.4,0.1),"cm")),
                       nrow = 3,
                       align = 'v', labels = c("A","B","C"), rel_heights = c(2,2,1))
## create a legend that is applicable to both 
# by using their combined stages legend
pl_combined_ex_diff <-plot_contribution(fitted_combined$contribution[,Id_Low_cosine],cosmic_signatures, type = "snv", mode= 'relative')+
  theme(legend.direction = "vertical", legend.title = element_text(size = 10), legend.text = element_text(size = 8))+
  guides(fill=guide_legend(ncol =2))
## extract the legend of this plot
pl_legend<- get_legend(pl_combined_ex_diff)
## inspect the combined plot with the legend added
plot_grid(lowcosine_72_plot,pl_legend, ncol = 2, rel_widths = c(2,.5),rel_heights = c(1,3))
## save for later purpose
# alter the rel_width/heights of the plot for your own needs, as well as the PDF options
pdf("~/path/to/images/low_cosine_betweenstages.pdf", paper = "a4")
plot_grid(lowcosine_72_plot,pl_legend, ncol = 2, rel_widths = c(2,.5),rel_heights = c(1,3))
dev.off()