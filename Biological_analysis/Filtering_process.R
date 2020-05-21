### Title: Filtering steps, metastatic cancer mutational signature analysis
### Author: Sander Vlugter
### Date: 17/May/2020
### Description:
### Using the HMF dataset on solid metastatic cancers single base substitutions of samples, splitting these 
### somatic mutations based on in which stage of the cancer they occured. During Primary tumour stage (clonal),
### during metastatic stage (sub-clonal). Here we filter out any samples of the HMF database that do not match our set standards.
### The Splitting of the somatic mutations into clonal and sub-clonal was done using PURPLE pipeline. This was done by other members before
### this study was conducted.
#######
library(devtools)
library(MutationalPatterns)
### load in the data
HMF_SNV_ALL<- readRDS("path/to/HMF_SNV_ALL.rds")
HMF_SNV_clonal<- readRDS("path/to/HMF_SNV_CLONAL.rds")
HMF_SNV_subclonal <-readRDS("path/to/HMF_SNV_SUBCLONAL.rds")

#### remove empty and or low mutation count samples from the subclonal
# threshold of atleast 50 SBS mutations in the subclonal stage
HMF_SNV_subclonal <- as.matrix(HMF_SNV_subclonal,drop=F)
subclonal <- HMF_SNV_subclonal[,which(colSums(HMF_SNV_subclonal)>=50)]
#### remove the same samples from both the combined and the clonal stage matrices
HMF_SNV_clonal <- as.matrix(HMF_SNV_clonal,drop=F)
clonal <- HMF_SNV_clonal[,match(colnames(subclonal),colnames(HMF_SNV_clonal))]
#### combined
HMF_SNV_ALL <- as.matrix(HMF_SNV_ALL,drop=F)
combined <- HMF_SNV_ALL[,match(colnames(subclonal),colnames(HMF_SNV_ALL))]
### one of the samples in the clonal stage seems to be empty, remove that sample from all of the matrices
clonal <- clonal[,-8] 
subclonal <- subclonal[,-8]
combined<- combined[,-8]
#### Filter based on the tumour purity
hmf_meta <- read.csv("~/path/to/metadata/upd_metadata.tsv", sep = '\t',header = T)
### threshold set to purity of atleast 0.2
hmf_meta_pure <- hmf_meta[which(hmf_meta$tumorPurity>=0.20),]
### select the samples that do not fit the purity standards
X<-colnames(subclonal)[match(hmf_meta_pure$sampleId,colnames(subclonal))]
Xfilter <- X[!is.na(X)]
### apply the filter to the matrices
subclonal<- subclonal[,match(Xfilter,colnames(subclonal))]
clonal <- clonal[,match(colnames(subclonal),colnames(clonal))]
combined <- combined[,match(colnames(subclonal),colnames(combined))]
############ save outputs for use in other scripts
saveRDS(subclonal, file = "~/path/to/hmfdata/subclonal.RDS")
saveRDS(clonal, file = "~/path/to/hmfdata/clonal.RDS")
saveRDS(combined, file = "~/path/to/hmfdata/combined.RDS")
saveRDS(hmf_meta_pure,file = "~/path/to/hmfdata/hmf_meta_pure.RDS")
