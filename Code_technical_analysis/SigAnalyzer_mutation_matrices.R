#### PCAWG data mutation matrices as by SigAnalyzer #####
### Author: Sander Vlugter
## Date: 15/May/2020
## Description: obtain the mutation matrices of the PCAWG samples as called by Siganalyzer

### data of SigAnalyzer retrieved from: 
# https://github.com/broadinstitute/getzlab-SignatureAnalyzer/tree/master/PCAWG/INPUT_SignatureAnalyzer
# contains PCAWG samples in mutation matrices as retrieved by SigAnalyzer
# download and load in these files:

load("~/path/to/downloaded_file/lego96.PAN.SNV.091217.RData")
load("~/path/to/downloaded_file/lego1536.INDEL.091217.RData")
load("~/path/to/downloaded_file/lego1536.DNP.091217.RData")

library(doParallel)
library(MutationalPatterns)
cosmic_signatures = readRDS(system.file("states/COSMIC_signatures.rds", package = "MutationalPatterns"))


#### create the mutation matrices as called by SigAnalyzer
rownames(lego96.SNV) = foreach(x = rownames(lego96.SNV), .combine = "c") %do% {
  return(sprintf("%s[%s>%s]%s", substr(x,3,3),substr(x,1,1),substr(x,2,2),substr(x,4,4)))
}
## set the roworder as that of mutationalpatterns
new_order = match(rownames(cosmic_signatures$snv), rownames(lego96.SNV))
lego96.SNV = lego96.SNV[new_order,]

rownames(lego1536.DNP) = sprintf("%s>%s", substr(rownames(lego1536.DNP), 1, 2), substr(rownames(lego1536.DNP), 4, 5))

rownames(lego1536.INDEL) = rownames(cosmic_signatures$indel)

sigana_mutationmatrix = list("snv" = lego96.SNV, "dbs" = lego1536.DNP, "indel"=lego1536.INDEL)
####
saveRDS(sigana_mutationmatrix, file = "~/path/to/place/save/sigana_mutation_matrices.RData")