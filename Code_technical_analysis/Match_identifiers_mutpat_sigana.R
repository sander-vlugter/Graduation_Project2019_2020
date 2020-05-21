### Matching the sampleID's of PCAWG samples of MutationalPatterns with Siganalyzer on basis of metadata ####
## Author: Sander Vlugter
## Date: 16/May/2020
## Description: as both tools have different identifiers we have to match the samples based on the metadata they have in common

##### load the mutation matrices as called by MutationalPatterns, created in "MutationalPatterns_vcf_to_mutation_matrix".
library(MutationalPatterns)
snv_mutpat<- readRDS("~/path/to/save_location/sbs_mutpat.RData")
dbs_mutpat<- readRDS("~/path/to/save_location/dbs_mutpat.RData")
indel_mutpat<- readRDS("~/path/to/save_location/indel_mutpat.RData")

#### load the mutation matrices as called by SigAnalyzer, created in "SigAnalyzer_mutation_matrices".
sigana_mutationmatrix<- readRDS("~/path/to/place/save/sigana_mutation_matrices.RData")

### load in the meta data files on the PCAWG data
# consists of 2 files, OverView_metadata & meta_data
library(readxl)
###### Loading in the meta data OVERVIEW (DO_id & numeral)####
overview_metadata <- read_excel("~/path/to/overview_metadata")
### keep only the needed information on the Do_id & numeral_id
overview_metadata <- cbind(overview_metadata$SAMPLE,overview_metadata$donor_ID_manifest_all)
colnames(overview_metadata) <- c("Numeral_do","DO_id")
## remove empty entries
overview_metadata <- na.omit(overview_metadata)
## make it a dataframe
overview_metadata<- as.data.frame(overview_metadata)

### loading in the meta data (DO_id & SP_id)
spid_doid<- read_excel("~/path/to/meta_data")
spid_doid <- cbind(spid_doid$...1,spid_doid$...9, spid_doid$...28)
colnames(spid_doid) <- c("SP_id","DO_id","library_type")
#remove first row containing colnames
# spid_doid<- spid_doid[-1,] 
# make it a dataframe
spid_doid <- as.data.frame(spid_doid)
# exclude the library type of RNA-seq
spid_doid <- spid_doid[spid_doid$library_type != "RNA-Seq",]
# find the possible dupes
spid_doid <- spid_doid[!duplicated(spid_doid),]
spid_doid$SP_id[which(duplicated(spid_doid$DO_id))]
spid_doid$DO_id[which(duplicated(spid_doid$DO_id))]

### As the samples as called by SigAnalyzer have a different identifier code than MutationalPatterns
# they are matched on the matches of these identifiers as found in the metadata files
# sigAnalyzer uses SP_id, whilst MutPat uses numeral_id, in both metadata files they have a matching other id: DO_id

# select the SP_id from siganalyzer 
SP_name_sigana <- colnames(sigana_mutationmatrix$snv)
# remove the additions
SP_name_sigana <- do.call(cbind,strsplit(SP_name_sigana,"__"))[2,]
SP_name_sigana <- as.data.frame(SP_name_sigana)
# keep in mind the dupes
SP_name_sigana$SP_name_sigana[spid_doid$SP_id[which(duplicated(spid_doid$DO_id))]]

#### now we have the "loose" SP-ids used by sigana, and a list of SP_id with the matching DO_id and
# a list of numeral_id with matchin do_id
#### match the DO_id "numeral version" with the SP_id"numeral version"

match(spid_doid$DO_id, overview_metadata$DO_id)
which(is.na(match(spid_doid$DO_id, overview_metadata$DO_id)))
### create a list that contains a full matching set between Numeral_id DO_id SP_id
SP_num_do_meta<- cbind(overview_metadata$Numeral_do[match(spid_doid$DO_id, overview_metadata$DO_id)], spid_doid)
## remove empty's
SP_num_do_meta <- na.omit(SP_num_do_meta)
### create a dataframe with the matching numeral_id colnames of mutpat snv matrix
match(SP_num_do_meta$`overview_metadata$Numeral_do[match(spid_doid$DO_id, overview_metadata$DO_id)]`, colnames(snv_mutpat))
snv_mutpat_sp_num_do_meta <- cbind(colnames(snv_mutpat)[match(SP_num_do_meta$`overview_metadata$Numeral_do[match(spid_doid$DO_id, overview_metadata$DO_id)]`, colnames(snv_mutpat))], SP_num_do_meta)

######### match the SP_id in the list to the SP_name of sigana 

SP_SPname <- cbind(SP_name_sigana, colnames(sigana_mutationmatrix$snv))
### create a list where everything is matched
SP_num_SPname_do_meta <- cbind(SP_SPname$`colnames(sigana_mutationmatrix$snv)`[match(SP_num_do_meta$SP_id,SP_SPname$SP_name_sigana)], SP_num_do_meta)

SP_num_SPname_do_meta2 <- SP_num_SPname_do_meta
colnames(SP_num_SPname_do_meta2) <- c("name_sp","numeral","SP_id","Do_id","library_type")

#### Rename the matrices of MutationalPatterns for easier comparison
snv_mutpat_rename <- snv_mutpat
colnames(snv_mutpat_rename)<-  SP_num_SPname_do_meta2$name_sp[match(colnames(snv_mutpat_rename),SP_num_SPname_do_meta2$numeral)]
snv_mutpat_rename<- as.data.frame(snv_mutpat_rename)

##### remove samples of which there is no matching Identification code
snv_mutpat_rename2 <- snv_mutpat_rename[!is.na(colnames(snv_mutpat_rename))]
snv_mutpat_rename2 <- as.matrix(snv_mutpat_rename2)
match(colnames(snv_mutpat_rename2),colnames(sigana_mutationmatrix$snv))
#### make the list of matrices from sigana into seperate matrix of SBS, with only matching samples
snv_sigana <- sigana_mutationmatrix$snv[,match(colnames(snv_mutpat_rename2),colnames(sigana_mutationmatrix$snv))]

########## finding all the duplicate DO_id's with the SP_id's ###########
all_do_multi<- spid_doid$DO_id[which(duplicated(spid_doid$DO_id))]           # 121
uni_do_multi<- unique(spid_doid$DO_id[which(duplicated(spid_doid$DO_id))])   # 63 unique (meaning only 2 or more samples )

multi_do_subset<- SP_num_SPname_do_meta2$name_sp[match(all_do_multi, SP_num_SPname_do_meta2$Do_id)]

dupli_do_ids <-SP_num_SPname_do_meta2[SP_num_SPname_do_meta2$Do_id %in% all_do_multi,] ## 182
match(dupli_do_ids$name_sp,colnames(sigana_mutationmatrix$snv)) # 182
match(dupli_do_ids$numeral,colnames(snv_mutpat))
snv_mutpat_mat <- as.matrix(snv_mutpat)
### select the ones that are correct based on cos_sim
dupli_snv_cos<- cos_sim_matrix(sigana_mutationmatrix$snv[,match(dupli_do_ids$name_sp,colnames(sigana_mutationmatrix$snv))],snv_mutpat_mat[,unique(match(dupli_do_ids$numeral,colnames(snv_mutpat)))])
library(tidyverse)
list_do_dupli <-dupli_do_ids %>%
  group_split(Do_id)

for(d in list_do_dupli){
  cos_dub_snv <- cos_sim_matrix(sigana_mutationmatrix$snv[,match(d$name_sp,colnames(sigana_mutationmatrix$snv))],snv_mutpat_mat[,(match(d$numeral,colnames(snv_mutpat_mat)))])
  cos_dub_snv <- as.data.frame(cos_dub_snv)
  print(cos_dub_snv)
}
snv_mutpat_rename3 <- snv_mutpat_rename2
############# replace the SP_names that were incorrectly matched###########

for(d in list_do_dupli){
  cos_dub_snv <- cos_sim_matrix(sigana_mutationmatrix$snv[,match(d$name_sp,colnames(sigana_mutationmatrix$snv))],snv_mutpat_mat[,(match(d$numeral,colnames(snv_mutpat_mat)))])
  cos_dub_snv <- as.data.frame(cos_dub_snv)
  #  print(cos_dub_snv[1])
  # print(rownames(cos_dub_snv)[which(cos_dub_snv[,1] >=0.9999999)])
  print(cos_dub_snv[cos_dub_snv[,1] >=0.9999999])
  print(rownames(cos_dub_snv)[which(cos_dub_snv[,1] >=0.9999999)])
  # print(colnames(cos_dub_snv)[1])
  colnames(snv_mutpat_rename3)[match(rownames(cos_dub_snv)[1],colnames(snv_mutpat_rename3))] <- rownames(cos_dub_snv)[which(cos_dub_snv[,1] >=0.9999999)]
}
###################################### create the correctly matched matrices ##########
# for both siganalyzer and mutationalpatterns
snv_sigana_2 <- sigana_mutationmatrix$snv[,match(colnames(snv_mutpat_rename3),colnames(sigana_mutationmatrix$snv))]


### make a list of the samples that have duplicated indentifiers and match the correct one
dupe_do <- NULL
for(d in list_do_dupli){
  cos_dub_snv <- cos_sim_matrix(sigana_mutationmatrix$snv[,match(d$name_sp,colnames(sigana_mutationmatrix$snv))],snv_mutpat_mat[,(match(d$numeral,colnames(snv_mutpat_mat)))])
  cos_dub_snv <- as.data.frame(cos_dub_snv)
  #    print(cos_dub_snv[1])
  #   print(rownames(cos_dub_snv)[which(cos_dub_snv[,1] >=0.9999999)])
  print(cos_dub_snv[cos_dub_snv[,1] >=0.9999999])
  print(rownames(cos_dub_snv)[which(cos_dub_snv[,1] >=0.9999999)])
  print(colnames(cos_dub_snv)[1])
  temp <- data.frame(rownames(cos_dub_snv)[which(cos_dub_snv[,1] >=0.9999999)],colnames(cos_dub_snv)[1])   
  # do.call(rbind,k) #### fix dit
  dupe_do <- rbind(dupe_do,temp)
}
colnames(dupe_do)<- c("SP_names","numeral")
############### Save the SBS matched mutation matrices of the PCAWG samples as called by MutationalPAtterns & SigAnalyzer

saveRDS(snv_mutpat_rename3, file = "~/path/to/saving/matched/sbs_matrix/mutationalpatterns_snv_countmatrices.RData")
saveRDS(snv_sigana_2, file = "~/path/to/saving/matched/sbs_matrix/siganalyzer_snv_countmatrices.Rdata")



#################################### DBS matrices matching ##############################
dbs_mutpat_rename <- dbs_mutpat
colnames(dbs_mutpat_rename)<- SP_num_SPname_do_meta2$name_sp[match(colnames(dbs_mutpat_rename),SP_num_SPname_do_meta2$numeral)]
dbs_mutpat_rename <- dbs_mutpat_rename[!is.na(colnames(dbs_mutpat_rename))]
dbs_mutpat_rename <- as.matrix(dbs_mutpat_rename)

#### set the matching identifiers based on cosine of snv
dbs_mutpat_rename2 <- dbs_mutpat_rename
for(d in list_do_dupli){
  cos_dub_snv <- cos_sim_matrix(sigana_mutationmatrix$snv[,match(d$name_sp,colnames(sigana_mutationmatrix$snv))],snv_mutpat_mat[,(match(d$numeral,colnames(snv_mutpat_mat)))])
  cos_dub_snv <- as.data.frame(cos_dub_snv)
  #  print(cos_dub_snv[1])
  # print(rownames(cos_dub_snv)[which(cos_dub_snv[,1] >=0.9999999)])
  print(cos_dub_snv[cos_dub_snv[,1] >=0.9999999])
  print(rownames(cos_dub_snv)[which(cos_dub_snv[,1] >=0.9999999)])
  # print(colnames(cos_dub_snv)[1])
  colnames(dbs_mutpat_rename2)[match(rownames(cos_dub_snv)[1],colnames(snv_mutpat_rename3))] <- rownames(cos_dub_snv)[which(cos_dub_snv[,1] >=0.99999999)]
}

### get the matchin matrices of siganalyzer 
sigana_dbs <- sigana_mutationmatrix$dbs[match(colnames(dbs_mutpat_rename2),colnames(sigana_mutationmatrix$dbs))]
#### SigAnalyzer uses a different annotation method than MutationalPatterns, even though the context is the same
# SigAnalyzer doesn't use the correct strand-agnostic annotations (displays the same DBS type but annotation isnt as COSMIC)
match(rownames(sigana_dbs),rownames(dbs_mutpat_rename2))
## fix the contexts dbs

row_sig_ana_dbs<- sigana_dbs
rownames(row_sig_ana_dbs)[25] <- "CG>TT"
rownames(row_sig_ana_dbs)[26] <- "CG>GT"
rownames(row_sig_ana_dbs)[28] <- "CG>TC"
rownames(row_sig_ana_dbs)[46] <- "TA>GT"
rownames(row_sig_ana_dbs)[47] <- "TA>CT"
rownames(row_sig_ana_dbs)[49] <- "TA>GG"

row_sig_ana_dbs <- row_sig_ana_dbs[c(1:24,27,29,26,30,28,25,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,48,50,47,51,49,46,52:78),]
### fixed version looks like this:
match(rownames(row_sig_ana_dbs),rownames(dbs_mutpat_rename2))

######################### Save the matched DBS mutation matrices of MutationalPatterns and Siganalyzer
saveRDS(dbs_mutpat_rename2, file = "~/path/to/saving/matched/sbs_matrix/mutationalpatterns_dbs_countmatrices.RData")
saveRDS(row_sig_ana_dbs, file ="~/path/to/saving/matched/sbs_matrix/siganalyzer_dbs_countmatrices.Rdata" )

############################## Indel matching of Siganalyzer & MutationalPatterns#####################
indel_mutpat_rename <- indel_mutpat
colnames(indel_mutpat_rename)<- SP_num_SPname_do_meta2$name_sp[match(colnames(indel_mutpat_rename),SP_num_SPname_do_meta2$numeral)]
indel_mutpat_rename <- indel_mutpat_rename[!is.na(colnames(indel_mutpat_rename))]
indel_mutpat_rename <- as.matrix(indel_mutpat_rename)

#### set the matching identifiers based on cosine of snv
indel_mutpat_rename2 <- indel_mutpat_rename
for(d in list_do_dupli){
  cos_dub_snv <- cos_sim_matrix(sigana_mutationmatrix$snv[,match(d$name_sp,colnames(sigana_mutationmatrix$snv))],snv_mutpat_mat[,(match(d$numeral,colnames(snv_mutpat_mat)))])
  cos_dub_snv <- as.data.frame(cos_dub_snv)
  #  print(cos_dub_snv[1])
  # print(rownames(cos_dub_snv)[which(cos_dub_snv[,1] >=0.9999999)])
  print(cos_dub_snv[cos_dub_snv[,1] >=0.9999999])
  print(rownames(cos_dub_snv)[which(cos_dub_snv[,1] >=0.9999999)])
  # print(colnames(cos_dub_snv)[1])
  colnames(indel_mutpat_rename2)[match(rownames(cos_dub_snv)[1],colnames(snv_mutpat_rename3))] <- rownames(cos_dub_snv)[which(cos_dub_snv[,1] >=0.99999999)]
}
### get the matchin matrices of siganalyzer 
sigana_indel <- sigana_mutationmatrix$indel[match(colnames(indel_mutpat_rename2),colnames(sigana_mutationmatrix$indel))]

### save the correct matching Indel matrices of both signanalyzer and mutationalpatterns
saveRDS(indel_mutpat_rename2, file = "~/path/to/saving/matched/sbs_matrix/mutationalpatterns_indel_countmatrices.RData")
saveRDS(sigana_indel, file ="~/path/to/saving/matched/sbs_matrix/siganalyzer_indel_countmatrices.Rdata" )

##### Create a Dataframe without any duplicates and with the correct match based on sbs cosine's of the dupe samples #####
correct_match_df <- NULL
correct_match <- NULL
## create a dataframe with only the right matches of identifiets
for(d in list_do_dupli){
  cos_dub_snv <- cos_sim_matrix(sigana_mutationmatrix$snv[,match(d$name_sp,colnames(sigana_mutationmatrix$snv))],snv_mutpat_mat[,(match(d$numeral,colnames(snv_mutpat_mat)))])
  cos_dub_snv <- as.data.frame(cos_dub_snv)
  correct_match<- as.data.frame(SP_num_SPname_do_meta2[match(rownames(cos_dub_snv)[which(cos_dub_snv[,1]>0.999999)], SP_num_SPname_do_meta2$name_sp),])
  correct_match_df <- rbind(correct_match_df,correct_match)
  
}
### create a dataframe of the identifier code matches removing all duplicates
all_identifiers_correct <- SP_num_SPname_do_meta2[!SP_num_SPname_do_meta2$Do_id %in% all_do_multi,]
### add the correct versions of the duplicates back in 
all_identifiers_correct <- rbind(all_identifiers_correct, correct_match_df)
## end up with a dataframe with only the correct matches of samples without any duplicates
saveRDS(all_identifiers_correct, file = "~/path/to/folder/all_identifiercodes_correct.RDS")

### OUTPUT RESULTS
# 3 matrices (sbs, dbs, indel) of the PCAWG data as analysed by MutationalPatterns
# 3 matrices (sbs, dbs, indel) of the PCAWG data as analysed by SigAnalyzer
# 1 dataframe containing matching all the different identifier codes of the samples in PCAWG
