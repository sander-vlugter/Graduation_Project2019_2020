### Title: MutationalPatterns 'true' DBS mutation matrix
## Author: Sander Vlugter
## Date (final version): 17/May/2020
## Description: 
## After running "MutationalPatterns_vcf_excluded_mnvs_mutation_matrix", we have the following outputs:
## for each sample SBS mutationmatrix as .txt file
## for each sample 'true' DBS mutationmatrix as .txt file (true DBS meaning, MNVs are excluded)
## one .txt matrix containing information on all the samples, how many MNV positions there are excluded and
## the number of "would be DBSs" would be called as a result of the MNVs by SigAnalyzer.
########
## Here we create the true DBS mutationmatrix by MutationalPatterns, and correct SigAnalyzers version for number of "would be DBSs".
## outputs:
## true dbs matrix by MutationalPatterns
## plot for supplementary
## dataframe showing and explaining the differences caused by MNVs
#### Library
library(MutationalPatterns)
## function to create mutation matrix of of individual .txt versions of the matrices
loadContextVectors <- function(path){
  files <- list.files(path,pattern = '*.txt', full.names=T)
  
  message("Loading txt files")
  pb <- txtProgressBar(min=0, max=length(files), initial=0, style=3) 
  counter <- 0
  l <- lapply(files, function(i){ 
    counter <<- counter + 1
    setTxtProgressBar(pb, counter)

    tryCatch({
      # print(i)
      read.table(i, sep=' ', check.names=F, comment.char = "")
    }, error=function(err){ 
      df_error <- data.frame(NA)
      names(df_error) <- sub('[.]txt','',basename(i))
      return(df_error)
    })
    
  })
  
  do.call(cbind, l)
}
#### set path to location of the 'true' dbs .txt files 
path = "~/path/to/output/truedbs/DBS/"
mutationalpatterns_dbs_countmatrices <- loadContextVectors(path = path)
### load in the file with the correct identifiers, as created in "Match_identifiers_mutpat_sigana.R"
all_identifiers_correct <- readRDS("~/path/to/folder/all_identifiercodes_correct.RDS")

## rename the sample idfentiers of the true DBS mutation count matrix by MutationalPatterns
colnames(mutationalpatterns_dbs_countmatrices) <- all_identifiers_correct$name_sp[match(colnames(mutationalpatterns_dbs_countmatrices),
                                                                                        all_identifiers_correct$numeral)]
## remove samples with missing identifiers
mutationalpatterns_dbs_countmatrices <- mutationalpatterns_dbs_countmatrices[!is.na(colnames(mutationalpatterns_dbs_countmatrices))]
## load in the DBS mutationmatrix created by SigAnalyzer, as generated in: "Match_identifiers_mutpat_sigana.R"
sigana_dbs <- readRDS("~/path/to/saving/matched/sbs_matrix/siganalyzer_dbs_countmatrices.Rdata")

## create a dataframe of 2 columns displaying the DBS mutation load of each sample individualy as produced by both tools 
colsum_of_both_dbs <-as.data.frame(cbind("mutpat"=colSums(dbs_mutpat_rename2), "sigana"=colSums(sigana_dbs)))
## add column of the difference in numbers of DBS of the 2 tools
colsum_of_both_dbs$difference <- colsum_of_both_dbs$sigana - colsum_of_both_dbs$mutpat
## add identifier code as a column:
colsum_of_both_dbs$sample <- rownames(colsum_of_both_dbs)
### load in the excluded MNV positions file, created in: "MutationalPatterns_vcf_excluded_mnvs_mutation_matrix.R"
excluded_dbs <- read.delim("~/path/to/output/truedbs/excluded_pos/excluded_pos.txt", sep ="")
## create a dataframe matching the indetier codes of the samples in the excluded MNV positions file
excluded_dbs <- cbind(excluded_dbs,"sp_name"=all_identifiers_correct$name_sp[match(excluded_dbs$sample_name, all_identifiers_correct$numeral)])
## remove any missing values
excluded_dbs <- excluded_dbs[!is.na(excluded_dbs$sp_name),]
colnames(excluded_dbs)[4]<- "sample"

### create a dataframe of the number of DBS as called by both tools & excluded mnvs
dbs_diff_explained_table <- merge(colsum_of_both_dbs,excluded_dbs)
### shows the different number of DBSs called for each sample by both tools
### shows the number of excluded positions by MutationalPatterns due to these positions being MNVs
### shows the number of Would be DBSs called by SigAnalyzer
### the difference in number of DBSs called is the same number as the would be DBS calls by SigAnalyzer

### add rownames
rownames(dbs_diff_explained_table)<- dbs_diff_explained_table$sample
### calculate percentage difference in number of DBS calls
dbs_diff_explained_table<-cbind(dbs_diff_explained_table,"percent"=((dbs_diff_explained_table$sigana -dbs_diff_explained_table$mutpat)/ dbs_diff_explained_table$sigana)*100)
### certain samples simply dont have DBSs and or MNVs those should be set to 0%
dbs_diff_explained_table$percent[is.nan(dbs_diff_explained_table$percent)] <- 0
### create percentage categories
dbs_diff_explained_table_category <- dbs_diff_explained_table
dbs_diff_explained_table_category$category <- cut(dbs_diff_explained_table_category$percent,
                                                  breaks = c(0,0.00001,10,20,30,40,50,60,70,Inf),
                                                  labels = c("0","0+ - 9.999","10 - 19.999","20 - 29.999","30 - 39.999","40 - 49.999","50 - 59.999","60 - 69.999","70+"),
                                                  right = F)


### save the output of the diff_explained_table_category for later use
saveRDS(dbs_diff_explained_table_category, "~/path/to/folder/dbs_diff_explained_table.RDS")