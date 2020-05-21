#### PCAWG VCF data - creating mutation matrices for SBS DBS & Indels #########
## Author: Sander Vlugter
## of importance: this script is an adaptation of code provided to me by the Cuppen group
## Date (final version): 15/05/2020
## Description:
## creating mutation matrices of PCAWG samples by MutationalPatterns
## OUTPUTS: 3 matrices for the samples of the PCAWG database
# one matrix containing the SBS for the samples
# one matrix containing the DBS for the samples
# one matrix containing the Indels for the samples

#### Library of packages needed ####
## The installing steps for packages are not included.
library(devtools)
# For MutationalPatterns, a developmental version is used which can be found at:
# https://github.com/UMCUGenetics/MutationalPatterns/tree/dbs-indels
library(MutationalPatterns)
library(BSgenome)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(foreach)
library(doParallel)

#### Loading in the PCAWG SBS/DBS data for MutationalPatterns ####
# As the for every sample in our version of the PCAWG data consists of of 2 VCF files additional steps are taken
# For the sake of time the VCFs of SBS/DBS of every sample are loaded and processed by MutationalPatterns individually
# resulting in an output of .txt files for each sample of both DBS and SBS



# NOTE: that this produces the right number of DBS calls as by MutationalPatterns, as the older version of MutationalPatterns that called MNVs as DBSs no longer exists
# but does not report the number of MNVs in a sample, the script that reports the number of MNVs see: "MutationalPatterns_vcf_excluded_mnvs_mutationmatrix"

path_sbs_vcfs = "~/path/to/pcawg/sbs_dbs/vcfs"

# list the VCFs
vcf_files_SBS <- list.files(path_sbs_vcfs, pattern = "PASS.vcf.gz$",recursive = T,full.names = T)
# set output directory
output = "~/path/to/output/"
# set the sample names 
sample_names <- do.call(rbind, strsplit(vcf_files_SBS, "mnv/"))[,2]
sample_names <- do.call(rbind, strsplit(sample_names, "\\."))[,1]
# select reference genome & set type as required by MutationalPatterns
ref_genome = "BSgenome.Hsapiens.1000genomes.hs37d5"
type = "all"

print("Process vcf files")

registerDoParallel(cores = 16)

print(length(vcf_files_SBS))

foreach (i = seq_along(vcf_files_SBS)) %dopar% {
  print(i)
  vcf <- read_vcfs_as_granges(vcf_files_SBS[i], sample_names[i], ref_genome, check_alleles = T)
  mut_matrix <- mut_matrix(vcf, ref_genome, type = c("snv","dbs"))
  
  if (class(mut_matrix) == "list"){
    write.table(mut_matrix$snv, file = sprintf("%s/SNV/%s.txt", output, sample_names[i]), quote = F )
    write.table(mut_matrix$dbs, file = sprintf("%s/DBS/%s.txt", output, sample_names[i]), quote = F )
  } else if (class(mut_matrix) == "matrix"){
    if (rownames(mut_matrix)[1] == "A[C>A]A"){
      write.table(mut_matrix, file = sprintf("%s/SNV/%s.txt", output, sample_names[i]), quote = F )
    } else if (rownames(mut_matrix)[1] == "AC>CA"){
      write.table(mut_matrix, file = sprintf("%s/DBS/%s.txt", output, sample_names[i]), quote = F )
    } else if (rownames(mut_matrix)[1] == "del.1bp.homopol.C.len.1"){
      write.table(mut_matrix, file = sprintf("%s/INDEL/%s.txt", output, sample_names[i]), quote = F)
    }
  }
}


############### Indel VCFs of the PCAWG samples ######################
# Loading in the PCAWG INDEL data for MutationalPatterns 
# As the for every sample in our version of the PCAWG data consists of of 2 VCF files additional steps are taken
# For the sake of time the VCFs of Indel of every sample are loaded and processed by MutationalPatterns individually
# resulting in an output of .txt files for each sample for INDELS
path_indel_vcfs = "~/path/to/pcawg/indel/vcfs"

# list the VCFs
vcf_files_indels <- list.files(path_indel_vcfs, pattern = "PASS.vcf.gz$",recursive = T,full.names = T)
# set output directory
output = "~/path/to/output/"
# set the sample names 
sample_names <- do.call(rbind, strsplit(vcf_files_indels, "mnv/"))[,2]
sample_names <- do.call(rbind, strsplit(sample_names, "\\."))[,1]
# select reference genome & set type as required by MutationalPatterns
ref_genome = "BSgenome.Hsapiens.1000genomes.hs37d5"
type = "all"

print("Process vcf files")

registerDoParallel(cores = 16)

print(length(vcf_files_indels))

foreach (i = seq_along(vcf_files_indels)) %dopar% {
  print(i)
  vcf <- read_vcfs_as_granges(vcf_files_indels[i], sample_names[i], ref_genome, check_alleles = T)
  mut_matrix <- mut_matrix(vcf, ref_genome, type = c("snv","dbs"))
  
  if (class(mut_matrix) == "list"){
    write.table(mut_matrix$snv, file = sprintf("%s/SNV/%s.txt", output, sample_names[i]), quote = F )
    write.table(mut_matrix$dbs, file = sprintf("%s/DBS/%s.txt", output, sample_names[i]), quote = F )
  } else if (class(mut_matrix) == "matrix"){
    if (rownames(mut_matrix)[1] == "A[C>A]A"){
      write.table(mut_matrix, file = sprintf("%s/SNV/%s.txt", output, sample_names[i]), quote = F )
    } else if (rownames(mut_matrix)[1] == "AC>CA"){
      write.table(mut_matrix, file = sprintf("%s/DBS/%s.txt", output, sample_names[i]), quote = F )
    } else if (rownames(mut_matrix)[1] == "del.1bp.homopol.C.len.1"){
      write.table(mut_matrix, file = sprintf("%s/INDEL/%s.txt", output, sample_names[i]), quote = F)
    }
  }
}

#### outputs thus far :
# for each sample 3 seperate .txt files containing SBS DBS and Indels

# Creating the mutation matrices as by MutationalPatterns from the .txt files #
####### function creating the mutation matrices from .txt files. ######
loadContextVectors <- function(path){
  files <- list.files(path,pattern = '*.txt', full.names=T)
  
  message("Loading txt files")
  pb <- txtProgressBar(min=0, max=length(files), initial=0, style=3) 
  counter <- 0
  l <- lapply(files, function(i){ 
    counter <<- counter + 1
    setTxtProgressBar(pb, counter)
    
    tryCatch({
      read.table(i, sep=' ', check.names=F, comment.char = "")
    }, error=function(err){ 
      df_error <- data.frame(NA)
      names(df_error) <- sub('[.]txt','',basename(i))
      return(df_error)
    })
    
  })
  
  do.call(cbind, l)

}
################################
# set path to the mutation type SBS .txt files
path_to_sbs_txt = "~/path/to/sbstxt_files/"
# use path to create the MutationalPatterns SBS matrix of PCAWG samples
mutpat_sbs_matrix<- loadContextVectors(path_to_sbs_txt)

# set path to the mutation type DBS .txt files
path_to_dbs_txt = "~/path/to/dbstxt_files/"
# use path to create the MutationalPatterns DBS matrix of PCAWG samples
mutpat_dbs_matrix<- loadContextVectors(path_to_dbs_txt)

# set path to the mutation type Indel .txt files
path_to_indel_txt = "~/path/to/indeltxt_files/"
# use path to create the MutationalPatterns Indel matrix of PCAWG samples
mutpat_indel_matrix<- loadContextVectors(path_to_indel_txt)

#### Save these matrices for use in other scripts

saveRDS(mutpat_sbs_matrix,file = "~/path/to/save_location/sbs_mutpat.RData")
saveRDS(mutpat_dbs_matrix,file = "~/path/to/save_location/dbs_mutpat.RData")
saveRDS(mutpat_indel_matrix,file = "~/path/to/save_location/indel_mutpat.RData")



