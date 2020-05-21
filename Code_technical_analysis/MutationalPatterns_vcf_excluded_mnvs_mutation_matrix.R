### Title: MutationalPatterns extract 'true' DBS from VCF report MNV
## Author: Sander Vlugter
## Date (final version): 17/May/2020
## Description: Adaptation of MutationalPatterns function: "read_vcf_as_Granges", to find and report MNV's.
## this adaptation reports the number of MNV's in a sample & reports the number of DBSs SigAnalyzer would report as a result of this.
## this adaption of the function is called: "read_vcf_excluded_positions".

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
## Important note: first load "read_vcf_excluded_positions" into the global environment
source("~/path/to/read_vcf_excluded_positions.R")

#### Loading in the PCAWG SBS/DBS data for MutationalPatterns ####
## this is the version that reports the number of MNVs in a sample and reports the number of DBSs would be called by SigAnalyzer
path_sbs_vcfs = "~/path/to/pcawg/sbs_dbs/vcfs"

# list the VCFs
vcf_files_SBS <- list.files(path_sbs_vcfs, pattern = "PASS.vcf.gz$",recursive = T,full.names = T)
# set output directory
output = "~/path/to/output/truedbs/"
# set the sample names 
sample_names <- do.call(rbind, strsplit(vcf_files_SBS, "mnv/"))[,2]
sample_names <- do.call(rbind, strsplit(sample_names, "\\."))[,1]
# select reference genome & set type as required by MutationalPatterns
ref_genome = "BSgenome.Hsapiens.1000genomes.hs37d5"
type = "all"

print("Process vcf files")

registerDoParallel(cores = 16)

print(length(vcf_files_SBS))
### create an empty dataframe with 3 columns
# for sample id, number of excluded positions (the mnv positions), would be dbs called by siganalyzer
emptdf <-setNames(data.frame(matrix(ncol=3, nrow=0)),c("sample_name","excluded_positions","would_be_DBSs"))
write.table(emptdf, file =sprintf("%s/excluded_pos/excluded_pos.txt",output))

foreach (i = seq_along(vcf_files_SBS)) %dopar% {
  print(i)
  #vcf <- read_vcfs_as_granges(vcf_files[i], sample_names[i], ref_genome, check_alleles = T)
  vcf <- read_vcf_excluded_positions(vcf_files_SBS[i], sample_names[i], ref_genome, check_alleles = T)
  mut_matrix <- mut_matrix(vcf[[1]], ref_genome, type = c("snv","dbs"))
  
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
  write.table(vcf[[2]], file =sprintf("%s/excluded_pos/excluded_pos.txt",output) ,append = T, quote = F,col.names = F, row.names = F)
}
#### Outputs:
# SBS mutation matrix as .txt file for each sample
# DBS mutation matrix as .txt file for each sample
# one .txt file consisting of the number of MNVs in a sample as well as the number of DBSs would be called by SigAnalyzer
