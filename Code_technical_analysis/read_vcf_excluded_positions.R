read_vcf_excluded_positions <-function (vcf_files, sample_names, genome, group = "auto+sex", 
          check_alleles = TRUE, n_cores) 
{
  if (length(vcf_files) != length(sample_names)) 
    stop("Please provide the same number of sample names as VCF files")
  ref_genome <- base::get(genome)
  ref_organism <- GenomeInfoDb::organism(ref_genome)
  ref_style <- seqlevelsStyle(ref_genome)
  genome_name <- genome(ref_genome)[[1]]
  if (!(class(ref_genome) == "BSgenome")) 
    stop("Please provide the name of a BSgenome object.")
  if (missing(n_cores)) {
    n_cores = detectCores()
    if (!(.Platform$OS.type == "windows" || is.na(n_cores))) 
      n_cores <- detectCores()
    else n_cores = 1
  }
  original_warn_state = getOption("warn")
  options(warn = -1)
  warnings <- NULL
  if (!check_alleles) {
    warning(paste("check_alleles is set to FALSE.  Make sure your", 
                  "input VCF does not contain any positions with", 
                  "multiple alternative", "alleles, as these positions cannot be analysed", 
                  "with MutationalPatterns and cause obscure", "errors."), 
            immediate. = TRUE)
  }
  vcf_list <- mclapply(seq_along(vcf_files), function(index) {
    file <- vcf_files[index]
    vcf <- rowRanges(readVcf(file, genome_name))
    if (length(vcf) == 0) 
      stop(sprintf("Vcf file %s is empty", file))
    seqlevelsStyle(vcf) <- ref_style[1]
    groups <- c()
    if (group != "none") {
      if (group == "auto+sex") {
        groups <- c(extractSeqlevelsByGroup(species = ref_organism, 
                                            style = ref_style, group = "auto"), extractSeqlevelsByGroup(species = ref_organism, 
                                                                                                        style = ref_style, group = "sex"))
        groups_names <- names(groups)
        if (!is.null(groups_names)) {
          unique_names <- unique(groups_names)
          groups <- llply(unique_names, function(x) groups[groups_names == 
                                                             x])
          groups <- llply(groups, unlist, recursive = FALSE)
          groups <- unique(as.vector(groups[[1]]))
        }
      }
      else {
        groups <- extractSeqlevelsByGroup(species = ref_organism, 
                                          style = ref_style, group = group)
        groups <- unique(as.vector(t(groups)))
      }
      groups <- intersect(groups, seqlevels(vcf))
      vcf <- keepSeqlevels(vcf, groups, pruning.mode = "tidy")
    }
    if (check_alleles) {
      rem <- which(lengths(vcf$ALT) > 1)
      if (length(rem) > 0) {
        vcf = vcf[-rem]
        warnings$check_allele <- rbind(warnings$check_allele, 
                                       c(sample_names[[index]], length(rem)))
      }
    }
    dbs = intersect(which(diff(start(vcf)) == 1), which(nchar(as.character(vcf$REF)) == 
                                                          1 & nchar(as.character(unlist(vcf$ALT))) == 1))
    if (length(dbs) > 0) 
      if (dbs[length(dbs)] == length(vcf)) 
        dbs <- dbs[-length(dbs)]
    mnv = NULL
    if (1 %in% diff(dbs)) {
      rem = which(diff(dbs) == 1)
      mnv = unique(sort(c(dbs[rem], dbs[rem] + 1, dbs[rem] + 
                            2)))
      rem <- unique(c(rem, rem + 1))
      dbs <- dbs[-rem]
      warnings$dbs <- rbind(warnings$dbs, c(sample_names[[index]], 
                                            length(mnv),length(which(diff(start(vcf)[mnv])==1))))
    }else{
      warnings$dbs <- rbind(warnings$dbs, c(sample_names[[index]],
                                            0,0)) ### addition 
    }
    if (length(dbs) > 0) {
      end(vcf)[dbs] = start(vcf)[dbs] + 1
      vcf$REF[dbs] = DNAStringSet(paste0(as.character(vcf$REF[dbs]), 
                                         as.character(vcf$REF[dbs + 1])))
      vcf$ALT[dbs] = DNAStringSetList(lapply(dbs, function(i) {
        DNAStringSet(paste0(as.character(unlist(vcf$ALT[i])), 
                            as.character(unlist(vcf$ALT[i + 1]))))
      }))
      vcf = vcf[-c(mnv, (dbs + 1)), ]
    }
    indel = which((nchar(as.character(vcf$REF)) == 1 & nchar(as.character(unlist(vcf$ALT))) > 
                     1) | (nchar(as.character(vcf$REF)) > 1 & nchar(as.character(unlist(vcf$ALT))) == 
                             1))
    if (any(c(grepl("[^ACGT]", as.character(vcf$REF[indel])), 
              grepl("[^ACGT]", as.character(unlist(vcf$ALT[indel])))))) {
      vcf = vcf[-indel, ]
      warnings$indel <- rbind(warnings$indel, c(sample_names[[index]]))
    }
    return(list(vcf, warnings))
  }, mc.cores = n_cores)
  options(warn = original_warn_state)
  summ <- lapply(vcf_list, function(item) {
    if (class(item) == "try-error") 
      stop(item)
    ref = as.character(item[[1]]$REF)
    alt = as.character(unlist(item[[1]]$ALT))
    nsnvs = length(which(nchar(ref) == 1 & nchar(alt) == 
                           1))
    ndbs = length(which(nchar(ref) == 2 & nchar(alt) == 2))
    nindel = length(which((nchar(ref) == 1 & nchar(alt) != 
                             1) | (nchar(ref) != 1 & nchar(alt) == 1)))
    return(c(nsnvs, ndbs, nindel))
  })
  warnings <- do.call(rbind, vcf_list)[, 2]
  warnings <- sapply(warnings, function(item){
      if (is.null(item)) 
        item = list("check_alleles"=NULL,"dbs"=NULL,"indel"=NULL)
      item[c("check_alleles", "dbs", "indel")]
    })
  warns = NULL
  for (i in which(!(is.na(rownames(warnings))))) {
    warns[[rownames(warnings)[i]]] = do.call(rbind, warnings[i, 
                                                             ])
    colnames(warns[[rownames(warnings)[i]]]) <- c("Sample", 
                                                  "Position(s)","wouldbe_DBSs")
  }
  lapply(names(warns), function(item) {
    if (item == "check_alleles") {
      warning("Position(s) with multiple ", "alternative alleles are excluded\n", 
              paste0(capture.output(warns$check_allele), collapse = "\n"), 
              immediate. = TRUE)
    }
    else if (item == "dbs") {
      warning("Position(s) excluded that form ", "multiple nucleotide variants ", 
              "of length more than 2.\n", paste0(capture.output(warns$dbs), 
                                                 collapse = "\n"), immediate. = TRUE)
    }
    else if (item == "indel") {
      warning("Indels not according to VCF format ", "4.2 or higher, all indels are excluded from:\n ", 
              paste(warns$indel, collapse = ", "), immediate. = TRUE)
    }
  })
  vcf_list <- GRangesList(do.call(rbind, vcf_list)[, 1])
  names(vcf_list) <- sample_names
  summ <- do.call(rbind, summ)
  colnames(summ) <- c("Number of SNV", "Number of DBS", "Number of indel")
  rownames(summ) <- names(vcf_list)
  print(summ)
  ### list of vcf list and warns$dbs
  warns_n_vcf <- list(vcf_list, warns$dbs)
  return(warns_n_vcf)
}
