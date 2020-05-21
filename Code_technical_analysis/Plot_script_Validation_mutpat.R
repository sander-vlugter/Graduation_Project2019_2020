## Generation of plots and images for Validation of MutationalPatterns
## by: Sander V
## date: 16/May/2020
## Uses data created in other scripts
##
load("~/path/to/mutationalpatterns_snv_countmatrices.Rdata")
load("~/path/to/snv_sig_ana_matched.Rdata")

## get colsums
colsum_snv_mutpat <- colSums(mutationalpatterns_snv)
colsums_snv_sigana <- colSums(snv_sigana_matched) 
snv_comparison <- NULL
snv_comparison <- as.data.frame(cbind("mutpat"= colsum_snv_mutpat, "sigana"=colsums_snv_sigana))

library(ggplot2)
library(MutationalPatterns)
#### preview of plot without any Rsquare
snv_plot_difference<- ggplot(snv_comparison,aes(x=snv_comparison$mutpat,y=snv_comparison$sigana))+
  geom_point(shape=1)+
  geom_abline()+
  scale_y_continuous(trans = 'log10')+
  scale_x_continuous(trans = 'log10')+
  xlab("SBS Mutation load for each sample according to MutationalPatterns  (log-scale)")+
  ylab("SBS Mutation load for each sample according to SigAnalyzer  (log-scale)")+
  # ggtitle("SBS in the PCAWG database")+
  theme(panel.grid = element_line(colour = "gray"),
        panel.border = element_rect("black", fill = NA),
        panel.background = element_blank(), axis.title = element_text(size = 12))
plot(snv_plot_difference)

### get the R squared 

ttt<- lm(formula=mutpat ~ sigana,data = snv_comparison)
summary(ttt)$r.squared

snv_df <- snv_comparison

colnames(snv_df) <- c("x","y")
m <- lm(y ~ x, data = snv_df)

p_snv <- ggplot(data = snv_df, aes(x=x,y=y))+
  geom_smooth(method = "lm", formula=y~x)+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  geom_point()

p_snv


eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,
                 
                 list(        a = format(unname(coef(m)[1]), digits = 4),
                              
                              b = format(unname(coef(m)[2]), digits = 4),
                              
                              r2 = format(summary(m)$r.squared, digits = 3)))



dftext <- data.frame(x = 100, y = 1200000, eq = as.character(as.expression(eq)))

p_snv + geom_text(aes(label = eq), data = dftext, parse = TRUE)+
  theme(axis.title = element_text(size=14),plot.title = element_text(size=16))+
  xlab("SBS Mutation load for each sample according to MutationalPatterns  (log-scale)")+
  ylab("SBS Mutation load for each sample according to SigAnalyzer  (log-scale)")+
  theme(panel.grid = element_line(colour = "gray"),
        panel.border = element_rect("black", fill = NA),
        panel.background = element_blank())

############################### DBS plot unadjusted for MNV ##############

load("~/path/to/mutationalpatterns_dbs_countmatrices.Rdata")
load("~/path/to/siganalyzer_dbs_countmatrices.Rdata")

colsum_dbs_mutpat <- colSums(mutationalpatterns_dbs)
colsum_dbs_sigana <- colSums(dbs_sig_ana)

colsum_dbs_both <- as.data.frame(cbind("mutpat"=colsum_dbs_mutpat,"sigana"=colsum_dbs_sigana))

### preview of plot
dbs_diff_plot <- ggplot(colsum_dbs_both, aes(x=colsum_dbs_both$mutpat,y=colsum_dbs_both$sigana))+
  geom_point(shape=1)+
  geom_abline()+
  scale_y_continuous(trans = 'log10')+
  scale_x_continuous(trans = 'log10')+
  xlab("DBS Mutation load for each sample according to MutationalPatterns  (log-scale)")+
  ylab("DBS Mutation load for each sample according to SigAnalyzer  (log-scale)")+
  # ggtitle("SBS in the PCAWG database")+
  theme(panel.grid = element_line(colour = "gray"),
        panel.border = element_rect("black", fill = NA),
        panel.background = element_blank(), axis.title = element_text(size = 14))
plot(dbs_diff_plot)


######################################## DBS plot ADJUSTED for MNV ###################

load("~/path/to/dbs_diff_explained_matrice.Rdata") ## as created in "MutationalPatterns_create_true_dbs_matrix.R"

corrected_dbs <- dbs_diff_explained_table_category[,c(2,3,7)]
corrected_dbs$sigana <- corrected_dbs$sigana - corrected_dbs$would_be_DBSs
corrected_dbs<- corrected_dbs[,c(1,2)]
colnames(corrected_dbs)<- c("x","y")

m <- lm(y~x,data = corrected_dbs)

p_cor_dbs <- ggplot(corrected_dbs,aes(x=x,y=y))+
  geom_smooth(method = "lm", formula=y~x)+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  geom_point()
p_cor_dbs

eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,
                 list(        a = format(unname(coef(m)[1]), digits = 4),
                              b = format(unname(coef(m)[2]), digits = 4),
                              r2 = format(summary(m)$r.squared, digits = 5)))

dftext <- data.frame(x = 5, y = 9000, eq = as.character(as.expression(eq)))

p_cor_dbs + geom_text(aes(label = eq), data = dftext, parse = TRUE)+
  theme(plot.title = element_text(size=16))+
  xlab("Mutation load DBSs according to MutationalPatterns")+
  ylab("Mutation load DBSs according to SigAnalyzer (corrected for MNVs)")+
  theme(panel.grid = element_line(colour = "gray"),
        panel.border = element_rect("black", fill = NA),
        panel.background = element_blank(), axis.title = element_text(size = 14))

#################################### DBS percentile change plot ###############################

Dbs_percent_diff <- as.data.frame(table(dbs_diff_explained_table_category$category))
colnames(Dbs_percent_diff)[1]<- "range_percentage"

Dbs_percent_diff$perc_total <- (Dbs_percent_diff$Freq/ sum(Dbs_percent_diff$Freq))

Dbs_percent_diff$range_percentage <- c("0","0+ - 9.99","10 - 19.99","20 - 29.99","30 - 39.99","40 - 49.99",
                                       "50 - 59.99","60 - 69.99","70+")

DBS_percent_extracalled<- ggplot(Dbs_percent_diff, aes(x=Dbs_percent_diff$range_percentage, y=Dbs_percent_diff$perc_total,
                                                       label =round(Dbs_percent_diff$perc_total, digits = 4)))+
  geom_bar(stat = "identity", fill="steelblue")+
  geom_text(size=4, position = position_dodge(width = 0.9),vjust = -0.25 )+
  theme_minimal()+
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.551))+
  theme(axis.title = element_text(size = 14), plot.title = element_text(size = 16))+
  xlab("Percentage of the DBSs in a sample that are actually MNVs")+
  ylab("Relative number of the PCAWG samples")+
  theme(panel.grid = element_line(colour = "gray"),
        panel.border = element_rect("black", fill = NA),
        panel.background = element_blank(), axis.title = element_text(size = 14))

plot(DBS_percent_extracalled)

################################### Indel comparison plots ##########################################

load("~/path/to/mutationalpatterns_indel_countmatrices.Rdata")
load("~/path/to/siganalyzer_indel_countmatrices.Rdata")

colsum_indel_mutpat <- colSums(mutationalpatterns_indel)
colsum_indel_sigana <- colSums(indel_sigana_micro)

colsum_indel_both <- as.data.frame(cbind("mutpat"= colsum_indel_mutpat, "sigana"= colsum_indel_sigana))

indel_df <- colsum_indel_both

colnames(indel_df) <- c("x","y")
m <- lm(y ~ x, data = indel_df)

p_indel <- ggplot(data = indel_df, aes(x=x,y=y))+
  geom_smooth(method = "lm", formula=y~x)+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  geom_point()

p_indel


eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,
                 
                 list(        a = format(unname(coef(m)[1]), digits = 4),
                              
                              b = format(unname(coef(m)[2]), digits = 4),
                              
                              r2 = format(summary(m)$r.squared, digits = 5)))



dftext <- data.frame(x = 10, y = 7000, eq = as.character(as.expression(eq)))

p_indel + geom_text(aes(label = eq), data = dftext, parse = TRUE)+
  theme(axis.title = element_text(size=14),plot.title = element_text(size=16))+
  xlab("Indel Mutation load for each sample according to MutationalPatterns  (log-scale)")+
  ylab("Indel Mutation load for each sample according to SigAnalyzer  (log-scale)")+
  theme(panel.grid = element_line(colour = "gray"),
        panel.border = element_rect("black", fill = NA),
        panel.background = element_blank())

########################################### cosine ####################

library(MutationalPatterns)
##### SBS ######

matrix_sbs_mutpat <- as.matrix(mutationalpatterns_snv, drop =F)
matrix_sbs_sigana <- as.matrix(snv_sigana_2,drop =F)

cosi_sbs <- cos_sim_matrix(matrix_sbs_mutpat, matrix_sbs_sigana)
cosi_sbs_diag <- as.matrix(diag(cosi_sbs),drop =F)
### lowest cosine similarity between a sample being 0.999997


########## DBS cosine pre alterations to code (so all samples) ########

matrix_dbs_mutpat <- as.matrix(mutationalpatterns_dbs, drop=F)
matrix_dbs_sigana <- as.matrix(dbs_sig_ana, drop=F)

cosi_dbs_all_diag <- as.matrix(diag(cos_sim_matrix(matrix_dbs_mutpat,matrix_dbs_sigana)),drop=F)

length(which(colSums(matrix_dbs_mutpat)==0))
length(which(colSums(matrix_dbs_sigana)==0))
length(which(is.nan(cosi_dbs_all_diag)))

match(which(colSums(matrix_dbs_mutpat)==0),
      which(colSums(matrix_dbs_sigana)==0))
### nan is caused by values that have 0 mutations as defining cosine similarity for a vector of complete 0 results infinite
##  setting these to 1 as they are infact the same and come from the same samples.

cosi_dbs_all_diag_backup <- cosi_dbs_all_diag
cosi_dbs_all_diag_backup[is.nan(cosi_dbs_all_diag_backup)] <- 1

mean(cosi_dbs_all_diag_backup)
### averaging 0.9906 overall the samples



################################## Cosine DBS only exact samples (without MNV) ##################

same_dbs<- dbs_diff_explained_table_category[which(dbs_diff_explained_table_category$difference == 0),]


mutpat_dbs_equal<- mutationalpatterns_dbs[,match(same_dbs$sample, colnames(mutationalpatterns_dbs))]
sigana_dbs_equal<- dbs_sig_ana[,match(same_dbs$sample, colnames(dbs_sig_ana))]

mutpat_dbs_equal <- as.matrix(mutpat_dbs_equal)
sigana_dbs_equal <- as.matrix(sigana_dbs_equal)

cosi_dbs_adj_diag <- as.matrix(diag(cos_sim_matrix(mutpat_dbs_equal,sigana_dbs_equal)))
#### gets u cosine of 1 on average

dbs_context <- as.matrix(sigana_dbs_equal - mutpat_dbs_equal, drop =F)

any(dbs_context >0)
any(dbs_context <0)
################################### Cosine Indel ###############################


matrix_indel_mutpat <- as.matrix(mutationalpatterns_indel, drop=F)
matrix_indel_sigana <- as.matrix(indel_sigana_micro, drop = F)



cosi_indel_all_diag <-as.matrix(diag(cos_sim_matrix(matrix_indel_mutpat, matrix_indel_sigana)))

which(colSums(matrix_indel_mutpat)==0)
which(colSums(matrix_indel_sigana)==0)
which(is.nan(cosi_indel_all_diag)== T)
#### one sample has 0 indels called by both resulting in nan
cosi_indel_all_diag[2544]<- 1
### the others are caused by 1 of them being empty 
## the VCF's of these samples are empty (most likely corrupted)
## therefore will be removed before calculating the averga cosine
which(is.nan(cosi_indel_all_diag)== T)

cosi_indel_all_diag <- cosi_indel_all_diag[-c(which(is.nan(cosi_indel_all_diag)== T))]
### range of 0.57 to 1
mean(cosi_indel_all_diag)
## 0.996462
### now even thought the indel counts differ remember that the cosine isnt used for exact copy similiarity but rather,
## trend similarity


#################################### difference of contexts by subtracting ###################
## SBS 

sbs_subtract <- as.matrix(snv_sigana_2 - mutationalpatterns_snv, drop =F)
any(sbs_subtract > 0 ) # F
any(sbs_subtract < 0) # T
which(sbs_subtract <0)

sbs_subtract[ which(sbs_subtract <0) ]
## same was done for the DBS not containing MNVs 