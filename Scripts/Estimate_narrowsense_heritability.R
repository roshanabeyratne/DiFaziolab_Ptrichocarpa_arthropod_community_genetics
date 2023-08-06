################################################################################
# 
### Estimate narrow-sense heritability for arthropod community indices for 
### Clatskanie and Corvallis trials of Populus trichocarpa
### 06-11-2023, (Updated: 08-06-2023) RoshanAbeyratne DavidMacayaSanz & StephenDiFazio
# 
################################################################################
rm(list=ls())
workdir <- '/Volumes/WD_CRA/CRA_Projects/CRA_Arthropod_community_genetics/GitHub_upload'
setwd(workdir)
library(tidyverse)
library(sommer)
options(stringsAsFactors = FALSE)

# Read input phenotypes
Pheno_clats <- read.table('./Data/Pheno_clats.txt', header=T)
Pheno_corv <- read.table('./Data/Pheno_corv.txt', header=T)

# Read input genomic relationship matrices
GRM_clats <- read.table('./Data/GRM_clats.txt', header=T)
colnames(GRM_clats) <- rownames(GRM_clats) <- 
  colnames(GRM_clats) %>% str_replace("^X", '')
GRM_clats <- as.matrix(GRM_clats)
GRM_corv <- read.table('./Data/GRM_corv.txt', header=T)
colnames(GRM_corv) <- rownames(GRM_corv) <- 
  colnames(GRM_corv) %>% str_replace("^X", '')
GRM_corv <- as.matrix(GRM_corv)


# # Estimate narrow-sense heritability using the individual tree model
# h2_df <- NULL
# for (trial in c('clats', 'corv')) {
#   # trial <- 'clats'
#   for (Trait in c('abd', 'div', 'rich')) {
#     # Trait <- 'abd'
#     data <- get(paste0('Pheno_', trial))
#     GRM <- get(paste0('GRM_', trial))
#     tmp_data <- data[,c('seq_name', paste0(trial, '_', Trait))]
#     colnames(tmp_data)[2] <- 'Trait'
# 
#     ### Individual-tree-Model
#     indv_mod <- sommer::mmer(Trait~1,
#                   random=~ vs(seq_name, Gu=GRM),
#                   rcov=~units,
#                   data=tmp_data)
#     sum_indv_mod <- summary(indv_mod)
#     var_comp <- sum_indv_mod$varcomp
#     h2 <- vpredict(indv_mod, h2 ~ V1/(V1+V2))
#     h2_df <- rbind(h2_df, c(trial, Trait, h2$Estimate, h2$SE))
#   }
# }
# # colnames(h2_df) <- c('Trial', 'Triat', 'h2', 'h2_se')
# # h2_df <- as.data.frame(h2_df)
# # h2_df$h2 <- as.numeric(h2_df$h2)
# # h2_df$h2_se <- as.numeric(h2_df$h2_se)
# # (h2_df)
# # Trial Triat        h2     h2_se
# # clats   abd 0.2536612 0.2273321
# # clats   div 0.0000000 0.1410154
# # clats  rich 0.2436657 0.1436884
# # corv   abd 0.0000000 0.1901620
# # corv   div 0.3388757 0.3008195
# # corv  rich 0.5925099 0.3345800
# # NOTE: Since h2 is zero for certain traits, we will be carrying out permutation
# # just for h2 > 0


# # create permuted dataset for abd, div, rich (all datasets are permuted
# # irrespective of whether h2 zero or not)
# npmt <- 10000
# for (trial in c('clats', 'corv')) {
#   # trial <- 'clats'
#   cat("\n", trial)
#   for (Trait in c('abd', 'div', 'rich')) {
#     # Trait <- 'abd'
#     cat("\t", Trait)
#     data <- get(paste0('Pheno_', trial))
#     tmp_data <- data[,c('seq_name', paste0(trial, '_', Trait))]
#     for (p in 1:npmt) {
#       orig <- as.numeric(as.matrix(tmp_data[,2]))
#       tmp_pmt <- sample(orig, size=length(orig), replace=FALSE)
#       tmp_data <- cbind(tmp_data, tmp_pmt)
#     }
#     colnames(tmp_data)[2:(npmt+2)] <- paste0(paste0(trial, '_', Trait),'_', 0:npmt)
#     write.table(tmp_data, file=paste0('./Data/', trial, '_', Trait, '_with_10k_pmut.txt'))
#   }
# }


# Estimate h2 for each trial and phenotype which carries non-zero h2
trait_key <- data.frame(matrix(c(rep('clats', 2),
                                 rep('corv', 2),
                                 'abd', 'rich',
                                 'div', 'rich'),
                               ncol=2, byrow=F))
colnames(trait_key) <- c('trial', 'Trait')
# NOTE: Commenting out this section of the code as this may take a long time to run.
# Output files are saved to the data folder
# for (i in 1:nrow(trait_key)) {
#   # i <- 1
#   trial <- trait_key$trial[i]
#   Trait <- trait_key$Trait[i]
#   tmp_data <- read.table(file=paste0('./Data/permutations/',
#                                      trial, '_', Trait,
#                                      '_with_10k_pmut.txt'),
#                          header=TRUE)
#   GRM <- get(paste0('GRM_', trial))
# 
#   # start_time <- Sys.time()
#   h2_df <- NULL
#   npmt <- 10000
#   for (p in 0:npmt) {
#     cat("\t", p)
#     tmp_data_ss <- tmp_data[,c('seq_name', paste0(trial, '_', Trait, '_', p))]
#     colnames(tmp_data_ss)[2] <- 'Trait'
#     indv_mod <- sommer::mmer(Trait~1,
#                              random=~ vs(seq_name, Gu=GRM),
#                              rcov=~units,
#                              data=tmp_data)
#     sum_indv_mod <- summary(indv_mod)
#     var_comp <- sum_indv_mod$varcomp
#     h2 <- vpredict(indv_mod, h2 ~ V1/(V1+V2))
#     h2_df <- rbind(h2_df, c(trial, Trait, p, h2$Estimate, h2$SE))
#   }
#   # end_time <- Sys.time()
#   # end_time-start_time
#   h2_df <- as.data.frame(h2_df)
#   colnames(h2_df) <- c('trial', 'Trait', 'iter', 'h2', 'h2_se')
#   h2_df$h2 <- as.numeric(h2_df$h2)
#   h2_df$h2_se <- as.numeric(h2_df$h2_se)
#   write.table(h2_df, file=paste0('./data/h2_out/', trial, '_', Trait, '_h2_df.txt'))
# }

# h2 for all traits
herit_p <- NULL
for (i in 1:nrow(trait_key)) {
  # i=1
  trial <- trait_key$trial[i]
  Trait <- trait_key$Trait[i]
  p_h2 <- read.table(paste0('./Data/h2_out/', trial, '_', Trait, '_h2_df.txt'),
                     header=TRUE)
  emp_p <- sum(p_h2$h2 > p_h2$h2[1])/length(p_h2$h2)
  herit_p <- rbind(herit_p, c(p_h2[1, c(1:2, 4)], emp_p))
}
herit_p <- data.frame(herit_p)
colnames(herit_p)[4] <- 'p_permuted'
(herit_p)
# trial Trait        h2 p_permuted
# clats   abd 0.2536612  0.0989901
# clats  rich 0.2436657 0.09639036
# corv   div 0.3388757 0.08749125
# corv  rich 0.5925099 0.02049795
################################################################################
# END
# R version 3.6.3 (2020-02-29)

# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
#
# other attached packages:
#   [1] sommer_4.1.5    crayon_1.4.1    lattice_0.20-44 MASS_7.3-54     Matrix_1.3-3    forcats_0.5.1   stringr_1.4.0   dplyr_1.0.8    
# [9] purrr_0.3.4     readr_1.3.1     tidyr_1.2.0     tibble_3.1.2    ggplot2_3.1.1   tidyverse_1.2.1
#
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.6       cellranger_1.1.0 pillar_1.6.1     compiler_3.6.3   plyr_1.8.6       tools_3.6.3      lubridate_1.7.4  jsonlite_1.7.2  
# [9] lifecycle_1.0.1  nlme_3.1-152     gtable_0.3.0     pkgconfig_2.0.3  rlang_1.0.2      rstudioapi_0.13  DBI_1.0.0        cli_2.5.0       
# [17] yaml_2.2.1       haven_2.1.0      withr_2.4.2      xml2_1.3.2       httr_1.4.2       generics_0.0.2   vctrs_0.3.8      hms_0.4.2       
# [25] grid_3.6.3       tidyselect_1.1.1 glue_1.4.2       R6_2.5.0         fansi_0.4.2      readxl_1.3.1     modelr_0.1.8     magrittr_2.0.1  
# [33] backports_1.2.1  scales_1.1.1     ellipsis_0.3.2   rvest_0.3.3      colorspace_1.4-1 utf8_1.2.1       stringi_1.6.2    lazyeval_0.2.2  
# [41] munsell_0.5.0    broom_0.5.2  
