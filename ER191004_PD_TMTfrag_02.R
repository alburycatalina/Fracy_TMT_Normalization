# Functions for TMT IRS Normalization with a variable number of MS experiments

library(edgeR)
library(ggplot2)
library(dplyr)
library(magrittr)
library(testthat)
library(pastecs)
library(rpgm)

#load and subset data

PD1 <- read.delim('190805_865_097_TMTfrag_Proteins.txt')

PD1_names <- read.delim('new_col_names.txt')

#re-order columns in dataset so that accession and quant are front and centre
PD1_sub <- PD1[c(4,19:38,2,5:18)]

names(PD1_sub) <- gsub(pattern = "Abundances.Grouped.",replacement = "",x = names(PD1_sub),fixed = TRUE)

PD1_sub2 <- PD1_sub

#check that re-naming of column headers was correct before replacing column headers (there may be a more elegant way to do this)
colnames(PD1_sub2)
colnames(PD1_names)

colnames(PD1_sub2) <- colnames(PD1_names)

#Re-order columns to give: accession, TMTsetA, TMTsetB, others. empty channels removed
PD1_sub3 <- PD1_sub2 [c(1,2,4,6,8,10,12,14,20,3,5,7,9,11,13,15,17,19,21,22:36)]

#replace zeros in data with NAs

PD1_sub4A <- PD1_sub3[c(2:19)]
PD1_sub4B <- PD1_sub3[c(1,20:34)]

PD1_sub4A[PD1_sub4A == 0] <- NA

#Add back in columns from PD1_sub3 df

PD1_sub4 <- cbind(PD1_sub4A, PD1_sub4B)

PD1_sub4 <- PD1_sub4[c(19,1:18,20:34)]

colnames(PD1_sub4)

#basic stats for raw table

write.csv(PD1_sub4, file = "PD1_raw.csv")

options(scipen=100)
options(digits=2)
stats_PD1_sub4 <- stat.desc(PD1_sub4)

write.csv(stats_PD1_sub4, file = "stats_PD1_raw.csv")

#Reomove or filter low count experiments, experiments 13,14,15,22,23,24 (TMT channels C3,D3,E3,B5,C5,D5) removed
PD1_sub4_filt <- PD1_sub4[c(1,2:15,17:23,25,27:33,35,37:43,45:66)]

PD1_sub4_filt2 <- PD1_sub4[c(1,2:15,17:23,25,27,28,30:33,35,37:43,45:66)]

#check the correct experiments were removed
colnames(PD1_sub4_filt)
colnames(PD1_sub4_filt2)

# to delete and impute rows with greater than n NAs. and posibly replace (CAUTION with this!!!)
delete.na <- function(DF, n=0) {
DF[rowSums(is.na(DF)) <= n,]
}

PD1_sub4_filt_22nas<- delete.na(PD1_sub4_filt, 22)

PD1_sub4_filt2_22nas<- delete.na(PD1_sub4_filt2, 22)

write.csv(PD1_sub4_filt_22nas, file = "PD1_sub4_filt_22nas.csv")
write.csv(PD1_sub4_filt2_22nas, file = "PD1_sub4_filt2_22nas.csv")


#Remove incomplete rows for normalization

PD1_sub5 <- PD1_sub4[complete.cases(PD1_sub4),]

PD1_sub5_filt <- PD1_sub4_filt[complete.cases(PD1_sub4_filt),]

PD1_sub5_filt2 <- PD1_sub4_filt2[complete.cases(PD1_sub4_filt2),]

#imputing on NAs with row 1/2 min 

replace_na <- function(df){
  # df<- PD1_sub3_filt2_22imp
  for(unique_row in 1:nrow(df)){
    # unique_row <- 1
    sub_df <- df[unique_row,]
    df[unique_row, is.na(sub_df)] <- sub_df$rowmin
  }
  return(df)
}

#removal of irs channels, imputing and adding back in of irs
PD1_sub5A_filt_22nas <- PD1_sub4_filt_22nas[c(2:10,12:19,21:26,29:34,37:43,45,46)]

PD1_sub5B_filt_22nas <- PD1_sub4_filt_22nas[c(1,11,20,27,28,35,36,44,47:60)]

PD1_sub5A_filt_22imp <- replace_na(PD1_sub5A_filt_22nas)

PD1_sub5_filt_22imp <- cbind(PD1_sub5A_filt_22imp, PD1_sub5B_filt_22nas)


PD1_sub5A_filt2_22nas <- PD1_sub4_filt2_22nas[c(2:10,12:19,21:25,28:33,36:42,44,45)]

PD1_sub5B_filt2_22nas <- PD1_sub4_filt2_22nas[c(1,11,20,26,27,34,35,43,46:59)]

PD1_sub5A_filt2_22imp <- replace_na(PD1_sub5A_filt2_22nas)

PD1_sub5_filt2_22imp <- cbind(PD1_sub5A_filt2_22imp, PD1_sub5B_filt2_22nas)


## check all rows are complete

PD1_sub6_filt_22imp <- PD1_sub5_filt_22imp[complete.cases(PD1_sub5_filt_22imp),]

PD1_sub6_filt2_22imp <- PD1_sub5_filt2_22imp[complete.cases(PD1_sub5_filt2_22imp),]


write.csv(PD1_sub6_filt_22imp, file = "PD1_sub6_filt_22imp.csv")

write.csv(PD1_sub6_filt2_22imp, file = "PD1_sub6_filt2_22imp.csv")

write.csv(PD1_sub5, file = "PD1_sub5.csv")

write.csv(PD1_sub5_filt, file = "PD1_sub5_filt.csv")

write.csv(PD1_sub5_filt2, file = "PD1_sub5_filt2.csv")

# subset columns for normalization

PD1_sub6 <- PD1_sub5[c(2:51)]

PD1_sub6_filt <- PD1_sub5_filt[c(2:45)]

PD1_sub6_filt2 <- PD1_sub5_filt2[c(2:44)]

PD1_sub7_filt_22imp <- PD1_sub6_filt_22imp[c(1:9,40,10:17,41,18:23,42,43,24:29,44,45,30:36,46,37)]

PD1_sub7_filt2_22imp <- PD1_sub6_filt2_22imp[c(1:9,39,10:17,40,18:22,41,42,23:28,43,44,29:35,45,36)]

# sample loading normalization
sl_normalization <- function(protein_df, tmt_exp_columns){
  
  # protein_df <- pd_light_all
  # tmt_exp_columns <- list(c(4:11), c(1:3))
  
  expect_is(tmt_exp_columns, 'list')
  expect_is(protein_df, 'data.frame')
  
  number_tmt_exps <- length(tmt_exp_columns)
  
  separated_prot_vals <- list()
  
  for(i in 1:number_tmt_exps){
    separated_prot_vals[[i]] <- protein_df[,tmt_exp_columns[[i]]]
  }
  
  # exp1_vals <- protein_df[,tmt1]
  # exp2_vals <- protein_df[,tmt2]
  
  normalization_factor_list <- list()
  
  for(i in 1:number_tmt_exps){
    normalization_factor_list[[i]] <- mean(colSums(separated_prot_vals[[i]])) / colSums(separated_prot_vals[[i]])
  }
  
  sl_normalized_prot_list <- list()
  
  for(i in 1:number_tmt_exps){
    sl_normalized_prot_list[[i]] <- sweep(separated_prot_vals[[i]], 2, normalization_factor_list[[i]], FUN = "*")
  }
  
  # col_blank <- rep('firebrick', length(box_labels))
  # col_blank[which(grepl(pattern = "_2$", x = box_labels))] <- 'darkblue'
  
  plot_df <- do.call('cbind', sl_normalized_prot_list)
  
  par(mar = c(7, 5, 3, 3))
  boxplot(log2(plot_df),
          # col = col_blank,
          xaxt = 'n',
          xlab = '',
          main = 'SL Normalization',
          ylab = 'Intensity')
  axis(1,
       labels = FALSE)
  text(x = seq_along(names(plot_df)),
       y = par("usr")[3] - 0.5,
       srt = 40,
       adj = 0.75,
       cex = 0.8,
       labels = names(plot_df),
       xpd = TRUE)
  
  return(sl_normalized_prot_list)
  
}

#test <- sl_normalization(protein_df = PD1_sub4, tmt_exp_columns = list(c(1:10), c(11:19), c(20:27), c(28:35), c(36:44)))


#write.csv(test, file = "PD1_sub4.csv")

#PD1_sub4 <- read.csv('PD1_sub3_filt1_Norm.csv')
#PD1_sub4 <- PD1_sub4[c(2:48)]


# Internal reference standard normalization
irs_normalization <- function(sl_normalized_list, tmt_common_channel_names){
  
  # sl_normalized_list <- test
  # tmt_common_channel_names <- list(c("hl_highmn_highfe_1_1", "hl_highmn_highfe_2_1"), c("hl_highmn_highfe_1_2", "hl_highmn_highfe_2_2"))
  # tmt1_common_channel <- c("hl_highmn_highfe_1_1", "hl_highmn_highfe_2_1")
  # tmt2_common_channel <- c("hl_highmn_highfe_1_2", "hl_highmn_highfe_2_2")
  
  number_tmt_exps <- length(sl_normalized_list)
  
  irs_list <- list()
  irs_rowmeans_matrix <- matrix(nrow = length(rowMeans(sl_normalized_list[[1]])))
  
  for(i in 1:number_tmt_exps){
    # i <- 1
    sl_normalized_df <- sl_normalized_list[[i]]
    irs_list[[i]] <- sl_normalized_df[tmt_common_channel_names[[i]]]
    irs_rowmeans_matrix <- cbind(irs_rowmeans_matrix, rowMeans(irs_list[[i]]))
  }
  
  rowmean_vector <- rowMeans(irs_rowmeans_matrix, na.rm = TRUE)
  
  scaling_factor_list <- list()
  irs_sl_normalized_list <- list()
  
  for(i in 1:number_tmt_exps){
    scaling_factor_list[[i]] <- rowmean_vector / rowMeans(irs_list[[i]])
    sl_normalized_df <- sl_normalized_list[[i]]
    irs_sl_normalized_list[[i]] <- sl_normalized_df * scaling_factor_list[[i]]
  }
  
  irs_sl_df <- matrix(nrow = nrow(irs_sl_normalized_list[[1]])) %>% as.data.frame()
  
  for(i in 1:number_tmt_exps){
    irs_sl_df <- cbind(irs_sl_df, irs_sl_normalized_list[[i]])
  }
  
  return_df <- irs_sl_df[,-1]
  # 
  # col_blank <- rep('firebrick', length(box_labels))
  # col_blank[which(grepl(pattern = "_2$", x = box_labels))] <- 'darkblue'
  # 
  par(mar = c(7, 5, 3, 3))
  boxplot(log2(return_df),
          # col = col_blank,
          xaxt = 'n',
          xlab = '',
          main = 'SL, IRS Normalization',
          ylab = 'Intensity')
  axis(1,
       labels = FALSE)
  text(x = seq_along(names(return_df)),
       y = par("usr")[3] - 0.5,
       srt = 40,
       adj = 0.75,
       cex = 0.8,
       labels = names(return_df),
       xpd = TRUE)
  
  return(return_df)
  
}

# edgeR tmm normalization
tmm_normalization <- function(irs_normalized_df){
  
  irs_tmm <- edgeR::calcNormFactors(irs_normalized_df)
  data_irs_tmm <- sweep(irs_normalized_df, 2, irs_tmm, FUN = "/")
  
  
  boxplot(log2(data_irs_tmm), 
          # col = c(rep('firebrick', length(tmt1)),
          # rep('darkblue', length(tmt2))),
          # col = col_blank,
          xaxt = 'n', 
          xlab = '',
          ylab = 'Intensity',
          main = 'SL, IRS, TMM Normalization')
  axis(1, 
       labels = FALSE)
  text(x = seq_along(names(data_irs_tmm)), 
       y = par("usr")[3] - 0.5, 
       srt = 45, 
       adj = 1,
       labels = names(data_irs_tmm), 
       xpd = TRUE)
  
  return(data_irs_tmm)
}

sl_irs_tmm_normalization <- function(protein_df, tmt_exp_columns, tmt_common_channel_names){
  
  # check that the inputs are of the correct type
  expect_is(protein_df, 'data.frame')
  expect_is(tmt_exp_columns, 'list')
  expect_is(tmt_common_channel_names, 'list')
  
  # checking input files
  is.contained <- function(x, y) {
    z <- x[x %in% setdiff(x, y)]
    length(z) == length(x) - length(y)
  }
  
  # checking that the common channel names are within the protein df
  testing_tmt_names <- is.contained(x = names(protein_df), y = unlist(tmt_common_channel_names))
  
  if(!testing_tmt_names){
    stop('It seems like the tmt_common_channel_names contains a name that is not in your protein_df; double check that!')
  }
  
  complete_protein_df <- protein_df[complete.cases(protein_df),]
  
  if(nrow(complete_protein_df) != nrow(protein_df)){
    warning('Your protein_df has some missing values. This normalization can only consider proteins observed across all TMT channels. The output of this function subsets only the normalized proteins observed across all channels!')
    protein_df <- complete_protein_df
  }
  
  for(i in 1:length(tmt_exp_columns)){
    if(i == 1) print('These are the tmt experiment column labels youve designated, just to be sure:')
    print('---------------------------------')
    print(paste0('TMT Experiment ', i))
    names(protein_df)[tmt_exp_columns[[i]]] %>% print()
    print('---------------------------------')
  }
  
  par(mfrow = c(2, 2))
  
  boxplot(log2(protein_df), 
          # col = c(rep('firebrick', length(tmt1)),
          # rep('darkblue', length(tmt2))),
          # col = col_blank,
          xaxt = 'n', 
          xlab = '',
          ylab = 'Intensity',
          main = 'No Normalization')
  axis(1, 
       labels = FALSE)
  text(x = seq_along(names(protein_df)), 
       y = par("usr")[3] - 0.5, 
       srt = 45, 
       adj = 1,
       labels = names(protein_df), 
       xpd = TRUE)
  
  sl_norm_file <- sl_normalization(protein_df = protein_df, tmt_exp_columns = tmt_exp_columns)
  irs_norm_file <- irs_normalization(sl_normalized_list = sl_norm_file, tmt_common_channel_names = tmt_common_channel_names)
  tmm_norm_file <- tmm_normalization(irs_normalized_df = irs_norm_file)
  return(tmm_norm_file)
}



#### Running the normalization:


PD1_Norm <- sl_irs_tmm_normalization(protein_df = PD1_sub6, 
                                           tmt_exp_columns = list(c(1:10), c(11:20), c(21:30), c(31:40), c(41:50)), 
                                           tmt_common_channel_names = list(c("A10_pool2"), c("B10_pool2"), c("C9_pool1", "C10_pool2"), c("D9_pool1", "D10_pool2"), c("E9_pool1")))


PD1_Norm_filt <- sl_irs_tmm_normalization(protein_df = PD1_sub6_filt, 
                                 tmt_exp_columns = list(c(1:10), c(11:19), c(20:27), c(28:35), c(36:44)), 
                                 tmt_common_channel_names = list(c("A10_pool2"), c("B10_pool2"), c("C9_pool1", "C10_pool2"), c("D9_pool1", "D10_pool2"), c("E9_pool1")))


PD1_Norm_filt2 <- sl_irs_tmm_normalization(protein_df = PD1_sub6_filt2, 
                                           tmt_exp_columns = list(c(1:10), c(11:19), c(20:26), c(27:34), c(35:43)), 
                                           tmt_common_channel_names = list(c("A10_pool2"), c("B10_pool2"), c("C9_pool1", "C10_pool2"), c("D9_pool1", "D10_pool2"), c("E9_pool1")))


PD1_Norm_filt_22imp <- sl_irs_tmm_normalization(protein_df = PD1_sub7_filt_22imp, 
                                                 tmt_exp_columns = list(c(1:10), c(11:19), c(20:27), c(28:35), c(36:44)), 
                                                 tmt_common_channel_names = list(c("A10_pool2"), c("B10_pool2"), c("C9_pool1", "C10_pool2"), c("D9_pool1", "D10_pool2"), c("E9_pool1")))

PD1_Norm_filt2_22imp <- sl_irs_tmm_normalization(protein_df = PD1_sub7_filt2_22imp, 
                                                tmt_exp_columns = list(c(1:10), c(11:19), c(20:26), c(27:34), c(35:43)), 
                                                tmt_common_channel_names = list(c("A10_pool2"), c("B10_pool2"), c("C9_pool1", "C10_pool2"), c("D9_pool1", "D10_pool2"), c("E9_pool1")))


#add back in protein accessions

PD1_sub7 <- PD1_sub5[c(1,53:66)]

PD1_sub7_filt <- PD1_sub5_filt[c(1,47:60)]

PD1_sub7_filt2 <- PD1_sub5_filt2[c(1,46:59)]

PD1_sub8_filt_22imp <- PD1_sub6_filt_22imp[c(38,47:60)]

PD1_sub8_filt2_22imp <- PD1_sub6_filt2_22imp[c(38,46:59)]


PD1_Norm_final <- cbind(PD1_Norm, PD1_sub7)

PD1_Norm_filt_final <- cbind(PD1_Norm_filt, PD1_sub7_filt)

PD1_Norm_filt2_final <- cbind(PD1_Norm_filt2, PD1_sub7_filt2)

PD1_Norm_impute_final <- cbind(PD1_Norm_filt_22imp, PD1_sub8_filt_22imp)

PD1_Norm_impute2_final <- cbind(PD1_Norm_filt2_22imp, PD1_sub8_filt2_22imp)


#output normalized data

write.csv(PD1_Norm_final, file = "PD1_Norm.csv")

write.csv(PD1_Norm_filt_final, file = "PD1_Norm_filt.csv")

write.csv(PD1_Norm_filt2_final, file = "PD1_Norm_filt2.csv")

write.csv(PD1_Norm_impute_final, file = "PD1_Norm_impute.csv")

write.csv(PD1_Norm_impute2_final, file = "PD1_Norm_impute2.csv")

