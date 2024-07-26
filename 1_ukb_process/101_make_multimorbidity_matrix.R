# Generate ICD10 blocks multimorbidity matrix

# Load libraries

print("Loading libraries...")

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(data.table)
library(magrittr)

# Variables

# ICD10 blocks dictionary
icd10_blocks_file = "icd10_blocks.tsv"

# table containing UK Biobank ICD10 diagnoses table
ukb_icd10_diag_file = ""

# column names of ICD10 diagnoses columns
icd10_cols = sapply(0:212, function(x) paste0("diagnoses_icd10_f41270_0_", x))

# table containing sex of UK Biobank participants; coded as '0' for females and '1' for males.
ukb_sex_file = ""

# column name of sex column
sex_col <- ""

# choose level of analysis: block or chapter
icd10_type = "block"

# output files

# multimorbidity matrix
multimorbidity_file = "ukb_multimorbidity_ICD10_block_table.tsv"

# female multimorbidity matrix
fem_mm_file = "ukb_fem_multimorbidity_ICD10_block_table.tsv"

# male multimorbidity matrix
male_mm_file = "ukb_male_multimorbidity_ICD10_block_table.tsv"

# Main Functions

make_icd10_block <- function(lett, first, second) {
  nums <- first:second
  single <- nums < 10
  nums <- as.character(nums)
  nums[single] <- paste0("0", nums[single])
  lett_num <- paste0(lett, nums)
  return(lett_num)
}

make_icd10_block2 <- function(parsed_vec) {
  lett <- parsed_vec[1]
  first <- as.numeric(parsed_vec[2])
  second <- as.numeric(parsed_vec[3])
  icd10_block <- make_icd10_block(lett = lett, first = first, second = second)
  return(icd10_block)
}

icd10_pad <- function(num) {
  return(str_pad(num, width = 2, pad = "0"))
}

chapters <- list(c(paste0("A", icd10_pad(0:99)), icd10_pad(paste0("B", 0:99))),
                 c(paste0("C", icd10_pad(0:99)), paste0("D", icd10_pad(0:49))),
                 c(paste0("D", icd10_pad(50:89))),
                 c(paste0("E", icd10_pad(0:89))),
                 c(paste0("F", icd10_pad(1:99))),
                 c(paste0("G", icd10_pad(0:99))),
                 c(paste0("H", icd10_pad(0:59))),
                 c(paste0("H", icd10_pad(60:95))),
                 c(paste0("I", icd10_pad(0:99))),
                 c(paste0("J", icd10_pad(0:99))),
                 c(paste0("K", icd10_pad(0:95))),
                 c(paste0("L", icd10_pad(0:99))),
                 c(paste0("M", icd10_pad(0:99))),
                 c(paste0("N", icd10_pad(0:99))))

# Read in ICD10

ukb_data_icd10 <- fread(ukb_icd10_diag_file)

# Filter for ICD10 all col

ukb_data_icd10 %<>% select(icd10_cols)

# Read in ICD10 block file

icd10_blocks_df <- read_tsv(icd10_blocks_file)

# get icd10 block codes

if(icd10_type == "block") {
  
  # Convert ICD10 codes for participants to 3 character codes
  
  ukb_data_icd10 <- as.matrix(ukb_data_icd10)
  
  ukb_data_icd10 <- apply(ukb_data_icd10, 2, function(x) str_trunc(x, 3, ellipsis = ""))
  
  # convert to ICD10 block data
  
  ukb_data_icd10 <- apply(ukb_data_icd10, 2, function(x) icd10_blocks_df$block[match(x, icd10_blocks_df$code)])
  
  # Create multimorbidity matrix
  
  multimorbidity_matrix <- matrix(nrow = nrow(ukb_data_icd10), ncol = length(unique(icd10_blocks_df$block)))
  
  colnames(multimorbidity_matrix) <- unique(icd10_blocks_df$block)

  # Loop through ICD10 codes and assign blocks to multimorbidity block matrix
  
  for(i in 1:nrow(ukb_data_icd10)) {
    participant_i <- ukb_data_icd10[i,]
    add_in <- which(colnames(multimorbidity_matrix) %in% participant_i)
    if(!is_empty(add_in)) {
      multimorbidity_matrix[i,add_in] <- 1
    } else {
      next
    }
    if(i %% 1000 == 0) {
      print(i)
    }
  }
  
  # write out block multimorbidity matrix
  
  write_tsv(as.data.frame(multimorbidity_matrix), multimorbidity_file)
  
  # sex filter
  
  # read in sex file
  
  ukb_sex <- as.data.frame(fread(ukb_sex_file))
  
  # get sex vector
  
  ukb_sex %<>% pull(sex_col)
  
  # get sex indices
  
  fem_ind <- which(ukb_sex == 0)
  male_ind <- which(ukb_sex == 1)
  
  # filter multimorbidity matrix
  
  fem_mm <- multimorbidity_matrix[fem_ind,]
  male_mm <- multimorbidity_matrix[male_ind,]
  
  # write out sex multimorbidity matrices
  
  write_tsv(as.data.frame(fem_mm), fem_mm_file)
  
  write_tsv(as.data.frame(male_mm), male_mm_file)
  
} else if(icd10_type == "chapter") {
  
  icd10_chapt_codes <- unique(icd10_blocks_df$chapter)[1:14] # only physiological chapters
  
  # Convert ICD10 codes for participants to 3 character codes
  
  ukb_data_icd10 <- as.matrix(ukb_data_icd10)
  
  ukb_data_icd10 <- apply(ukb_data_icd10, 2, function(x) str_trunc(x, 3, ellipsis = ""))
  
  # Create multimorbidity matrix
  
  multimorbidity_matrix <- matrix(nrow = nrow(ukb_data_icd10), ncol = length(icd10_chapt_codes))
  
  colnames(multimorbidity_matrix) <- icd10_chapt_codes
  
  # Loop through ICD10 codes and assign blocks to multimorbidity block matrix
  
  for(j in 1:length(icd10_chapt_codes)) {
    chapt <- chapters[[j]]
    for(i in 1:nrow(ukb_data_icd10)) {
      participant_i <- ukb_data_icd10[i,]
      matrix_add_logi <- any(participant_i %in% chapt)
      if(matrix_add_logi) {
        multimorbidity_matrix[i,j] <- 1
      } else {
        multimorbidity_matrix[i,j] <- 0
      }
    }
    print(j)
  }
  
  # write out block multimorbidity matrix
  
  write_tsv(as.data.frame(multimorbidity_matrix), multimorbidity_file)
  
  # sex filter
  
  # read in sex file
  
  ukb_sex <- as.data.frame(fread(ukb_sex_file))
  
  # get sex vector
  
  ukb_sex %<>% pull(sex_col)
  
  # get sex indices
  
  fem_ind <- which(ukb_sex == 0)
  male_ind <- which(ukb_sex == 1)
  
  # filter multimorbidity matrix
  
  fem_mm <- multimorbidity_matrix[fem_ind,]
  male_mm <- multimorbidity_matrix[male_ind,]
  
  # write out sex multimorbidity matrices
  
  write_tsv(as.data.frame(fem_mm), fem_mm_file)
  
  write_tsv(as.data.frame(male_mm), male_mm_file)
  
}




