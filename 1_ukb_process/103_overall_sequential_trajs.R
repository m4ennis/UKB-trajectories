# create ordered trajectories of diseases 

rm(list = ls())

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(data.table)
library(magrittr)

# Variables

# input files

# ICD10 blocks dictionary
icd10_blocks_file = "icd10_blocks.tsv"

# file containing ICD10 diagnoses
ukb_icd10_diag_file = ""

# file containing ICD10 diagnosis dates
ukb_icd10_diag_date_file = ""

# file containing sex of UK Biobank participants
ukb_sex_file = ""

# column names of ICD10 diagnoses columns
icd10_cols = sapply(0:211, function(x) paste0("diagnoses_icd10_f41270_0_", x))

# column names of ICD10 dates columns
icd10_dates_cols = sapply(0:211, function(x) paste0("date_of_first_inpatient_diagnosis_icd10_f41280_0_", x))

# column name of sex column
sex_col <- ""

# file containing female diseases passing prevalence filter
fem_disease_file = ""

# file containing male diseases passing prevalence filter
male_disease_file = ""

# output files

# female diagnosis trajectories
fem_traj_file = "fem_disease_block_1pct_traj.tsv"
# female diagnosis dates
fem_traj_dates_file = "fem_diag_block_1pct_dates.tsv"

# male diagnosis trajectories
male_traj_file = "male_disease_block_1pct_traj.tsv"
# male diagnosis dates
male_traj_dates_file = "male_diag_block_1pct_dates.tsv"

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

# Load data

ukb_data_icd10 <- as.data.frame(fread(ukb_icd10_diag_file))

# read in temporal data

temporal <- as.data.frame(fread(ukb_icd10_diag_date_file))

# do eids match?

if(any(!temporal$eid == ukb_data_icd10$eid)) {
  stop("eIDs do not match up")
}

# Filter for ICD10 all col

ukb_data_icd10 %<>% select(all_of(icd10_cols))

# diagnosis dates

dates <- temporal %>% select(all_of(icd10_dates_cols))

dates <- as.matrix(dates)

# read in sex file

ukb_sex <- as.data.frame(fread(ukb_sex_file))

# get sex vector

ukb_sex %<>% pull(sex_col)

# get sex indices

fem_ind <- which(ukb_sex == 0)
male_ind <- which(ukb_sex == 1)

# read in icd10 blocks

blocks <- read_tsv(icd10_blocks_file)

block_codes <- blocks$block
block_names <- blocks$name

# read in female/male disease use

fem_disease <- pull(read_tsv(fem_disease_file), 1)

male_disease <- pull(read_tsv(male_disease_file), 1)

# order codes by date of diagnosis

icd10_order <- as.matrix(ukb_data_icd10)
colnames(icd10_order) <- NULL

dates_order <- dates
colnames(dates_order) <- NULL

for(i in 1:nrow(icd10_order)) {
  icd10_i <- icd10_order[i,]
  order_i <- order(dates[i,])
  dates_order[i,] <- dates_order[i, order_i]
  icd10_order_i <- icd10_i[order_i]
  icd10_order[i,] <- icd10_order_i
}

# parse block codes

fem_parsed_block <- lapply(fem_disease, function(x) {
  lett <- substr(x, 1, 1)
  first_num <- substr(x, 2, 3)
  second_num <- substr(x, 6, 7)
  return(c(lett, first_num, second_num))
})

fem_blocks <- lapply(fem_parsed_block, function(x) make_icd10_block2(x))
names(fem_blocks) <- fem_disease

male_parsed_block <- lapply(male_disease, function(x) {
  lett <- substr(x, 1, 1)
  first_num <- substr(x, 2, 3)
  second_num <- substr(x, 6, 7)
  return(c(lett, first_num, second_num))
})

male_blocks <- lapply(male_parsed_block, function(x) make_icd10_block2(x))
names(male_blocks) <- male_disease

# create code-block map

fem_block_lens <- lengths(fem_blocks)
fem_blocks_rep <- unlist(map2(fem_disease, fem_block_lens, function(x, y) rep(x, y)))
fem_blocks_map <- cbind(unlist(fem_blocks), fem_blocks_rep)
colnames(fem_blocks_map) <- c("code", "block")

male_block_lens <- lengths(male_blocks)
male_blocks_rep <- unlist(map2(male_disease, male_block_lens, function(x, y) rep(x, y)))
male_blocks_map <- cbind(unlist(male_blocks), male_blocks_rep)
colnames(male_blocks_map) <- c("code", "block")

# Convert ICD10 codes for participants to 3 character codes

icd10_order <- apply(icd10_order, 2, function(x) str_trunc(x, 3, ellipsis = ""))

# filter for sex

fem_icd10_order <- icd10_order[fem_ind,]
fem_dates_order <- dates_order[fem_ind,]

male_icd10_order <- icd10_order[male_ind,]
male_dates_order <- dates_order[male_ind,]

# convert to blocks

fem_icd10_order2 <- fem_icd10_order
fem_dates_order2 <- fem_dates_order

for(i in 1:nrow(fem_icd10_order2)) {
  participant <- fem_icd10_order2[i,]
  block_ind <- match(participant, fem_blocks_map[,1])
  block_i <- fem_blocks_map[block_ind, 2]
  fem_icd10_order2[i, ] <- block_i
  fem_dates_order2[i, which(is.na(block_i))] <- NA
}

male_icd10_order2 <- male_icd10_order
male_dates_order2 <- male_dates_order

for(i in 1:nrow(male_icd10_order2)) {
  participant <- male_icd10_order2[i,]
  block_ind <- match(participant, male_blocks_map[,1])
  block_i <- male_blocks_map[block_ind, 2]
  male_icd10_order2[i, ] <- block_i
  male_dates_order2[i, which(is.na(block_i))] <- NA
}

# collapse NAs

fem_icd10_order3 <- fem_icd10_order2
fem_dates_order3 <- fem_dates_order2

for(i in 1:nrow(fem_icd10_order3)) {
  indiv_i <- fem_icd10_order3[i,]
  dates_i <- fem_dates_order3[i,]
  na_ind <- is.na(indiv_i)
  no_na <- indiv_i[!na_ind]
  dates_no_na <- dates_i[!na_ind]
  clpse_i <- c(no_na, rep(NA, ncol(fem_icd10_order3) - length(no_na)))
  clpse_date_i <- c(dates_no_na, rep(NA, ncol(fem_dates_order3) - length(dates_no_na)))
  fem_icd10_order3[i,] <- clpse_i
  fem_dates_order3[i,] <- clpse_date_i
}

male_icd10_order3 <- male_icd10_order2
male_dates_order3 <- male_dates_order2

for(i in 1:nrow(male_icd10_order3)) {
  indiv_i <- male_icd10_order3[i,]
  dates_i <- male_dates_order3[i,]
  na_ind <- is.na(indiv_i)
  no_na <- indiv_i[!na_ind]
  dates_no_na <- dates_i[!na_ind]
  clpse_i <- c(no_na, rep(NA, ncol(male_icd10_order3) - length(no_na)))
  clpse_date_i <- c(dates_no_na, rep(NA, ncol(male_dates_order3) - length(dates_no_na)))
  male_icd10_order3[i,] <- clpse_i
  male_dates_order3[i,] <- clpse_date_i
}

# remove duplicates

for(i in 1:nrow(fem_icd10_order3)) {
  indiv_i <- fem_icd10_order3[i,]
  dates_i <- fem_dates_order3[i,]
  dup_ind <- duplicated(indiv_i)
  no_dup <- indiv_i[!dup_ind]
  dates_no_dup <- dates_i[!dup_ind]
  no_dup_o <- c(no_dup, rep(NA, ncol(fem_icd10_order3) - length(no_dup)))
  date_no_dup_o <- c(dates_no_dup, rep(NA, ncol(fem_dates_order3) - length(dates_no_dup)))
  fem_icd10_order3[i,] <- no_dup_o
  fem_dates_order3[i,] <- date_no_dup_o
}

for(i in 1:nrow(male_icd10_order3)) {
  indiv_i <- male_icd10_order3[i,]
  dates_i <- male_dates_order3[i,]
  dup_ind <- duplicated(indiv_i)
  no_dup <- indiv_i[!dup_ind]
  dates_no_dup <- dates_i[!dup_ind]
  no_dup_o <- c(no_dup, rep(NA, ncol(male_icd10_order3) - length(no_dup)))
  date_no_dup_o <- c(dates_no_dup, rep(NA, ncol(male_dates_order3) - length(dates_no_dup)))
  male_icd10_order3[i,] <- no_dup_o
  male_dates_order3[i,] <- date_no_dup_o
}


# remove columns with all NA

fem_na_cols <- which(apply(fem_icd10_order3, 2, function(x) sum(is.na(x))) == nrow(fem_icd10_order3))
fem_icd10_order3 <- fem_icd10_order3[, -fem_na_cols]

male_na_cols <- which(apply(male_icd10_order3, 2, function(x) sum(is.na(x))) == nrow(male_icd10_order3))
male_icd10_order3 <- male_icd10_order3[, -male_na_cols]

fem_na_cols2 <- which(apply(fem_dates_order3, 2, function(x) sum(is.na(x))) == nrow(fem_dates_order3))
fem_dates_order3 <- fem_dates_order3[, -fem_na_cols2]

male_na_cols2 <- which(apply(male_dates_order3, 2, function(x) sum(is.na(x))) == nrow(male_dates_order3))
male_dates_order3 <- male_dates_order3[, -male_na_cols2]

# write out trajectories

write_tsv(as.data.frame(fem_icd10_order3), fem_traj_file)

write_tsv(as.data.frame(male_icd10_order3), male_traj_file)

# write out trajectory dates

fem_dates_order3 <- as.data.frame(fem_dates_order3)

write.table(fem_dates_order3, fem_traj_dates_file, sep = "\t", row.names = F, col.names = F)

male_dates_order3 <- as.data.frame(male_dates_order3)

write.table(male_dates_order3, male_traj_dates_file, sep = "\t", row.names = F, col.names = F)