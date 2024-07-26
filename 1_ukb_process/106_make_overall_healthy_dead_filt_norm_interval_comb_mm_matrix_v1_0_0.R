# create time filtered multimorbidity matrices based on time of diagnoses from first diagnoses
rm(list = ls())

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Variables

# ICD10 blocks dictionary
icd10_blocks_file = "icd10_blocks.tsv"

# diagnosis trajectory matrix
traj_file = "fem_disease_block_1pct_traj.tsv"

# diagnosis time matrix
traj_dates_file = "fem_diag_block_1pct_norm_dates.tsv"

# death time matrix
death_interval_file = "fem_norm_death.tsv"

# follow-up time

fu_file = "fem_norm_follow_up.tsv"

# diagnosis time distribution quartiles
interval_file = "fem_1pct_block_quartile_interval_norm_dates.tsv"

# output

# population frequencies of filtered multimorbidity matrices
seq_mm_sub_cnts_file = "fem_disease_block_1pct_traj_comb_mm_sub_cnts.tsv"

# total combined multimorbidity matrix of 4 quartile-filtered multimorbidity matrices
combined_mm_file = "fem_disease_block_1pct_traj_comb_mm.tsv"

# 4 uncombined quartile-filtered multimorbidity matrices
uncombined_mm_files = paste0("fem_disease_block_1pct_traj_", 1:4, "_mm.tsv")

# functions

make_temporal <- function(interval) {
  int_lt <- matrix(nrow = nrow(dates), ncol = ncol(dates))
  for(i in 1:nrow(int_lt)) {
    indiv_i <- dates[i,]
    int_lt_i <- indiv_i <= interval
    int_lt[i,] <- int_lt_i
  }
  
  interval_diag <- traj
  
  for(i in 1:nrow(int_lt)) {
    row_i <- int_lt[i,]
    int_row_i <- which(row_i)
    if(length(int_row_i) == 0 | length(int_row_i) == ncol(interval_diag)) {
      next
    }
    interval_diag[i, -int_row_i] <- NA
  }
  
  interval_mm <- matrix(nrow = nrow(traj), ncol = length(dis_uni))
  colnames(interval_mm) <- dis_uni
  interval_diag <- as.matrix(interval_diag)
  
  for(i in 1:nrow(interval_diag)) {
    participant <- interval_diag[i,]
    block_ind <- which(colnames(interval_mm) %in% participant)
    if(length(block_ind) == 0) {
      interval_mm[i,] <- 0
    } else {
      interval_mm[i, block_ind] <- 1
      interval_mm[i, -block_ind] <- 0
    }
  }
  return(interval_mm)
}

# read in data

traj <- read_tsv(traj_file)
traj <- as.matrix(traj)

dates <- read_tsv(traj_dates_file)
dates <- as.matrix(dates)

deaths <- read_tsv(death_interval_file)
deaths <- pull(deaths, 1)

fu <- fread(fu_file)
fu <- fu %>% pull(1)

interval <- read_tsv(interval_file)
interval <- pull(interval, 1)

dis_uni <- sort(unique(as.vector(traj)))

int_1_mm <- make_temporal(interval[1])
int_2_mm <- make_temporal(interval[2])
int_3_mm <- make_temporal(interval[3])
int_4_mm <- make_temporal(interval[4])

# remove those who only have 1 disease by the end

rm_ind <- which(rowSums(int_4_mm) <= 1)

deaths <- deaths[-rm_ind]
fu <- fu[-rm_ind]

int_1_mm_filt <- int_1_mm[-rm_ind,]
int_1_h <- rowSums(int_1_mm_filt) == 0
int_1_d <- deaths <= interval[1]
int_1_fu <- fu <= interval[1]
int_1_hd <- which(int_1_h | int_1_d | int_1_fu)
int_1_mm_filt <- int_1_mm_filt[-int_1_hd,]

int_2_mm_filt <- int_2_mm[-rm_ind,]
int_2_h <- rowSums(int_2_mm_filt) == 0
int_2_d <- deaths <= interval[2]
int_2_fu <- fu <= interval[2]
int_2_hd <- which(int_2_h | int_2_d | int_2_fu)
int_2_mm_filt <- int_2_mm_filt[-int_2_hd,]

int_3_mm_filt <- int_3_mm[-rm_ind,]
int_3_h <- rowSums(int_3_mm_filt) == 0
int_3_d <- deaths <= interval[3]
int_3_fu <- fu <= interval[3]
int_3_hd <- which(int_3_h | int_3_d | int_3_fu)
int_3_mm_filt <- int_3_mm_filt[-int_3_hd,]

int_4_mm_filt <- int_4_mm[-rm_ind,]
int_4_h <- rowSums(int_4_mm_filt) == 0
int_4_d <- deaths <= interval[4]
int_4_fu <- fu <= interval[4]
int_4_hd <- which(int_4_h | int_4_d | int_4_fu)
int_4_mm_filt <- int_4_mm_filt[-int_4_hd,]


# combine mm

row_cnts <- c(nrow(int_1_mm_filt),
              nrow(int_2_mm_filt),
              nrow(int_3_mm_filt),
              nrow(int_4_mm_filt))

# write out row counts

write_tsv(as.data.frame(row_cnts), seq_mm_sub_cnts_file)

# write out combined mm

int_comb_mm <- rbind(int_1_mm_filt,
                     int_2_mm_filt,
                     int_3_mm_filt,
                     int_4_mm_filt)

write_tsv(as.data.frame(int_comb_mm), combined_mm_file)

# write out individual mm's

int_uncomb_mm <- list(int_1_mm_filt,
                      int_2_mm_filt,
                      int_3_mm_filt,
                      int_4_mm_filt)

for(i in 1:length(int_uncomb_mm)) {
  write_tsv(as.data.frame(int_uncomb_mm[[i]]), uncombined_mm_files[i])
}

