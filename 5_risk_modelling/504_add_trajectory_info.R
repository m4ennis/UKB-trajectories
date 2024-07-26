# add information such as population frequency and time length of trajectories

library(tidyverse)
library(data.table)
library(magrittr)
library(parallel)
library(foreach)
library(doParallel)

# variables

# female diagnosis trajectory matrix
fem_traj_file = "fem_disease_block_1pct_traj.tsv"
# male diagnosis trajectory matrix
male_traj_file = "male_disease_block_1pct_traj.tsv"

# female time to diagnosis matrix
fem_dates_file = "fem_diag_block_1pct_norm_dates.tsv"
# male time to diagnosis matrix
male_dates_file = "male_diag_block_1pct_norm_dates.tsv"

# female dates of diagnosis matrix
fem_unnorm_dates_file = "fem_diag_block_1pct_dates.tsv"
# male dates of diagnosis matrix
male_unnorm_dates_file = "male_diag_block_1pct_dates.tsv"

# female significant built trajectories
fem_sig_file = "female_significant_trajectory_info.tsv"
# male significant built trajectories
male_sig_file = "male_significant_trajectory_info.tsv"

# ukb file containing age at attending UK Biobank assessment centre and participant sex
ukb_file = "C:/Users/matth/OneDrive - Ulster University/paper_1/analysis/data/ukb_raw/ukb_df_general_col_filtered_all.tsv"
# ukb file containing date of attending assessment centre
ukb_file2 = "C:/Users/matth/OneDrive - Ulster University/paper_1/analysis/data/ukb_raw/ukb_df_t2d_date_attendence_filtered.tsv"

min_comp = 20 # minimum number of individuals required to have a trajectory
n_days = 14 # minimum number of days between diagnoses to be considered a trajectory

# read in trajs

fem_traj <- fread(fem_traj_file)
male_traj <- fread(male_traj_file)

# read in traj times

fem_dates <- fread(fem_dates_file)
male_dates <- fread(male_dates_file)

# read in trajectory dates

fem_unnorm_dates <- fread(fem_unnorm_dates_file)
male_unnorm_dates <- fread(male_unnorm_dates_file)

# read in UKB file

ukb <- fread(ukb_file)

age <- ukb$age_when_attended_assessment_centre_f21003_0_0

fem_age <- age[ukb$sex_f31_0_0 == 0]
male_age <- age[ukb$sex_f31_0_0 == 1]

# read in UKB file 2

ukb2 <- fread(ukb_file2)

bl_date <- ukb2$date_of_attending_assessment_centre_f53_0_0
fem_bl <- bl_date[ukb$sex_f31_0_0 == 0]
male_bl <- bl_date[ukb$sex_f31_0_0 == 1]

# read in significant trajectories

fem_sig <- fread(fem_sig_file)
male_sig <- fread(male_sig_file)

# get maximum length of trajectories

fem_max <- max(str_count(fem_sig$direct, "->")) + 1
male_max <- max(str_count(male_sig$direct, "->")) + 1

# parse trajectories

fem_sig %<>%
  select(direct) %>%
  separate(direct, into = paste0("d", 1:fem_max), sep = "->")

male_sig %<>%
  select(direct) %>%
  separate(direct, into = paste0("d", 1:male_max), sep = "->")

# get morbidity

fem_morb <- apply(fem_traj, 1, function(x) sum(!is.na(x)))
male_morb <- apply(male_traj, 1, function(x) sum(!is.na(x)))

# filter for morbidity > 1

fem_traj <- fem_traj[fem_morb > 1,]
male_traj <- male_traj[male_morb > 1,]

fem_dates <- fem_dates[fem_morb > 1,]
male_dates <- male_dates[male_morb > 1,]

fem_unnorm_dates <- fem_unnorm_dates[fem_morb > 1, ]
male_unnorm_dates <- male_unnorm_dates[male_morb > 1, ]

fem_age <- fem_age[fem_morb > 1]
male_age <- male_age[male_morb > 1]

fem_bl <- as.Date(fem_bl[fem_morb > 1])
male_bl <- as.Date(male_bl[fem_morb > 1])

# convert dates to days from birth (roughly)

fem_unnorm_dates <- as.matrix(fem_unnorm_dates)
male_unnorm_dates <- as.matrix(male_unnorm_dates)

fem_age_dates <- t(sapply(1:nrow(fem_unnorm_dates), function(i) (as.Date(fem_unnorm_dates[i,]) - fem_bl[i]) + (365*(fem_age[i]-1) + 366/2)))
male_age_dates <- t(sapply(1:nrow(male_unnorm_dates), function(i) (as.Date(male_unnorm_dates[i,]) - male_bl[i]) + (365*(male_age[i]-1) + 366/2)))

# convert to matrix

fem_traj <- as.matrix(fem_traj)
male_traj <- as.matrix(male_traj)

fem_dates <- as.matrix(fem_dates)
male_dates <- as.matrix(male_dates)

# get full list of trajectories

fem_traj_all <- rbind(cbind(fem_sig$d1, fem_sig$d2, NA, NA, NA),
                      cbind(fem_sig$d1, fem_sig$d2, fem_sig$d3, NA, NA),
                      cbind(fem_sig$d1, fem_sig$d2, fem_sig$d3, fem_sig$d4, NA))

fem_traj_all <- fem_traj_all[!duplicated.array(fem_traj_all),]

male_traj_all <- rbind(cbind(male_sig$d1, male_sig$d2, NA, NA, NA),
                      cbind(male_sig$d1, male_sig$d2, male_sig$d3, NA, NA),
                      cbind(male_sig$d1, male_sig$d2, male_sig$d3, male_sig$d4, NA),
                      cbind(male_sig$d1, male_sig$d2, male_sig$d3, male_sig$d4, male_sig$d5))

male_traj_all <- male_traj_all[!duplicated.array(male_traj_all),]

# get trajectory counts

fem_counts <- numeric(length = nrow(fem_traj_all))
male_counts <- numeric(length = nrow(male_traj_all))

fem_t <- numeric(length = nrow(fem_traj_all))
male_t <- numeric(length = nrow(male_traj_all))

for(i in 1:nrow(fem_traj_all)) {
  traj_all_i <- unlist(fem_traj_all[i, ])
  traj_all_i <- traj_all_i[!is.na(traj_all_i)]
  ind_i <- apply(fem_traj, 1, function(x) all(traj_all_i %in% x))
  traj_filt <- fem_traj[ind_i,]
  dates_filt <- fem_dates[ind_i,]
  ord <- apply(traj_filt, 1, function(x) all(x[which(x %in% traj_all_i)] == traj_all_i))
  traj_filt <- traj_filt[ord,]
  dates_filt <- dates_filt[ord,]
  dis_ind <- t(apply(traj_filt, 1, function(x) which(x %in% traj_all_i)))
  dis_t <- dis_ind
  for(j in 1:nrow(dates_filt)) {
    dis_t[j,] <- dates_filt[j, dis_ind[j,]]
  }
  same_t <- apply(dis_t, 1, function(x) any(duplicated(x)))
  dis_t <- dis_t[!same_t,]
  dates_filt <- dates_filt[!same_t,]
  dis_ind <- dis_ind[!same_t,]
  # ord_i <- apply(dis_t, 1, function(x) all(x == sort(x)))
  # dis_t <- dis_t[ord_i,]
  # dates_filt <- dates_filt[ord_i,]
  # dis_ind <- dis_ind[ord_i,]
  # morb_i <- numeric(length = nrow(dates_filt))
  # for(j in 1:nrow(dates_filt)) {
  #   t_j <- dates_filt[j,]
  #   ind_j <- dis_ind[j,]
  #   t_j <- t_j[min(which(t_j == t_j[ind_j[1]])):max(which(t_j == max(t_j[ind_j[length(ind_j)]])))]
  #   morb_j <- length(t_j) - length(ind_j)
  #   morb_i[j] <- morb_j
  # }
  t_diff <- apply(dis_t, 1, function(x) x[length(x)] - x[1])
  lt_days <- apply(dis_t, 1, function(x) any(x - lag(x) < n_days, na.rm = T))
  t_diff <- t_diff[!lt_days]
  # morb_i <- morb_i[!lt_days]
  fem_counts[i] <- sum(!lt_days)
  fem_t[i] <- median(t_diff)
  # fem_len[i] <- median(morb_i)
  # fem_mt[i] <- median(morb_i/(t_diff/365))
  print(i)
}

for(i in 1:nrow(male_traj_all)) {
  traj_all_i <- unlist(male_traj_all[i, ])
  traj_all_i <- traj_all_i[!is.na(traj_all_i)]
  ind_i <- apply(male_traj, 1, function(x) all(traj_all_i %in% x))
  traj_filt <- male_traj[ind_i,]
  dates_filt <- male_dates[ind_i,]
  ord <- apply(traj_filt, 1, function(x) all(x[which(x %in% traj_all_i)] == traj_all_i))
  traj_filt <- traj_filt[ord,]
  dates_filt <- dates_filt[ord,]
  dis_ind <- t(apply(traj_filt, 1, function(x) which(x %in% traj_all_i)))
  dis_t <- dis_ind
  for(j in 1:nrow(dates_filt)) {
    dis_t[j,] <- dates_filt[j, dis_ind[j,]]
  }
  same_t <- apply(dis_t, 1, function(x) any(duplicated(x)))
  dis_t <- dis_t[!same_t,]
  dates_filt <- dates_filt[!same_t,]
  dis_ind <- dis_ind[!same_t,]
  # ord_i <- apply(dis_t, 1, function(x) all(x == sort(x)))
  # dis_t <- dis_t[ord_i,]
  # dates_filt <- dates_filt[ord_i,]
  # dis_ind <- dis_ind[ord_i,]
  # morb_i <- numeric(length = nrow(dates_filt))
  # for(j in 1:nrow(dates_filt)) {
  #   t_j <- dates_filt[j,]
  #   ind_j <- dis_ind[j,]
  #   t_j <- t_j[min(which(t_j == t_j[ind_j[1]])):max(which(t_j == max(t_j[ind_j[length(ind_j)]])))]
  #   morb_j <- length(t_j) - length(ind_j)
  #   morb_i[j] <- morb_j
  # }
  t_diff <- apply(dis_t, 1, function(x) x[length(x)] - x[1])
  lt_days <- apply(dis_t, 1, function(x) any(x - lag(x) < n_days, na.rm = T))
  t_diff <- t_diff[!lt_days]
  # morb_i <- morb_i[!lt_days]
  male_counts[i] <- sum(!lt_days)
  male_t[i] <- median(t_diff)
  # male_len[i] <- median(morb_i)
  # male_mt[i] <- median(morb_i/(t_diff/365))
  print(i)
}

# add to df

fem_traj_all %<>%
  as.data.frame() %>%
  mutate(count = fem_counts,
         t = fem_t)

male_traj_all %<>%
  as.data.frame() %>%
  mutate(count = male_counts,
         t = male_t)

# convert to years

fem_traj_all %<>%
  mutate(t = t/365)

male_traj_all %<>%
  mutate(t = t/365)

# write out 

write_tsv(fem_traj_all, "fem_significant_trajectory_info.tsv")
write_tsv(male_traj_all, "male_significant_trajectory_info.tsv")