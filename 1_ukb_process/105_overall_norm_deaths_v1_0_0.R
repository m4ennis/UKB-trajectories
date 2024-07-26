# normalise death dates

print("Loading libraries...")

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(data.table)

# Variables

# file containing date of death for individuals who have died
ukb_icd10_diag_date_file = ""

# column name of column containing date of death
death_col = ""

# female diagnosis dates matrix
fem_traj_dates_file = "fem_diag_block_1pct_dates.tsv"

# male diagnosis dates matrix
male_traj_dates_file = "male_diag_block_1pct_dates.tsv"

# file containing sex of participants
ukb_sex_file = ""

# column name of column containing sex of participants
sex_col <- ""

# output files

# female time of death relative to first diagnosis matrix
fem_norm_death_file = "fem_norm_death.tsv"
# male time of death relative to first diagnosis matrix
male_norm_death_file = "male_norm_death.tsv"

# female time of follow-up end relative to first diagnosis matrix
fem_norm_fu_file = "fem_norm_follow_up.tsv"
# male time of follow-up end relative to first diagnosis matrix
male_norm_fu_file = "male_norm_follow_up.tsv"

# read in temporal data

temporal <- read_tsv(ukb_icd10_diag_date_file)

ukb_sex <- read_tsv(ukb_sex_file) %>% pull(sex_col)

# get sex indices

fem_ind <- which(ukb_sex == 0)
male_ind <- which(ukb_sex == 1)

# get deaths

deaths <- temporal %>% pull(death_col)

# sex deaths

fem_deaths <- deaths[fem_ind]
male_deaths <- deaths[male_ind]

# read in dates

fem_dates <- fread(fem_traj_dates_file)
male_dates <- fread(male_traj_dates_file)

# get first diagnosis

fem_dates_1st <- pull(fem_dates, 1)
male_dates_1st <- pull(male_dates, 1)

# normalise deaths to first diagnosis

fem_norm_deaths <- fem_deaths - as.Date(fem_dates_1st)
male_norm_deaths <- male_deaths - as.Date(male_dates_1st)

# get max follow-up date

max_fu <- as.Date(max(c(unlist(fem_dates), unlist(male_dates)), na.rm = T))
  
# get individual follow-up

fem_fu <- rep(max_fu, nrow(fem_dates))
fem_fu[!is.na(fem_deaths)] <- fem_deaths[!is.na(fem_deaths)]

male_fu <- rep(max_fu, nrow(male_dates))
male_fu[!is.na(male_deaths)] <- male_deaths[!is.na(male_deaths)]

fem_norm_fu <- fem_fu - as.Date(fem_dates_1st)
male_norm_fu <- male_fu - as.Date(male_dates_1st)

fem_norm_fu[is.na(fem_norm_fu)] <- 0
male_norm_fu[is.na(male_norm_fu)] <- 0

# female median follow-up
median(as.numeric(fem_norm_fu))/365.5

# male median follow-up
median(as.numeric(male_norm_fu))/365.5

# treat those with no diagnoses and die as dead at 0

fem_dead_0_ind <- which(is.na(fem_dates_1st) & !is.na(fem_deaths))
male_dead_0_ind <- which(is.na(male_dates_1st) & !is.na(male_deaths))

fem_norm_deaths[fem_dead_0_ind] <- min(fem_norm_deaths[!is.na(fem_norm_deaths)])
male_norm_deaths[male_dead_0_ind] <- min(male_norm_deaths[!is.na(male_norm_deaths)])

# write out time normalised deaths

write_tsv(as.data.frame(fem_norm_deaths), fem_norm_death_file)
write_tsv(as.data.frame(male_norm_deaths), male_norm_death_file)

# write out normalised follow-up

write_tsv(as.data.frame(fem_norm_fu), fem_norm_fu_file)
write_tsv(as.data.frame(male_norm_fu), male_norm_fu_file)
