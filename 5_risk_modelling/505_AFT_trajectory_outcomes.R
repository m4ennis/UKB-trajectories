# calculate rates of death + hospitalisation of significant built trajectories using Accelerated Failure Time (AFT) models in parallel

library(tidyverse)
library(data.table)
library(magrittr)
library(parallel)
library(foreach)
library(doParallel)
library(survival)
library(broom)

# variables

# female diagnosis trajectory matrix
fem_traj_file = "fem_disease_block_1pct_traj.tsv"
# male diagnosis trajectory matrix
male_traj_file = "male_disease_block_1pct_traj.tsv"

# female time to diagnosis matrix
fem_dates_file = "fem_diag_block_1pct_norm_dates.tsv"
# male time to diagnosis matrix
male_dates_file = "male_diag_block_1pct_norm_dates.tsv"

# female date of diagnosis matrix
fem_unnorm_dates_file = "fem_diag_block_1pct_dates.tsv"
# male date of diagnosis matrix
male_unnorm_dates_file = "male_diag_block_1pct_dates.tsv"

# ukb file containing age at attending UK Biobank assessment centre and participant sex
ukb_file = "ukb_df_general_col_filtered_all.tsv"
# ukb file containing ICD10 diagnoses (f.41270) and dates (41280)
ukb_file2 = "ukb_df_t2d_col_filtered_all.tsv"
# ukb file containing date of attending assessment centre
ukb_file3 = "ukb_df_t2d_date_attendence_filtered.tsv"
# ukb file containing dates of death of participants
death_file = "ukb_df_t2d_diag_date_death_filtered.tsv"

# output files

# female significant built trajectories with added info
fem_sig_file = "female_significant_trajectory_info.tsv"
# male significant built trajectories with added info
male_sig_file = "male_significant_trajectory_info.tsv"

# ICD10 block dictionary
icd10_blocks_file = "icd10_blocks.tsv"

# functions

traj_risk_hosp <- function(traj, traj_mat, dates_mat, dates_all_mat, age_mat, death, th, fu, n_hosp) {
  inds <- sapply(traj, function(x) apply(traj_mat, 1, function(y) x %in% y)) # get indices of diseases
  ind <- inds[, ncol(inds)] # get those with presenting disease
  
  # filter for those with presenting disease
  d <- traj_mat[ind,]
  t <- dates_mat[ind,]
  t_all <- dates_all_mat[ind,]
  age <- age_mat[ind,]
  dead <- death[ind]
  traj_fu <- fu[ind]
  
  dis_t <- sapply(1:nrow(d), function(x) t[x, which(d[x,] == traj[length(traj)])])
  
  # get time of trajectory
  
  traj_t <- sapply(traj[-length(traj)], function(x) apply(d, 1, function(y) ifelse(x %in% y, which(y == x), NA)))
  if(length(traj) == 2) {
    traj_t <- sapply(1:nrow(d), function(i) t[i, traj_t[i,]])
    traj_t <- as.matrix(traj_t)
    
    # get time difference between presenting disease
    
    t_diff <- sapply(1:nrow(d), function(i) traj_t[i,] - dis_t[i])
    t_diff <- as.matrix(t_diff)
    
  } else {
    traj_t <- t(sapply(1:nrow(d), function(i) t[i, traj_t[i,]]))
    t_diff <- t(sapply(1:nrow(d), function(i) traj_t[i,] - dis_t[i]))
  }
 
  # get age at last in comparison
  
  age <- age[cbind(1:nrow(age), apply(d, 1, function(x) which(x == traj[length(traj)])))]
  
  t_start <- t[cbind(1:nrow(t), apply(d, 1, function(x) which(x == traj[length(traj)])))]
  event_dt <- sapply(1:length(t_start), function(i) min(t_all[i,][which(t_all[i,] == unique(t_all[i,][t_all[i,] > t_start[i]], na.rm = T)[n_hosp])]))
  event_dt <- event_dt - t_start
  dead <- dead - t_start
  dead <- ifelse(dead <= 0, 1, dead)
  traj_fu <- traj_fu - t_start
  # traj_fu <- ifelse(traj_fu <= 0, 1, traj_fu)
  status <- ifelse(is.infinite(event_dt) | event_dt > (th*365), 0, 1)
  dead_ind <- which(dead < event_dt)
  event_dt[dead_ind] <- dead[dead_ind]
  event_dt[is.infinite(event_dt)] <- th*365
  event_dt[event_dt > (th*365)] <- th*365
  event_dt[event_dt > traj_fu] <- traj_fu[event_dt > traj_fu]
  
  surv_df <- data.frame(time = c(event_dt), status = c(status))
  surv_df$time <- ifelse(surv_df$time == 0, surv_df$time + 1, surv_df$time)
  surv_df$time <- ifelse(surv_df$time < 0, 1, surv_df$time)
  
  prior_dis <- apply(t_diff, 2, function(x) ifelse(x > 0 | is.na(x), 0, 1))
  prior_dis <- as.matrix(prior_dis)
  
  prior_traj <- d
  for(i in 1:nrow(d)) {
    prior_traj[i, which(t[i,] > dis_t[i])] <- NA
  }
  
  prior_morb <- apply(prior_traj, 1, function(x) sum(!is.na(x)))
  age <- age/365.25
  
  surv_df %<>% mutate(morb = prior_morb - 1, # correct for presenting disease 
                      age = age)
  
  colnames(prior_dis) <- c(paste0("hist", 1:ncol(prior_dis)))
  morb_adj <- apply(prior_dis, 1, sum)
  prior_dis <- apply(prior_dis, 2, as.factor)
  surv_df <- cbind(surv_df, prior_dis) %>%
    mutate(morb = morb - morb_adj) # correct for prior histories
  model_form <- as.formula(paste0("Surv(time = time, event = status) ~ age + morb + ", paste0("(", paste0("hist", 1:ncol(prior_dis), collapse = " + "), ")^3")))
  
  survreg_model <- survreg(model_form, 
                           data = surv_df, 
                           dist = "lognormal")
  survreg_out <- tidy(survreg_model)
  survreg_out %<>%
    mutate(traj = paste0(traj, collapse = " -> "),
           th = th) %>%
    select(traj, everything())
  
  # add in population sizes
  
  survreg_out %<>%
    mutate(pres_size = nrow(surv_df),
           hist_size = sum(apply(prior_dis, 1, function(x) all(x == 1))),
           traj_size = sum(apply(t_diff, 1, function(x) all(x < 0 & order(x) == 1:length(x))), na.rm = T))
  
  return(survreg_out)
}

morb_risk_hosp <- function(disease, traj_mat, dates_mat, dates_all_mat, age_mat, death, th, fu, n_hosp) {
  ind <- apply(traj_mat, 1, function(y) disease %in% y) # get indices of disease

  # filter for those with presenting disease
  d <- traj_mat[ind,]
  t <- dates_mat[ind,]
  t_all <- dates_all_mat[ind,]
  age <- age_mat[ind,]
  dead <- death[ind]
  traj_fu <- fu[ind]
  
  dis_t <- sapply(1:nrow(d), function(x) t[x, which(d[x,] == disease)])
  
  # get age at last in comparison
  
  age <- age[cbind(1:nrow(age), apply(d, 1, function(x) which(x == disease)))]

  t_start <- t[cbind(1:nrow(t), apply(d, 1, function(x) which(x == disease)))]
  event_dt <- sapply(1:length(t_start), function(i) min(t_all[i,][which(t_all[i,] == unique(t_all[i,][t_all[i,] > t_start[i]], na.rm = T)[n_hosp])]))
  event_dt <- event_dt - t_start
  dead <- dead - t_start
  dead <- ifelse(dead <= 0, 1, dead)
  traj_fu <- traj_fu - t_start
  # traj_fu <- ifelse(traj_fu <= 0, 1, traj_fu)
  status <- ifelse(is.infinite(event_dt) | event_dt > (th*365), 0, 1)
  dead_ind <- which(dead < event_dt)
  event_dt[dead_ind] <- dead[dead_ind]
  event_dt[is.infinite(event_dt)] <- th*365
  event_dt[event_dt > (th*365)] <- th*365
  event_dt[event_dt > traj_fu] <- traj_fu[event_dt > traj_fu]
  
  surv_df <- data.frame(time = c(event_dt), status = c(status))
  surv_df$time <- ifelse(surv_df$time == 0, surv_df$time + 1, surv_df$time)
  surv_df$time <- ifelse(surv_df$time < 0, 1, surv_df$time)

  prior_traj <- d
  for(i in 1:nrow(d)) {
    prior_traj[i, which(t[i,] > dis_t[i])] <- NA
  }
  
  prior_morb <- apply(prior_traj, 1, function(x) sum(!is.na(x)))
  age <- age/365.25
  
  surv_df %<>% mutate(morb = prior_morb - 1, # correct for presenting disease
                      age = age)

  model_form <- as.formula(paste0("Surv(time = time, event = status) ~ age + morb"))
  
  survreg_model <- survreg(model_form, 
                           data = surv_df, 
                           dist = "lognormal")
  survreg_out <- tidy(survreg_model)
  survreg_out %<>%
    mutate(disease = disease,
           th = th) %>%
    select(disease, everything())
  
  # add in population size
  
  survreg_out %<>%
    mutate(pres_size = nrow(surv_df))
  
  return(survreg_out)
}


traj_risk_death <- function(traj, traj_mat, dates_mat, dates_all_mat, age_mat, death, th, fu) {
 
  inds <- sapply(traj, function(x) apply(traj_mat, 1, function(y) x %in% y)) # get indices of diseases
  ind <- inds[, ncol(inds)] # get those with presenting disease
  
  # filter for those with presenting disease
  d <- traj_mat[ind,]
  t <- dates_mat[ind,]
  t_all <- dates_all_mat[ind,]
  age <- age_mat[ind,]
  dead <- death[ind]
  traj_fu <- fu[ind]
  
  dis_t <- sapply(1:nrow(d), function(x) t[x, which(d[x,] == traj[length(traj)])])
  
  # get time of trajectory
  
  traj_t <- sapply(traj[-length(traj)], function(x) apply(d, 1, function(y) ifelse(x %in% y, which(y == x), NA)))
  if(length(traj) == 2) {
    traj_t <- sapply(1:nrow(d), function(i) t[i, traj_t[i,]])
    traj_t <- as.matrix(traj_t)
    
    # get time difference between presenting disease
    
    t_diff <- sapply(1:nrow(d), function(i) traj_t[i,] - dis_t[i])
    t_diff <- as.matrix(t_diff)
    
  } else {
    traj_t <- t(sapply(1:nrow(d), function(i) t[i, traj_t[i,]]))
    t_diff <- t(sapply(1:nrow(d), function(i) traj_t[i,] - dis_t[i]))
  }
  
  # get age at last in comparison
  
  age <- age[cbind(1:nrow(age), apply(d, 1, function(x) which(x == traj[length(traj)])))]
  
  t_start <- t[cbind(1:nrow(t), apply(d, 1, function(x) which(x == traj[length(traj)])))]
  event_dt <- dead - t_start
  status <- ifelse(event_dt > (th*365) | is.na(event_dt), 0, 1)
  event_dt[event_dt > (th*365) & !is.na(event_dt)] <- (th*365)
  traj_fu <- traj_fu - t_start
  event_dt[traj_fu < (th*365) & is.na(event_dt)] <- traj_fu[traj_fu < (th*365) & is.na(event_dt)]
  event_dt[traj_fu >= (th*365) & is.na(event_dt)] <- (th*365)
  
  surv_df <- data.frame(time = c(event_dt), status = c(status))
  surv_df$time <- ifelse(surv_df$time == 0, surv_df$time + 1, surv_df$time)
  surv_df$time <- ifelse(surv_df$time < 0, 1, surv_df$time)
  
  prior_dis <- apply(t_diff, 2, function(x) ifelse(x > 0 | is.na(x), 0, 1))
  prior_dis <- as.matrix(prior_dis)
  
  prior_traj <- d
  for(i in 1:nrow(d)) {
    prior_traj[i, which(t[i,] > dis_t[i])] <- NA
  }
  
  prior_morb <- apply(prior_traj, 1, function(x) sum(!is.na(x)))
  age <- age/365.25
  
  surv_df %<>% mutate(morb = prior_morb - 1, # correct for presenting disease
                      age = age)
  
  colnames(prior_dis) <- c(paste0("hist", 1:ncol(prior_dis)))
  morb_adj <- apply(prior_dis, 1, sum)
  prior_dis <- apply(prior_dis, 2, as.factor)
  surv_df <- cbind(surv_df, prior_dis) %>%
    mutate(morb = morb - morb_adj) # correct for prior histories
  model_form <- as.formula(paste0("Surv(time = time, event = status) ~ age + morb + ", paste0("(", paste0("hist", 1:ncol(prior_dis), collapse = " + "), ")^3")))
  
  survreg_model <- survreg(model_form, 
                           data = surv_df, 
                           dist = "lognormal")
  survreg_out <- tidy(survreg_model)
  survreg_out %<>%
    mutate(traj = paste0(traj, collapse = " -> "),
           th = th) %>%
    select(traj, everything())
  
  # add in population sizes
  
  survreg_out %<>%
    mutate(pres_size = nrow(surv_df),
           hist_size = sum(apply(prior_dis, 1, function(x) all(x == 1))),
           traj_size = sum(apply(t_diff, 1, function(x) all(x < 0 & order(x) == 1:length(x))), na.rm = T))
  
  return(survreg_out)
}

morb_risk_death <- function(disease, traj_mat, dates_mat, dates_all_mat, age_mat, death, th, fu) {
  
  ind <- apply(traj_mat, 1, function(y) disease %in% y) # get indices of disease

  # filter for those with presenting disease
  d <- traj_mat[ind,]
  t <- dates_mat[ind,]
  t_all <- dates_all_mat[ind,]
  age <- age_mat[ind,]
  dead <- death[ind]
  traj_fu <- fu[ind]
  
  dis_t <- sapply(1:nrow(d), function(x) t[x, which(d[x,] == disease)])
  
  # get age at last in comparison
  
  age <- age[cbind(1:nrow(age), apply(d, 1, function(x) which(x == disease)))]
  
  t_start <- t[cbind(1:nrow(t), apply(d, 1, function(x) which(x == disease)))]
  event_dt <- dead - t_start
  status <- ifelse(event_dt > (th*365) | is.na(event_dt), 0, 1)
  event_dt[event_dt > (th*365) & !is.na(event_dt)] <- (th*365)
  traj_fu <- traj_fu - t_start
  event_dt[traj_fu < (th*365) & is.na(event_dt)] <- traj_fu[traj_fu < (th*365) & is.na(event_dt)]
  event_dt[traj_fu >= (th*365) & is.na(event_dt)] <- (th*365)
  
  surv_df <- data.frame(time = c(event_dt), status = c(status))
  surv_df$time <- ifelse(surv_df$time == 0, surv_df$time + 1, surv_df$time)
  surv_df$time <- ifelse(surv_df$time < 0, 1, surv_df$time)
  
  prior_traj <- d
  for(i in 1:nrow(d)) {
    prior_traj[i, which(t[i,] > dis_t[i])] <- NA
  }
  
  prior_morb <- apply(prior_traj, 1, function(x) sum(!is.na(x)))
  age <- age/365.25
  
  surv_df %<>% mutate(morb = prior_morb - 1, # correct for presenting disease
                      age = age)

  model_form <- as.formula(paste0("Surv(time = time, event = status) ~ age + morb"))
  
  survreg_model <- survreg(model_form, 
                           data = surv_df, 
                           dist = "lognormal")
  survreg_out <- tidy(survreg_model)
  survreg_out %<>%
    mutate(disease = disease,
           th = th) %>%
    select(disease, everything())
  
  # add in population size
  
  survreg_out %<>%
    mutate(pres_size = nrow(surv_df))
  
  return(survreg_out)
}

# read in trajs

fem_traj <- fread(fem_traj_file)
male_traj <- fread(male_traj_file)

# read in traj times

fem_dates <- fread(fem_dates_file)
male_dates <- fread(male_dates_file)

fem_traj <- as.matrix(fem_traj)
male_traj <- as.matrix(male_traj)

fem_dates <- as.matrix(fem_dates)
male_dates <- as.matrix(male_dates)

# read in trajectory dates

fem_unnorm_dates <- fread(fem_unnorm_dates_file)
male_unnorm_dates <- fread(male_unnorm_dates_file)

fem_fu_date <- max(as.Date(as.vector(as.matrix(fem_unnorm_dates))), na.rm = T)
male_fu_date <- max(as.Date(as.vector(as.matrix(male_unnorm_dates))), na.rm = T)

# get maximum date of follow-up

# read in UKB file

ukb <- fread(ukb_file)

age <- ukb$age_when_attended_assessment_centre_f21003_0_0

fem_age <- age[ukb$sex_f31_0_0 == 0]
male_age <- age[ukb$sex_f31_0_0 == 1]

# read in baseline dates

ukb3 <- fread(ukb_file3)

bl_date <- ukb3$date_of_attending_assessment_centre_f53_0_0
fem_bl <- bl_date[ukb$sex_f31_0_0 == 0]
male_bl <- bl_date[ukb$sex_f31_0_0 == 1]

# read in death

death_all <- fread(death_file)
death_all <- death_all$date_of_death_f40000_0_0

# filter death

fem_death <- death_all[ukb$sex_f31_0_0 == 0]
male_death <- death_all[ukb$sex_f31_0_0 == 1]

# read in icd10 blocks

icd10_blocks <- fread(icd10_blocks_file)

# get morbidity

fem_morb <- apply(fem_traj, 1, function(x) sum(!is.na(x)))
male_morb <- apply(male_traj, 1, function(x) sum(!is.na(x)))

# filter for morbidity > 1

fem_traj <- fem_traj[fem_morb > 0,]
male_traj <- male_traj[male_morb > 0,]

fem_dates <- fem_dates[fem_morb > 0,]
male_dates <- male_dates[male_morb > 0,]

fem_death <- fem_death[fem_morb > 0]
male_death <- male_death[male_morb > 0]

fem_unnorm_dates <- fem_unnorm_dates[fem_morb > 0, ]
male_unnorm_dates <- male_unnorm_dates[male_morb > 0, ]

fem_age <- fem_age[fem_morb > 0]
male_age <- male_age[male_morb > 0]

fem_bl <- as.Date(fem_bl[fem_morb > 0])
male_bl <- as.Date(male_bl[fem_morb > 0])

# convert dates to days from birth (roughly)

fem_unnorm_dates <- as.matrix(fem_unnorm_dates)
male_unnorm_dates <- as.matrix(male_unnorm_dates)

fem_age_dates <- t(sapply(1:nrow(fem_unnorm_dates), function(i) (as.Date(fem_unnorm_dates[i,]) - fem_bl[i]) + (365*(fem_age[i]-1) + 366/2)))
male_age_dates <- t(sapply(1:nrow(male_unnorm_dates), function(i) (as.Date(male_unnorm_dates[i,]) - male_bl[i]) + (365*(male_age[i]-1) + 366/2)))

# get follow-up time

fem_fu <- sapply(1:nrow(fem_unnorm_dates), function(i) fem_fu_date - min(as.Date(fem_unnorm_dates[i, !is.na(fem_unnorm_dates[i,])])))
male_fu <- sapply(1:nrow(male_unnorm_dates), function(i) male_fu_date - min(as.Date(male_unnorm_dates[i, !is.na(male_unnorm_dates[i,])])))

# convert death dates to relative time

fem_death <- sapply(1:nrow(fem_unnorm_dates), function(i) as.Date(fem_death[i]) - as.Date(fem_unnorm_dates[i, 1]))
male_death <- sapply(1:nrow(male_unnorm_dates), function(i) as.Date(male_death[i]) - as.Date(male_unnorm_dates[i, 1]))

# convert to matrix

fem_traj <- as.matrix(fem_traj)
male_traj <- as.matrix(male_traj)

fem_dates <- as.matrix(fem_dates)
male_dates <- as.matrix(male_dates)

# read in significant trajs

fem_sig <- fread(fem_sig_file)
male_sig <- fread(male_sig_file)

# get diagnoses

# read in UKB file 2

ukb2 <- fread(ukb_file2)

diags <- ukb2 %>%
  select(contains("41270"))

dates <- ukb2 %>%
  select(contains("41280"))

# convert dates to numeric

dates %<>%
  mutate(across(.fns = as.numeric))

# convert to matrix

diags <- as.matrix(diags)

dates <- as.matrix(dates)

# convert to 3 letter code

diags <- str_trunc(diags, 3, ellipsis = "")

# convert to relative time since first diagnosis

dates <- t(apply(dates, 1, function(x) x - min(x, na.rm = T)))

# filter dates for sex + morbidity

fem_dates_all <- as.matrix(dates[ukb$sex_f31_0_0 == 0,][fem_morb > 0,])
male_dates_all <- as.matrix(dates[ukb$sex_f31_0_0 == 1,][male_morb > 0,])

fem_diags_all <- as.matrix(diags[ukb$sex_f31_0_0 == 0,][fem_morb > 0,])
male_diags_all <- as.matrix(diags[ukb$sex_f31_0_0 == 1,][male_morb > 0,])

# get diseases which we are looking at

dis_use <- sort(unique(c(unlist(fem_traj), unlist(male_traj))))

# convert diags to blocks

fem_diags_all <- t(apply(fem_diags_all, 1, function(x) icd10_blocks$block[match(x, icd10_blocks$code)]))
male_diags_all <- t(apply(male_diags_all, 1, function(x) icd10_blocks$block[match(x, icd10_blocks$code)]))

# get time of earliest diagnostic code

fem_min_t <- sapply(1:nrow(fem_dates_all), function(x) min(fem_dates_all[x,][which(fem_diags_all[x,] %in% dis_use)], na.rm = T))
male_min_t <- sapply(1:nrow(male_dates_all), function(x) min(male_dates_all[x,][which(male_diags_all[x,] %in% dis_use)], na.rm = T))

# order dates all

fem_dates_all <- t(apply(fem_dates_all, 1, function(x) x[order(x)]))
male_dates_all <- t(apply(male_dates_all, 1, function(x) x[order(x)]))

fem_dates_all <- t(sapply(1:nrow(fem_dates_all), function(x) fem_dates_all[x,] - fem_min_t[x]))
male_dates_all <- t(sapply(1:nrow(male_dates_all), function(x) male_dates_all[x,] - male_min_t[x]))

# 1 year relative risk of hospitalisation for each trajectory

fem_nodes <- rbind(cbind(fem_sig$`Disease 1`, fem_sig$`Disease 2`, NA, NA), cbind(fem_sig$`Disease 1`, fem_sig$`Disease 2`, fem_sig$`Disease 3`, NA),
                   cbind(fem_sig$`Disease 1`, fem_sig$`Disease 2`, fem_sig$`Disease 3`, fem_sig$`Disease 4`))

fem_nodes <- unique.array(fem_nodes)

# female prior morbidity risk hospitalisation

fem_pres <- sort(unique(as.vector(fem_traj)))

fem_morb_1_hosp_coefs <- lapply(fem_pres, function(disease) morb_risk_hosp(disease, fem_traj, fem_dates, fem_dates_all, fem_age_dates, fem_death, 1, fem_fu, 1))

fem_morb_1_hosp_coefs <- rbindlist(fem_morb_1_hosp_coefs)

fem_morb_1_hosp_coefs %<>%
  mutate(adj_p_val = p.adjust(p.value, method = "holm"))

# trajectory risk

fem_risk_1_hosp_coefs <- lapply(1:nrow(fem_nodes), function(i) {
  traj_i <- fem_nodes[i,] %>% unlist()
  traj_i <- traj_i[!is.na(traj_i)]
  risk_i <- traj_risk_hosp(traj_i, fem_traj, fem_dates, fem_dates_all, fem_age_dates, fem_death, 1, fem_fu, 1)
  print(i/nrow(fem_nodes))
  return(risk_i)
})

fem_risk_1_hosp_coefs <- rbindlist(fem_risk_1_hosp_coefs)

fem_risk_1_hosp_coefs %<>%
  mutate(adj_p_val = p.adjust(p.value, method = "holm"))

# 1 year risk of mortality of trajectories

# female prior morbidity risk death

fem_morb_death_coefs <- lapply(fem_pres, function(disease) morb_risk_death(disease, fem_traj, fem_dates, fem_dates_all, fem_age_dates, fem_death, 1, fem_fu))

fem_morb_death_coefs <- rbindlist(fem_morb_death_coefs)

fem_morb_death_coefs %<>%
  mutate(adj_p_val = p.adjust(p.value, method = "holm"))

# trajectory risk

fem_risk_death_coefs <- lapply(1:nrow(fem_nodes), function(i) {
  traj_i <- fem_nodes[i,] %>% unlist()
  traj_i <- traj_i[!is.na(traj_i)]
  risk_i <- traj_risk_death(traj_i, fem_traj, fem_dates, fem_dates_all, fem_age_dates, fem_death, 1, fem_fu)
  print(i/nrow(fem_nodes))
  return(risk_i)
})

fem_risk_death_coefs <- rbindlist(fem_risk_death_coefs)

fem_risk_death_coefs %<>%
  mutate(adj_p_val = p.adjust(p.value, method = "holm"))

# males

# 1 year relative risk of hospitalisation for each trajectory

male_nodes <- rbind(cbind(male_sig$`Disease 1`, male_sig$`Disease 2`, NA, NA), cbind(male_sig$`Disease 1`, male_sig$`Disease 2`, male_sig$`Disease 3`, NA),
                   cbind(male_sig$`Disease 1`, male_sig$`Disease 2`, male_sig$`Disease 3`, male_sig$`Disease 4`))

male_nodes <- unique.array(male_nodes)

# male prior morbidity risk hospitalisation

male_pres <- sort(unique(as.vector(male_traj)))

male_morb_1_hosp_coefs <- lapply(male_pres, function(disease) morb_risk_hosp(disease, male_traj, male_dates, male_dates_all, male_age_dates, male_death, 1, male_fu, 1))

male_morb_1_hosp_coefs <- rbindlist(male_morb_1_hosp_coefs)

male_morb_1_hosp_coefs %<>%
  mutate(adj_p_val = p.adjust(p.value, method = "holm"))

# trajectory risk

male_risk_1_hosp_coefs <- lapply(1:nrow(male_nodes), function(i) {
  traj_i <- male_nodes[i,] %>% unlist()
  traj_i <- traj_i[!is.na(traj_i)]
  risk_i <- traj_risk_hosp(traj_i, male_traj, male_dates, male_dates_all, male_age_dates, male_death, 1, male_fu, 1)
  print(i/nrow(male_nodes))
  return(risk_i)
})

male_risk_1_hosp_coefs <- rbindlist(male_risk_1_hosp_coefs)

male_risk_1_hosp_coefs %<>%
  mutate(adj_p_val = p.adjust(p.value, method = "holm"))

# 1 year risk of mortality of trajectories

# male prior morbidity risk death

male_morb_death_coefs <- lapply(male_pres, function(disease) morb_risk_death(disease, male_traj, male_dates, male_dates_all, male_age_dates, male_death, 1, male_fu))

male_morb_death_coefs <- rbindlist(male_morb_death_coefs)

male_morb_death_coefs %<>%
  mutate(adj_p_val = p.adjust(p.value, method = "holm"))

# trajectory risk

male_risk_death_coefs <- lapply(1:nrow(male_nodes), function(i) {
  traj_i <- male_nodes[i,] %>% unlist()
  traj_i <- traj_i[!is.na(traj_i)]
  risk_i <- traj_risk_death(traj_i, male_traj, male_dates, male_dates_all, male_age_dates, male_death, 1, male_fu)
  print(i/nrow(male_nodes))
  return(risk_i)
})

male_risk_death_coefs <- rbindlist(male_risk_death_coefs)

male_risk_death_coefs %<>%
  mutate(adj_p_val = p.adjust(p.value, method = "holm"))

fem_risk_1_df <- fem_risk_1_hosp_coefs
fem_risk_1_df %<>%
  mutate(n_hosp = 1)

fem_risk_death_df <- fem_risk_death_coefs

write_tsv(fem_risk_1_df, "fem_traj_risk_hosp_coefs_new.tsv")

write_tsv(fem_risk_death_df, "fem_traj_risk_death_coefs_new.tsv")

write_tsv(fem_morb_1_hosp_coefs, "fem_morb_risk_hosp_coefs_new.tsv")

write_tsv(fem_morb_death_coefs, "fem_morb_risk_death_coefs_new.tsv")

# male

male_risk_1_df <- male_risk_1_hosp_coefs
male_risk_1_df %<>%
  mutate(n_hosp = 1)

male_risk_death_df <- male_risk_death_coefs

write_tsv(male_risk_1_df, "male_traj_risk_hosp_coefs_new.tsv")

write_tsv(male_risk_death_df, "male_traj_risk_death_coefs_new.tsv")

write_tsv(male_morb_1_hosp_coefs, "male_morb_risk_hosp_coefs_new.tsv")

write_tsv(male_morb_death_coefs, "male_morb_risk_death_coefs_new.tsv")
