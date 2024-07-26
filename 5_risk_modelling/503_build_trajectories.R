# build significant, longer trajectories using significant pairwise directional trajectories

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

# female pairwise directionality test statistics
fem_sig_dir_file = "fem_sig_pairwise_directional_diseases.tsv"
# male pairwise directionality test statistics
male_sig_dir_file = "male_sig_pairwise_directional_diseases.tsv"

# female built trajectories
fem_out_file = "fem_sig_built_trajectories.tsv"
# male built trajectories
male_out_file = "male_sig_built_trajectories.tsv"

# functions

# trajectory functions

get_traj_counts <- function(x_sig, x_traj, x_dates, x_ind_list, test_len, min_size, n_cpu = NULL) {
  if(is.null(n_cpu)) {
    x_direct_df <- data.frame(matrix(NA, nrow = 0, ncol = 6))
    colnames(x_direct_df) <- c("traj_a", "traj_b", "no_overlap", "a_to_b", "b_to_a", "eq_ab")
    # get overlapping trajs
    tested <- vector(mode = "list")
    finish <- 0
    while(finish != 1) {
      for(i in 1:nrow(x_sig)) {
        cnt <- 1
        x_direct <- data.frame(matrix(NA, nrow = 0, ncol = 6))
        colnames(x_direct) <- c("traj_a", "traj_b", "no_overlap", "a_to_b", "b_to_a", "eq_ab")
        traj_i <- unlist(x_sig[i,])
        if(test_len == 1) {
          for(j in 1:nrow(x_sig)) {
            if(i == j) {
              next
            }
            traj_j <- unlist(x_sig[j,])
            is_in <- any(sapply(tested, function(x) all(x %in% c(traj_i, traj_j))))
            if(is_in) {
              next
            } else {
              tested <- append(tested, list(c(traj_i, traj_j)))
              # filter for traj 
              test_ind <- which(x_ind_list[[which(names(x_ind_list) == traj_i)]] & x_ind_list[[which(names(x_ind_list) == traj_j)]])
              traj_filt <- x_traj[test_ind,]
              dates_filt <- x_dates[test_ind,]
              if(length(test_ind) < min_size) {
                x_direct[cnt,] <- data.frame(traj_a = paste0(traj_i, collapse = "->"), traj_b = paste0(traj_j, collapse = "->"), no_overlap = NA, NA, NA, NA)
                next
              }
              ind_a <- apply(traj_filt, 1, function(x) which(x == traj_i))
              ind_b <- apply(traj_filt, 1, function(x) which(x == traj_j))
              t_a <- dates_filt[cbind(1:nrow(dates_filt), ind_a)]
              t_b <- dates_filt[cbind(1:nrow(dates_filt), ind_b)]
              x_direct[cnt,] <- data.frame(traj_a = paste0(traj_i, collapse = "->"), traj_b = paste0(traj_j, collapse = "->"), no_overlap = paste0("a_", 1, "_b_", 1), a_to_b = sum(t_a < t_b), b_to_a = sum(t_a > t_b), eq_ab = sum(t_a == t_b))
              cnt <- cnt + 1
            }
          }
        } else {
          ol_i <- apply(x_sig, 1, function(x) any(x %in% traj_i))
          ol_i <- which(ol_i)
          ol_i <- ol_i[ol_i != i]
          ol_use_i <- x_sig[ol_i, ]
          for(j in 1:nrow(ol_use_i)) {
            traj_j <- unlist(ol_use_i[j,])
            is_in <- any(sapply(tested, function(x) all(x %in% c(traj_i, traj_j))))
            if(is_in) {
              next
            } else {
              tested <- append(tested, list(c(traj_i, traj_j)))
              # where is overlap?
              where_ol_i <- which(traj_i %in% traj_j)
              where_ol_j <- which(traj_j %in% traj_i)
              if(length(where_ol_j) < test_len - 1 | all(where_ol_i != sort(where_ol_i))) {
                next
              }
              # where is non-overlapping?
              where_not_i <- which(!traj_i %in% traj_j)
              where_not_j <- which(!traj_j %in% traj_i)
              # filter for traj 
              test_ind <- which(apply(sapply(unique(c(traj_i, traj_j)), function(x) x_ind_list[[which(names(x_ind_list) == x)]]), 1, all))
              traj_filt <- x_traj[test_ind,]
              dates_filt <- x_dates[test_ind,]
              if(length(test_ind) < min_size) {
                x_direct[cnt,] <- data.frame(traj_a = paste0(traj_i, collapse = "->"), traj_b = paste0(traj_j, collapse = "->"), no_overlap = NA, NA, NA, NA)
                next
              }
              # order filter
              traj_ord <- traj_i[where_ol_i]
              traj_ord <- sapply(traj_ord, function(x) dates_filt[cbind(1:nrow(dates_filt), apply(traj_filt, 1, function(y) which(y == x)))])
              t_filt <- apply(traj_ord, 1, function(x) all(x == sort(x) & !duplicated(x)))
              if(sum(t_filt) < min_size) {
                x_direct[cnt,] <- data.frame(traj_a = paste0(traj_i, collapse = "->"), traj_b = paste0(traj_j, collapse = "->"), no_overlap = NA, NA, NA, NA)
                next
              }
              traj_filt <- traj_filt[t_filt,]
              dates_filt <- dates_filt[t_filt,]
              dis_a <- traj_i[where_not_i]
              dis_b <- traj_j[where_not_j]
              ind_a <- apply(traj_filt, 1, function(x) which(x == dis_a))
              ind_b <- apply(traj_filt, 1, function(x) which(x == dis_b))
              t_a <- dates_filt[cbind(1:nrow(dates_filt), ind_a)]
              t_b <- dates_filt[cbind(1:nrow(dates_filt), ind_b)]
              x_direct[cnt,] <- data.frame(traj_a = paste0(traj_i, collapse = "->"), traj_b = paste0(traj_j, collapse = "->"), no_overlap = paste0("a_", where_not_i, "_b_", where_not_j), a_to_b = sum(t_a < t_b), b_to_a = sum(t_a > t_b), eq_ab = sum(t_a == t_b))
              cnt <- cnt + 1
            }
          }
        }
        x_direct_df <- rbind(x_direct_df, x_direct)
        print(i)
      }
      finish <- 1
    }
    return(x_direct_df)
  } else {
    x_direct_df <- data.frame(matrix(NA, nrow = 0, ncol = 6))
    colnames(x_direct_df) <- c("traj_a", "traj_b", "no_overlap", "a_to_b", "b_to_a", "eq_ab")
    # get overlapping trajs
    tested <- vector(mode = "list")
    finish <- 0
    registerDoParallel(cl <- makeCluster(n_cpu))
    while(finish != 1) {
      x_direct_df <- foreach(i = 1:nrow(x_sig)) %dopar% {
        cnt <- 1
        x_direct <- data.frame(matrix(NA, nrow = 0, ncol = 6))
        colnames(x_direct) <- c("traj_a", "traj_b", "no_overlap", "a_to_b", "b_to_a", "eq_ab")
        traj_i <- unlist(x_sig[i,])
        if(test_len == 1) {
          for(j in 1:nrow(x_sig)) {
            if(i == j) {
              next
            }
            traj_j <- unlist(x_sig[j,])
            is_in <- any(sapply(tested, function(x) all(x %in% c(traj_i, traj_j))))
            if(is_in) {
              next
            } else {
              tested <- append(tested, list(c(traj_i, traj_j)))
              # filter for traj 
              test_ind <- which(x_ind_list[[which(names(x_ind_list) == traj_i)]] & x_ind_list[[which(names(x_ind_list) == traj_j)]])
              traj_filt <- x_traj[test_ind,]
              dates_filt <- x_dates[test_ind,]
              if(length(test_ind) < min_size) {
                x_direct[cnt,] <- data.frame(traj_a = paste0(traj_i, collapse = "->"), traj_b = paste0(traj_j, collapse = "->"), no_overlap = NA, NA, NA, NA)
                next
              }
              ind_a <- apply(traj_filt, 1, function(x) which(x == traj_i))
              ind_b <- apply(traj_filt, 1, function(x) which(x == traj_j))
              t_a <- dates_filt[cbind(1:nrow(dates_filt), ind_a)]
              t_b <- dates_filt[cbind(1:nrow(dates_filt), ind_b)]
              x_direct[cnt,] <- data.frame(traj_a = paste0(traj_i, collapse = "->"), traj_b = paste0(traj_j, collapse = "->"), no_overlap = paste0("a_", 1, "_b_", 1), a_to_b = sum(t_a < t_b), b_to_a = sum(t_a > t_b), eq_ab = sum(t_a == t_b))
              cnt <- cnt + 1
            }
          }
        } else {
          ol_i <- apply(x_sig, 1, function(x) any(x %in% traj_i))
          ol_i <- which(ol_i)
          ol_i <- ol_i[ol_i != i]
          ol_use_i <- x_sig[ol_i, ]
          for(j in 1:nrow(ol_use_i)) {
            traj_j <- unlist(ol_use_i[j,])
            is_in <- any(sapply(tested, function(x) all(x %in% c(traj_i, traj_j))))
            if(is_in) {
              next
            } else {
              tested <- append(tested, list(c(traj_i, traj_j)))
              # where is overlap?
              where_ol_i <- which(traj_i %in% traj_j)
              where_ol_j <- which(traj_j %in% traj_i)
              if(length(where_ol_j) < test_len - 1 | all(where_ol_i != sort(where_ol_i))) {
                next
              }
              # where is non-overlapping?
              where_not_i <- which(!traj_i %in% traj_j)
              where_not_j <- which(!traj_j %in% traj_i)
              # filter for traj 
              test_ind <- which(apply(sapply(unique(c(traj_i, traj_j)), function(x) x_ind_list[[which(names(x_ind_list) == x)]]), 1, all))
              traj_filt <- x_traj[test_ind,]
              dates_filt <- x_dates[test_ind,]
              if(length(test_ind) < min_size) {
                x_direct[cnt,] <- data.frame(traj_a = paste0(traj_i, collapse = "->"), traj_b = paste0(traj_j, collapse = "->"), no_overlap = NA, NA, NA, NA)
                next
              }
              # order filter
              traj_ord <- traj_i[where_ol_i]
              traj_ord <- sapply(traj_ord, function(x) dates_filt[cbind(1:nrow(dates_filt), apply(traj_filt, 1, function(y) which(y == x)))])
              t_filt <- apply(traj_ord, 1, function(x) all(x == sort(x) & !duplicated(x)))
              if(sum(t_filt) < min_size) {
                x_direct[cnt,] <- data.frame(traj_a = paste0(traj_i, collapse = "->"), traj_b = paste0(traj_j, collapse = "->"), no_overlap = NA, NA, NA, NA)
                next
              }
              traj_filt <- traj_filt[t_filt,]
              dates_filt <- dates_filt[t_filt,]
              dis_a <- traj_i[where_not_i]
              dis_b <- traj_j[where_not_j]
              ind_a <- apply(traj_filt, 1, function(x) which(x == dis_a))
              ind_b <- apply(traj_filt, 1, function(x) which(x == dis_b))
              t_a <- dates_filt[cbind(1:nrow(dates_filt), ind_a)]
              t_b <- dates_filt[cbind(1:nrow(dates_filt), ind_b)]
              x_direct[cnt,] <- data.frame(traj_a = paste0(traj_i, collapse = "->"), traj_b = paste0(traj_j, collapse = "->"), no_overlap = paste0("a_", where_not_i, "_b_", where_not_j), a_to_b = sum(t_a < t_b), b_to_a = sum(t_a > t_b), eq_ab = sum(t_a == t_b))
              cnt <- cnt + 1
            }
          }
        }
        return(x_direct)
      }
      finish <- 1
    }
    x_direct_df <- rbindlist(x_direct_df)
    colnames(x_direct_df) <- c("traj_a", "traj_b", "no_overlap", "a_to_b", "b_to_a", "eq_ab")
    return(x_direct_df)
  }
}

get_dir_stats <- function(direct_df) {
  # remove NAs
  
  direct_df <- na.omit(direct_df)
  
  # calculate stats
  
  dir_stat <- data.frame(matrix(nrow = nrow(direct_df), ncol = 10))
  colnames(dir_stat) <- c("traj_a", "traj_b", "no_overlap", "a_to_b", "b_to_a", "eq_ab", "est_a_to_b", 
                          "est_b_to_a", "p_a_to_b", "p_b_to_a")
  
  for(i in 1:nrow(direct_df)) {
    dir_i <- unlist(direct_df[i,])
    t1 <- binom.test(as.numeric(c(dir_i[4], sum(as.numeric(dir_i[5:6])))), alternative = "greater")
    t2 <- binom.test(as.numeric(c(dir_i[5], sum(as.numeric(dir_i[c(4, 6)])))), alternative = "greater")
    p1 <- t1$p.value
    p2 <- t2$p.value
    est_1 <- t1$estimate
    est_2 <- t2$estimate
    dir_stat[i,] <- c(dir_i, est_1, est_2, p1, p2)
  }
  
  dir_stat %<>%
    mutate(a_to_b = as.numeric(a_to_b),
           b_to_a = as.numeric(b_to_a),
           eq_ab = as.numeric(eq_ab),
           est_a_to_b = as.numeric(est_a_to_b),
           est_b_to_a = as.numeric(est_b_to_a),
           p_a_to_b = as.numeric(p_a_to_b),
           p_b_to_a = as.numeric(p_b_to_a))
  
  # adjust p-values
  
  dir_stat %<>%
    mutate(adj_p_a_to_b = p.adjust(p_a_to_b, method = "holm"),
           adj_p_b_to_a = p.adjust(p_b_to_a, method = "holm"))
  
  # filter for significant
  
  dir_stat %<>%
    filter(adj_p_a_to_b < 0.05 | adj_p_b_to_a < 0.05) %>%
    mutate(adj_p_value = ifelse(adj_p_a_to_b < 0.05, adj_p_a_to_b, adj_p_b_to_a)) %>%
    select(-p_a_to_b, -p_b_to_a, -adj_p_a_to_b, -adj_p_b_to_a) %>%
    select(traj_a, traj_b, no_overlap, a_to_b, b_to_a, adj_p_value)
  
  return(dir_stat)
}


get_direct <- Vectorize(function(no_overlap, traj_a, traj_b, a_to_b, b_to_a, test_len) {
  no_overlap <- as.numeric(str_split(no_overlap, "_")[[1]][c(2, 4)])
  traj_a <- str_split(traj_a, "->")[[1]]
  traj_b <- str_split(traj_b, "->")[[1]]
  traj_ol <- traj_a[-no_overlap[1]]
  dis_a <- traj_a[no_overlap[1]]
  dis_b <- traj_b[no_overlap[2]]
  if(a_to_b < b_to_a) {
    traj_o <- rep(NA, test_len + 1)
    traj_o[no_overlap[2]] <- dis_b
    traj_o[no_overlap[1]+1] <- dis_a
    traj_o[is.na(traj_o)] <- traj_ol
    traj_o <- paste0(traj_o, collapse = "->")
    return(traj_o)
  } else {
    traj_o <- rep(NA, test_len+1)
    traj_o[no_overlap[1]] <- dis_a
    traj_o[no_overlap[2]+1] <- dis_b
    traj_o[is.na(traj_o)] <- traj_ol
    traj_o <- paste0(traj_o, collapse = "->")
    return(traj_o)
  }
})


# read in significant pairwise directionalities

fem_sig_dir <- fread(fem_sig_dir_file)
male_sig_dir <- fread(male_sig_dir_file)

# read in trajectories

fem_traj <- fread(fem_traj_file)
male_traj <- fread(male_traj_file)

# read in dates

fem_dates <- fread(fem_dates_file)
male_dates <- fread(male_dates_file)

# get morbidity

fem_morb <- apply(fem_traj, 1, function(x) sum(!is.na(x)))
male_morb <- apply(male_traj, 1, function(x) sum(!is.na(x)))

# morbidity filter for > 1

fem_traj <- fem_traj[fem_morb > 1,]
male_traj <- male_traj[male_morb > 1,]

fem_dates <- fem_dates[fem_morb > 1,]
male_dates <- male_dates[male_morb > 1,]

# convert to matrix

fem_traj <- as.matrix(fem_traj)
male_traj <- as.matrix(male_traj)

fem_dates <- as.matrix(fem_dates)
male_dates <- as.matrix(male_dates)

# pre-compute indices

fem_ind_list <- lapply(unique(c(pull(fem_sig_dir, 1), pull(fem_sig_dir, 2))), function(x) apply(fem_traj, 1, function(y) x %in% y))
names(fem_ind_list) <- unique(c(pull(fem_sig_dir, 1), pull(fem_sig_dir, 2)))

male_ind_list <- lapply(unique(c(pull(male_sig_dir, 1), pull(male_sig_dir, 2))), function(x) apply(male_traj, 1, function(y) x %in% y))
names(male_ind_list) <- unique(c(pull(male_sig_dir, 1), pull(male_sig_dir, 2)))

# get triplet counts

fem_trip_cnts <- get_traj_counts(fem_sig_dir[, 1:2], fem_traj, fem_dates, fem_ind_list, 2, 20)
male_trip_cnts <- get_traj_counts(male_sig_dir[, 1:2], male_traj, male_dates, male_ind_list, 2, 20)

# get triplet stats

fem_trip_stats <- get_dir_stats(fem_trip_cnts)
male_trip_stats <- get_dir_stats(male_trip_cnts)

# get directionality

fem_trip_stats %<>%
  mutate(direct = get_direct(no_overlap, traj_a, traj_b, a_to_b, b_to_a, 2))

male_trip_stats %<>%
  mutate(direct = get_direct(no_overlap, traj_a, traj_b, a_to_b, b_to_a, 2))

# split directions

fem_trip_stats %<>%
  separate(direct, into = c("traj_a", "traj_b", "traj_c"), sep = "->") %>%
  select(traj_a, traj_b, traj_c, no_overlap, a_to_b, b_to_a, adj_p_value)

male_trip_stats %<>%
  separate(direct, into = c("traj_a", "traj_b", "traj_c"), sep = "->") %>%
  select(traj_a, traj_b, traj_c, no_overlap, a_to_b, b_to_a, adj_p_value)

# get quartet counts

fem_quart_cnts <- get_traj_counts(fem_trip_stats[, 1:3], fem_traj, fem_dates, fem_ind_list, 3, 20)
male_quart_cnts <- get_traj_counts(male_trip_stats[, 1:3], male_traj, male_dates, male_ind_list, 3, 20)

# get quartet stats

fem_quart_stats <- get_dir_stats(fem_quart_cnts)
male_quart_stats <- get_dir_stats(male_quart_cnts)

# get directionality

fem_quart_stats %<>%
  mutate(direct = get_direct(no_overlap, traj_a, traj_b, a_to_b, b_to_a, 3))

male_quart_stats %<>%
  mutate(direct = get_direct(no_overlap, traj_a, traj_b, a_to_b, b_to_a, 3))

# split directions

fem_quart_stats %<>%
  separate(direct, into = c("traj_a", "traj_b", "traj_c", "traj_d"), sep = "->") %>%
  select(traj_a, traj_b, traj_c, traj_d, no_overlap, a_to_b, b_to_a, adj_p_value)

male_quart_stats %<>%
  separate(direct, into = c("traj_a", "traj_b", "traj_c", "traj_d"), sep = "->") %>%
  select(traj_a, traj_b, traj_c, traj_d, no_overlap, a_to_b, b_to_a, adj_p_value)

# get quintet counts

fem_quint_cnts <- get_traj_counts(fem_quart_stats[, 1:4], fem_traj, fem_dates, fem_ind_list, 4, 20)
male_quint_cnts <- get_traj_counts(male_quart_stats[, 1:4], male_traj, male_dates, male_ind_list, 4, 20)

# get quintet stats

fem_quint_stats <- get_dir_stats(fem_quint_cnts)
# no more
male_quint_stats <- get_dir_stats(male_quint_cnts)

male_quint_stats %<>%
  mutate(direct = get_direct(no_overlap, traj_a, traj_b, a_to_b, b_to_a, 4))

male_quint_stats %<>%
  separate(direct, into = c("traj_a", "traj_b", "traj_c", "traj_d", "traj_e"), sep = "->") %>%
  select(traj_a, traj_b, traj_c, traj_d, traj_e, no_overlap, a_to_b, b_to_a, adj_p_value)


# combine all trajectories

fem_sig_pair <- fem_sig_dir %>%
  unite("traj", dis_a, dis_b, sep = "->") %>%
  transmute(traj = traj, p_value = adj_p_val_dir)

male_sig_pair <- male_sig_dir %>%
  unite("traj", dis_a, dis_b, sep = "->") %>%
  transmute(traj = traj, p_value = adj_p_val_dir)

fem_sig_trip <- fem_trip_stats %>%
  unite("traj", traj_a, traj_b, traj_c, sep = "->") %>%
  transmute(traj = traj, p_value = adj_p_value)

male_sig_trip <- male_trip_stats %>%
  unite("traj", traj_a, traj_b, traj_c, sep = "->") %>%
  transmute(traj = traj, p_value = adj_p_value)

fem_sig_quart <- fem_quart_stats %>%
  unite("traj", traj_a, traj_b, traj_c, traj_d, sep = "->") %>%
  transmute(traj = traj, p_value = adj_p_value)

male_sig_quart <- male_quart_stats %>%
  unite("traj", traj_a, traj_b, traj_c, traj_d, sep = "->") %>%
  transmute(traj = traj, p_value = adj_p_value)

male_sig_quint <- male_quint_stats %>%
  unite("traj", traj_a, traj_b, traj_c, traj_d, traj_e, sep = "->") %>%
  transmute(traj = traj, p_value = adj_p_value)


fem_sig_traj <- rbind(fem_sig_pair,
                      fem_sig_trip,
                      fem_sig_quart)

male_sig_traj <- rbind(male_sig_pair,
                      male_sig_trip,
                      male_sig_quart,
                      male_sig_quint)

# set column names

colnames(fem_sig_traj) <- c("direct", "adj_p_value")
colnames(male_sig_traj) <- c("direct", "adj_p_value")

# write out

write_tsv(fem_sig_traj, "fem_sig_trajs.tsv")
write_tsv(male_sig_traj, "male_sig_trajs.tsv")
