# runs tests of pairwise directionality 

library(tidyverse)
library(data.table)
library(magrittr)

# Variables

# female Jaccard coincidence statistical test statistics 
fem_stat_file = "ukb_fem_multimorbidity_ICD10_physio_block_1pct_stats.tsv"
# male Jaccard coincidence statistical test statistics
male_stat_file = "ukb_male_multimorbidity_ICD10_physio_block_1pct_stats.tsv"

# female Jaccard coincidence matrix
fem_jacc_file = "ukb_fem_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_table.tsv"
# male Jaccard coincidence matrix
male_jacc_file = "ukb_male_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_table.tsv"

# female diagnosis trajectory matrix
fem_traj_file = "fem_disease_block_1pct_traj.tsv"
# male diagnosis trajectory matrix
male_traj_file = "male_disease_block_1pct_traj.tsv"

# female time to diagnosis matrix
fem_dates_file = "fem_diag_block_1pct_norm_dates.tsv"
# male time to diagnosis matrix
male_dates_file = "male_diag_block_1pct_norm_dates.tsv"

# female prevalence filtered multimorbidity matrix
fem_mm_file = "ukb_fem_multimorbidity_ICD10_physio_block_1pct_table.tsv"
# male prevalence filtered multimorbidity matrix
male_mm_file = "ukb_male_multimorbidity_ICD10_physio_block_1pct_table.tsv"

# Load Jaccard stats

fem_stat <- fread(fem_stat_file)
male_stat <- fread(male_stat_file)

# filter for significant combinations

fem_sig <- fem_stat %>% filter(p_value < 0.05)
male_sig <- male_stat %>% filter(p_value < 0.05)

# filter for positive effects

fem_sig %<>%
  filter(stat > 0)

male_sig %<>%
  filter(stat > 0)

# read in jaccards

fem_jacc <- as.matrix(fread(fem_jacc_file))
male_jacc <- as.matrix(fread(male_jacc_file))

# set rownames

rownames(fem_jacc) <- colnames(fem_jacc)
rownames(male_jacc) <- colnames(male_jacc)

# read in trajectories

fem_traj <- fread(fem_traj_file)
male_traj <- fread(male_traj_file)

# morbidity filter

fem_morb <- apply(fem_traj, 1, function(x) sum(!is.na(x)))
male_morb <- apply(male_traj, 1, function(x) sum(!is.na(x)))

fem_traj <- fem_traj[which(fem_morb >= 2),]
male_traj <- male_traj[which(male_morb >= 2),]

# read in dates

fem_dates <- as.matrix(fread(fem_dates_file))
male_dates <- as.matrix(fread(male_dates_file))

# filter 

fem_dates <- fem_dates[fem_morb > 1,]
male_dates <- male_dates[male_morb > 1,]

# get directionality counts

# females

fem_direct <- fem_jacc %>%
  as.data.frame() %>%
  mutate(dis_a = rownames(fem_jacc)) %>%
  pivot_longer(-dis_a, names_to = "dis_b", values_to = "jacc") %>%
  select(-jacc) %>%
  mutate(a_to_b = NA, b_to_a = NA, eq_ab = NA) %>%
  as.data.frame()

fem_traj <- as.matrix(fem_traj)

# pre-create subsets

unique_dis <- unique(fem_direct$dis_a)

ind_list <- lapply(unique_dis, function(x) which(apply(fem_traj, 1, function(y) x %in% y)))
names(ind_list) <- unique_dis

for(i in 1:nrow(fem_direct)) {
  dis_a <- fem_direct$dis_a[i]
  dis_b <- fem_direct$dis_b[i]
  ind <- intersect(ind_list[[which(names(ind_list) == dis_a)]], ind_list[[which(names(ind_list) == dis_b)]])
  traj_filt <- fem_traj[ind,]
  dates_filt <- fem_dates[ind,]
  ind_a <- apply(traj_filt, 1, function(x) which(x == dis_a))
  ind_b <- apply(traj_filt, 1, function(x) which(x == dis_b))
  t_a <- dates_filt[cbind(1:nrow(dates_filt), ind_a)]
  t_b <- dates_filt[cbind(1:nrow(dates_filt), ind_b)]
  fem_direct[i, 3:5] <- c(sum(t_a < t_b), sum(t_a > t_b), sum(t_a == t_b))
}

# males

male_direct <- male_jacc %>%
  as.data.frame() %>%
  mutate(dis_a = rownames(male_jacc)) %>%
  pivot_longer(-dis_a, names_to = "dis_b", values_to = "jacc") %>%
  select(-jacc) %>%
  mutate(a_to_b = NA, b_to_a = NA, eq_ab = NA) %>%
  as.data.frame()

male_traj <- as.matrix(male_traj)

# pre-create subsets

unique_dis <- unique(male_direct$dis_a)

ind_list <- lapply(unique_dis, function(x) which(apply(male_traj, 1, function(y) x %in% y)))
names(ind_list) <- unique_dis

for(i in 1:nrow(male_direct)) {
  dis_a <- male_direct$dis_a[i]
  dis_b <- male_direct$dis_b[i]
  ind <- intersect(ind_list[[which(names(ind_list) == dis_a)]], ind_list[[which(names(ind_list) == dis_b)]])
  traj_filt <- male_traj[ind,]
  dates_filt <- male_dates[ind,]
  ind_a <- apply(traj_filt, 1, function(x) which(x == dis_a))
  ind_b <- apply(traj_filt, 1, function(x) which(x == dis_b))
  t_a <- dates_filt[cbind(1:nrow(dates_filt), ind_a)]
  t_b <- dates_filt[cbind(1:nrow(dates_filt), ind_b)]
  male_direct[i, 3:5] <- c(sum(t_a < t_b), sum(t_a > t_b), sum(t_a == t_b))
}

# checkpoint save

saveRDS(fem_direct, "fem_direct_jacc.rds")
saveRDS(male_direct, "male_direct_jacc.rds")

fem_direct <- readRDS("fem_direct_jacc.rds")
male_direct <- readRDS("male_direct_jacc.rds")

# run binomial tests to identify significantly directional combinations

# filter for significantly coincident disease combinations

fem_direct %<>%
  right_join(fem_sig, by = c("dis_a", "dis_b")) %>%
  select(dis_a, dis_b, a_to_b, b_to_a, eq_ab)

male_direct %<>%
  right_join(male_sig, by = c("dis_a", "dis_b")) %>%
  select(dis_a, dis_b, a_to_b, b_to_a, eq_ab)

# calculate stats

fem_dir_stat <- data.frame(matrix(nrow = nrow(fem_direct), ncol = 9))
colnames(fem_dir_stat) <- c("dis_a", "dis_b", "a_to_b", "b_to_a", "eq_ab", "est_a_to_b", 
                            "est_b_to_a", "p_a_to_b", "p_b_to_a")

for(i in 1:nrow(fem_direct)) {
  dir_i <- unlist(fem_direct[i,])
  t1 <- binom.test(as.numeric(c(dir_i[3], sum(as.numeric(dir_i[4:5])))), alternative = "greater")
  t2 <- binom.test(as.numeric(c(dir_i[4], sum(as.numeric(dir_i[c(3, 5)])))), alternative = "greater")
  p1 <- t1$p.value
  p2 <- t2$p.value
  est_1 <- t1$estimate
  est_2 <- t2$estimate
  # if(p1 > 0.05 & p2 > 0.05) {
  #   direct_i <- paste0("Not significant")
  # } else if(p1 < 0.05 & p2 > 0.05) {
  #   direct_i <- paste0(dir_i[1], " -> ", dir_i[2])
  # } else if(p1 > 0.05 & p2 < 0.05) {
  #   direct_i <- paste0(dir_i[2], " -> ", dir_i[1])
  # }
  fem_dir_stat[i,] <- c(dir_i, est_1, est_2, p1, p2)
  print(i)
}

male_dir_stat <- data.frame(matrix(nrow = nrow(male_direct), ncol = 9))
colnames(male_dir_stat) <- c("dis_a", "dis_b", "a_to_b", "b_to_a", "eq_ab", "est_a_to_b", 
                             "est_b_to_a", "p_a_to_b", "p_b_to_a")

for(i in 1:nrow(male_direct)) {
  dir_i <- unlist(male_direct[i,])
  t1 <- binom.test(as.numeric(c(dir_i[3], sum(as.numeric(dir_i[4:5])))), alternative = "greater")
  t2 <- binom.test(as.numeric(c(dir_i[4], sum(as.numeric(dir_i[c(3, 5)])))), alternative = "greater")
  p1 <- t1$p.value
  p2 <- t2$p.value
  est_1 <- t1$estimate
  est_2 <- t2$estimate
  # if(p1 > 0.05 & p2 > 0.05) {
  #   direct_i <- paste0("Not significant")
  # } else if(p1 < 0.05 & p2 > 0.05) {
  #   direct_i <- paste0(dir_i[1], " -> ", dir_i[2])
  # } else if(p1 > 0.05 & p2 < 0.05) {
  #   direct_i <- paste0(dir_i[2], " -> ", dir_i[1])
  # }
  male_dir_stat[i,] <- c(dir_i, est_1, est_2, p1, p2)
}

# adjust p-values

fem_dir_stat %<>%
  mutate(adj_p_a_to_b = p.adjust(p_a_to_b, method = "BH"),
         adj_p_b_to_a = p.adjust(p_b_to_a, method = "BH"))

male_dir_stat %<>%
  mutate(adj_p_a_to_b = p.adjust(p_a_to_b, method = "BH"),
         adj_p_b_to_a = p.adjust(p_b_to_a, method = "BH"))

# filter for adjusted p-values

fem_sig_dir <- fem_dir_stat %>%
  filter(adj_p_a_to_b < 0.05 | adj_p_b_to_a < 0.05) %>%
  mutate(direct = ifelse(adj_p_a_to_b < 0.05, paste0(dis_a, " -> ", dis_b), paste0(dis_b, " -> ", dis_a))) %>%
  mutate(adj_p_value = ifelse(adj_p_a_to_b < 0.05, adj_p_a_to_b, adj_p_b_to_a)) %>%
  select(-adj_p_a_to_b, -adj_p_b_to_a, -p_a_to_b, -p_b_to_a) %>%
  select(dis_a, dis_b, a_to_b, b_to_a, direct, adj_p_value)

male_sig_dir <- male_dir_stat %>%
  filter(adj_p_a_to_b < 0.05 | adj_p_b_to_a < 0.05) %>%
  mutate(direct = ifelse(adj_p_a_to_b < 0.05, paste0(dis_a, " -> ", dis_b), paste0(dis_b, " -> ", dis_a))) %>%
  mutate(adj_p_value = ifelse(adj_p_a_to_b < 0.05, adj_p_a_to_b, adj_p_b_to_a)) %>%
  select(-adj_p_a_to_b, -adj_p_b_to_a, -p_a_to_b, -p_b_to_a) %>%
  select(dis_a, dis_b, a_to_b, b_to_a, direct, adj_p_value)

# filter for adjusted p-values

fem_sig_dir <- fem_dir_stat %>%
  filter(adj_p_a_to_b < 0.05 | adj_p_b_to_a < 0.05) %>%
  mutate(direct = ifelse(adj_p_a_to_b < 0.05, paste0(dis_a, " -> ", dis_b), paste0(dis_b, " -> ", dis_a))) %>%
  mutate(adj_p_value = ifelse(adj_p_a_to_b < 0.05, adj_p_a_to_b, adj_p_b_to_a)) %>%
  select(-adj_p_a_to_b, -adj_p_b_to_a, -p_a_to_b, -p_b_to_a) %>%
  select(dis_a, dis_b, a_to_b, b_to_a, direct, adj_p_value)

male_sig_dir <- male_dir_stat %>%
  filter(adj_p_a_to_b < 0.05 | adj_p_b_to_a < 0.05) %>%
  mutate(direct = ifelse(adj_p_a_to_b < 0.05, paste0(dis_a, " -> ", dis_b), paste0(dis_b, " -> ", dis_a))) %>%
  mutate(adj_p_value = ifelse(adj_p_a_to_b < 0.05, adj_p_a_to_b, adj_p_b_to_a)) %>%
  select(-adj_p_a_to_b, -adj_p_b_to_a, -p_a_to_b, -p_b_to_a) %>%
  select(dis_a, dis_b, a_to_b, b_to_a, direct, adj_p_value)

# join with jaccards

fem_sig_dir %<>%
  full_join(fem_sig, by = c("dis_a", "dis_b")) %>%
  filter(!is.na(p_value) & !is.na(adj_p_value)) %>%
  transmute(dis_a, dis_b, p_val_jacc = p_value, adj_p_val_dir = adj_p_value, direct) %>%
  arrange(p_val_jacc, adj_p_val_dir)

male_sig_dir %<>%
  full_join(male_sig, by = c("dis_a", "dis_b")) %>%
  filter(!is.na(p_value) & !is.na(adj_p_value)) %>%
  transmute(dis_a, dis_b, p_val_jacc = p_value, adj_p_val_dir = adj_p_value, direct) %>%
  arrange(p_val_jacc, adj_p_val_dir)


fem_sig_dir <- fem_jacc %>%
  as.data.frame() %>%
  mutate(dis_a = rownames(fem_jacc)) %>%
  pivot_longer(-dis_a, names_to = "dis_b", values_to = "jacc") %>%
  right_join(fem_sig_dir, by = c("dis_a", "dis_b"))

male_sig_dir <- male_jacc %>%
  as.data.frame() %>%
  mutate(dis_a = rownames(male_jacc)) %>%
  pivot_longer(-dis_a, names_to = "dis_b", values_to = "jacc") %>%
  right_join(male_sig_dir, by = c("dis_a", "dis_b"))

# convert to significant order

fem_sig_dir %<>%
  select(-dis_a, -dis_b) %>%
  separate(direct, into = c("dis_a", "dis_b"), sep = "->") %>%
  select(dis_a, dis_b, everything())

male_sig_dir %<>%
  select(-dis_a, -dis_b) %>%
  separate(direct, into = c("dis_a", "dis_b"), sep = "->") %>%
  select(dis_a, dis_b, everything())

write_tsv(fem_sig_dir, "fem_sig_pairwise_directional_diseases.tsv")
write_tsv(male_sig_dir, "male_sig_pairwise_directional_diseases.tsv")
