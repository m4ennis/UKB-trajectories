# run relative time clustering with optimal parameters

library(tidyverse)
library(data.table)
library(FactoMineR)
library(cluster)
library(fpc)
library(magrittr)
library(ComplexHeatmap)
library(parallel)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidymodels)

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

# female time of follow-up end relative to first diagnosis matrix
fem_norm_fu_file = "fem_norm_follow_up.tsv"
# male time of follow-up end relative to first diagnosis matrix
male_norm_fu_file = "male_norm_follow_up.tsv"

# file containing sex of participants
sex_file = ""

# file containing date of death of participants
death_file = ""

# ICD10 block dictionary
icd10_blocks_file = "icd10_blocks.tsv"

# ICD10 all chapters matrix

all_chapt_file = "ukb_multimorbidity_ICD10_all_chapter_table.tsv"

# unfiltered block multimorbidity matrix

block_mm_file = "ukb_multimorbidity_ICD10_block_table.tsv"

# read in trajectories

fem_traj <- fread(fem_traj_file)
male_traj <- fread(male_traj_file)

fem_dates <- fread(fem_dates_file)
male_dates <- fread(male_dates_file)

# read in blocks dictionary

icd10_blocks <- fread(icd10_blocks_file)

# read in death

death <- fread(death_file)
death <- death$date_of_death_f40000_0_0

# filter for sex

sex <- fread(sex_file)

sex <- sex$sex_f31_0_0

# filter death for sex

fem_death <- death[sex == 0]
male_death <- death[sex == 1]

# read in unnormalised trajectory dates

fem_unnorm <- fread(fem_unnorm_dates_file)
male_unnorm <- fread(male_unnorm_dates_file)

# read in follow-up

fem_fu <- fread(fem_norm_fu_file) %>% pull(1)
male_fu <- fread(male_norm_fu_file) %>% pull(1)

# get date of first diagnosis

fem_first <- fem_unnorm %>% pull(1)
male_first <- male_unnorm %>% pull(1)

# get morb

fem_morb <- apply(fem_traj, 1, function(x) sum(!is.na(x)))
male_morb <- apply(male_traj, 1, function(x) sum(!is.na(x)))

# filter for morbidity > 1

fem_traj <- fem_traj[fem_morb > 1,]
male_traj <- male_traj[male_morb > 1,]

fem_dates <- fem_dates[fem_morb > 1,]
male_dates <- male_dates[male_morb > 1,]

fem_death <- fem_death[fem_morb > 1]
male_death <- male_death[male_morb > 1]

fem_first <- fem_first[fem_morb > 1]
male_first <- male_first[male_morb > 1]

fem_fu <- fem_fu[fem_morb > 1]
male_fu <- male_fu[male_morb > 1]

# convert death date to relative time

fem_death <- fem_death - fem_first
male_death <- male_death - male_first

# get quartiles

fem_qnts <- sapply(c(0.25, 0.5, 0.75, 0.9), function(x) quantile(unlist(fem_dates)[!is.na(unlist(fem_dates))], x))
male_qnts <- sapply(c(0.25, 0.5, 0.75, 0.9), function(x) quantile(unlist(male_dates)[!is.na(unlist(male_dates))], x))

# get death labels

fem_deaths <- sapply(fem_qnts, function(x) fem_death < x)
male_deaths <- sapply(male_qnts, function(x) male_death < x)

# convert NAs to false

fem_deaths[is.na(fem_deaths)] <- F
male_deaths[is.na(male_deaths)] <- F

# get follow-up labels

fem_fu <- sapply(fem_qnts, function(x) fem_fu < x)
male_fu <- sapply(male_qnts, function(x) male_fu < x)

# convert to vector

fem_deaths <- c(fem_deaths[, 1], fem_deaths[, 2], fem_deaths[, 3], fem_deaths[, 4])
male_deaths <- c(male_deaths[, 1], male_deaths[, 2], male_deaths[, 3], male_deaths[, 4])

fem_fu <- c(fem_fu[, 1], fem_fu[, 2], fem_fu[, 3], fem_fu[, 4])
male_fu <- c(male_fu[, 1], male_fu[, 2], male_fu[, 3], male_fu[, 4])

# convert to data.frame

fem_dates <- as.data.frame(fem_dates)
male_dates <- as.data.frame(male_dates)

# create filtered sets of trajs and dates

fem_q_trajs <- lapply(1:4, function(x) as.matrix(fem_traj))
for(i in 1:length(fem_qnts)) {
  for(j in 1:ncol(fem_dates)) {
    fem_q_trajs[[i]][, j][fem_dates[, j] > fem_qnts[i]] <- NA
  }
}

male_q_trajs <- lapply(1:4, function(x) as.matrix(male_traj))
for(i in 1:length(male_qnts)) {
  for(j in 1:ncol(male_dates)) {
    male_q_trajs[[i]][, j][male_dates[, j] > male_qnts[i]] <- NA
  }
}

# convert to multimorbidity matrices

fem_dis <- unique(as.vector(fem_q_trajs[[4]]))
fem_dis <- sort(fem_dis[!is.na(fem_dis)])

male_dis <- unique(as.vector(male_q_trajs[[4]]))
male_dis <- sort(male_dis[!is.na(male_dis)])

fem_mm <- lapply(fem_q_trajs, function(x) {
  mm_x <- matrix(0, nrow = nrow(x), ncol = length(fem_dis))
  colnames(mm_x) <- fem_dis
  for(i in 1:nrow(x)) {
    ind_i <- which(colnames(mm_x) %in% x[i,])
    mm_x[i, ind_i] <- 1
  }
  return(mm_x)
})

male_mm <- lapply(male_q_trajs, function(x) {
  mm_x <- matrix(0, nrow = nrow(x), ncol = length(male_dis))
  colnames(mm_x) <- male_dis
  for(i in 1:nrow(x)) {
    ind_i <- which(colnames(mm_x) %in% x[i,])
    mm_x[i, ind_i] <- 1
  }
  return(mm_x)
})

# add quantile label

for(i in 1:length(fem_qnts)) {
  # colnames(fem_mm[[i]]) <- paste0(colnames(fem_mm[[i]]), "_", round(fem_qnts[i], 1), "_yrs")
  fem_mm[[i]] <- cbind(fem_mm[[i]], i)
}

for(i in 1:length(male_qnts)) {
  # colnames(male_mm[[i]]) <- paste0(colnames(male_mm[[i]]), "_", round(male_qnts[i], 1), "_yrs")
  male_mm[[i]] <- cbind(male_mm[[i]], i)
}

# reduce mm to single mm

fem_mm <- reduce(fem_mm, rbind)
male_mm <- reduce(male_mm, rbind)

# get quantile label

fem_lab <- fem_mm[, ncol(fem_mm)]
male_lab <- male_mm[, ncol(male_mm)]

mm_lab <- c(fem_lab, male_lab)

# get sex label

sex_lab <- c(rep("f", nrow(fem_mm)), rep("m", nrow(male_mm)))

# remove label

fem_mm <- fem_mm[, -ncol(fem_mm)]
male_mm <- male_mm[, -ncol(male_mm)]

# only look at diseases common to both

dis_all <- intersect(colnames(fem_mm), colnames(male_mm))

fem_mm <- fem_mm[, which(colnames(fem_mm) %in% dis_all)]
male_mm <- male_mm[, which(colnames(male_mm) %in% dis_all)]

# combine

mm <- rbind(fem_mm, male_mm)

# remove those who have died or lost follow-up

fm_deaths <- c(fem_deaths, male_deaths)
fm_fu <- c(fem_fu, male_fu)

mm <- mm[!fm_fu & !fm_deaths,]

# convert to logical

mm <- apply(mm, 2, as.logical)

# apply MCA to the mm

mm_mca <- MCA(mm)

# plot screes

plot <- data.frame(eigenvalue = mm_mca$eig[, 1], comp = 1:length(mm_mca$eig[, 1])) %>%
  ggplot(aes(x = comp, y = eigenvalue)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, length(mm_mca$eig[, 1]), 3)) +
  xlab("MCA Component") +
  ylab("Eigenvalue")
ggsave("mm_mca_scree.png", plot, units = "px", width = 1800, height = 1600)


# 5 for females
# 6 for males
# 5 for all

# fem_ncp <- 5
# male_ncp <- 6

mm_ncp <- 6

mm_mca <- MCA(mm, ncp = mm_ncp)

# get coordinates

mm_coor <- mm_mca$ind$coord

# cluster

k_max <- 15

set.seed(42)

mm_crit <- sapply(2:k_max, function(k) {
  print(k)
  cl <- kmeans(mm_coor, k, iter.max = 30, nstart = 100)
  rss <- cl$tot.withinss
  ch <- calinhara(mm_coor, cl$cluster)
  return(c(rss = rss, ch = ch))
})

plot <- data.frame(rss = mm_crit[1, ], k = 2:k_max) %>%
  ggplot(aes(x = k, y = rss)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, k_max, 2)) +
  xlab("Number of clusters") +
  ylab("Total within cluster sum of squares")
ggsave("mm_kmeans_k_rss.png", plot, units = "px", width = 1800, height = 1600)

plot <- data.frame(ch = mm_crit[2, ], k = 2:k_max) %>%
  ggplot(aes(x = k, y = ch)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, k_max, 2)) +
  xlab("Number of clusters") +
  ylab("Calinski-Harabasz score")
ggsave("mm_kmeans_k_ch.png", plot, units = "px", width = 1800, height = 1600)


# decide k

mm_k <- 6

set.seed(42)

mm_ks <- lapply(2:mm_k, function(k) kmeans(mm_coor, k, iter.max = 30, nstart = 100))
mm_k <- mm_ks[[length(mm_ks)]]
mm_cl <- mm_k$cluster

# get cluster sizes

mm_sizes <- data.frame(qnt = paste0("t_", mm_lab), sex = sex_lab, cl = mm_cl) %>%
  group_by(qnt, sex) %>%
  count(cl) %>%
  pivot_wider(names_from = qnt, values_from = n)

# add total

mm_sizes %<>%
  group_by(cl) %>%
  mutate(tot = sum(t_1, t_2, t_3, t_4))

# format column names

colnames(mm_sizes) <- c("Sex", "Cluster", "Size at T1", "Size at T2", "Size at T3", "Size at T4", "Total size")

# write out

write_tsv(mm_sizes, "mm_cl_sizes.tsv")

# add clusters

mm %<>%
  as.data.frame() %>%
  mutate(cl = mm_cl)

# get disease prevalences

mm_prev <- mm %>% group_by(cl) %>% summarise(across(.fns = function(x) sum(x)/length(cl)))

# get those with >15% in at least one cluster

mm_use <- mm_prev %>%
  pivot_longer(-cl, names_to = "block", values_to = "prev") %>%
  filter(prev > 0.15) %>%
  pull(block) %>%
  unique()

# write out table of top prevalences

mm_prev %>%
  pivot_longer(-cl, names_to = "block", values_to = "prev") %>%
  filter(prev > 0.15) %>%
  arrange(cl, desc(prev)) %>%
  mutate(name = icd10_blocks$name[match(block, icd10_blocks$block)]) %>%
  select(cl, block, name, prev) %>%
  write_tsv("mm_cl_top_dis_prev.tsv")


mm_mat <- mm_prev %>%
  select(mm_use) %>%
  as.matrix() %>%
  scale()

# format row and columns names

rownames(mm_mat) <- paste0("Cluster ", mm_prev$cl)
colnames(mm_mat) <- str_trunc(icd10_blocks$name[match(colnames(mm_mat), icd10_blocks$block)], 30)

# add in death labels

mm_deaths <- c(fem_deaths, male_deaths)

mm_state <- character(length = length(mm_lab))
mm_state[fm_fu] <- "End of follow-up"
mm_state[fm_deaths] <- "Death"
mm_state[!fm_fu & !fm_deaths] <- mm_cl
mm_cl <- mm_state

fem_snk <- data.frame(qnt_1 = paste0("qnt_1_cl_", mm_cl[sex_lab == "f"][fem_lab == 1]),
                      qnt_2 = paste0("qnt_2_cl_", mm_cl[sex_lab == "f"][fem_lab == 2]),
                      qnt_3 = paste0("qnt_3_cl_", mm_cl[sex_lab == "f"][fem_lab == 3]),
                      qnt_4 = paste0("qnt_4_cl_", mm_cl[sex_lab == "f"][fem_lab == 4]))

fem_links <- vector(mode = "list", length = 3)
for(i in 1:3) {
  fem_links[[i]] <- as.data.frame(table(fem_snk[, i:(i+1)]))
  colnames(fem_links[[i]]) <- c("source", "target", "value")
}

fem_links <- reduce(fem_links, rbind)

# write out table

write_tsv(fem_links, "fem_cl_flows.tsv")

fem_links %<>%
  filter(value > 0) %>%
  unite("source", -target, sep = " [") %>% unite("final", sep = "] ")

write_tsv(fem_links, "fem_cluster_sankey_flows.tsv")


male_snk <- data.frame(qnt_1 = paste0("qnt_1_cl_", mm_cl[sex_lab == "m"][male_lab == 1]),
                       qnt_2 = paste0("qnt_2_cl_", mm_cl[sex_lab == "m"][male_lab == 2]),
                       qnt_3 = paste0("qnt_3_cl_", mm_cl[sex_lab == "m"][male_lab == 3]),
                       qnt_4 = paste0("qnt_4_cl_", mm_cl[sex_lab == "m"][male_lab == 4]))

male_links <- vector(mode = "list", length = 3)
for(i in 1:3) {
  male_links[[i]] <- as.data.frame(table(male_snk[, i:(i+1)]))
  colnames(male_links[[i]]) <- c("source", "target", "value")
}

male_links <- reduce(male_links, rbind)

# write out table

write_tsv(male_links, "male_cl_flows.tsv")

male_links %<>%
  filter(value > 0) %>%
  unite("source", -target, sep = " [") %>% unite("final", sep = "] ")
  
write_tsv(male_links, "male_cluster_sankey_flows.tsv")
