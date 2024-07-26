# run cluster solutions in parallel and decide parameters for relative time clustering

library(tidyverse)
library(data.table)
library(FactoMineR)
library(cluster)
library(fpc)
library(magrittr)
library(parallel)
library(foreach)
library(doParallel)

# female diagnosis trajectory matrix
fem_traj_file = "fem_disease_block_1pct_traj.tsv"
# male diagnosis trajectory matrix
male_traj_file = "male_disease_block_1pct_traj.tsv"

# female time to diagnosis matrix
fem_dates_file = "fem_diag_block_1pct_norm_dates.tsv"
# male time to diagnosis matrix
male_dates_file = "male_diag_block_1pct_norm_dates.tsv"

# ICD10 block dictionary
icd10_blocks_file = "icd10_blocks.tsv"

# number of CPUs to use in parallel
n_cpu = 1

# read in trajectories

fem_traj <- fread(fem_traj_file)
male_traj <- fread(male_traj_file)

fem_dates <- fread(fem_dates_file)
male_dates <- fread(male_dates_file)

# read in blocks dictionary

icd10_blocks <- fread(icd10_blocks_file)

# get morb

fem_morb <- apply(fem_traj, 1, function(x) sum(!is.na(x)))
male_morb <- apply(male_traj, 1, function(x) sum(!is.na(x)))

# filter for morbidity > 1

fem_traj <- fem_traj[fem_morb > 1,]
male_traj <- male_traj[male_morb > 1,]

fem_dates <- fem_dates[fem_morb > 1,]
male_dates <- male_dates[male_morb > 1,]

# get quartiles

fem_qnts <- sapply(1:4, function(x) quantile(unlist(fem_dates)[!is.na(unlist(fem_dates))], x/4))
male_qnts <- sapply(1:4, function(x) quantile(unlist(male_dates)[!is.na(unlist(male_dates))], x/4))

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

# add sex-specific diseases to other gender
# 
# fem_add <- colnames(male_mm)[which(!colnames(male_mm) %in% colnames(fem_mm))]
# male_add <- colnames(fem_mm)[which(!colnames(fem_mm) %in% colnames(male_mm))]
# 
# fem_mm <- cbind(fem_mm, matrix(0, nrow = nrow(fem_mm), ncol = length(fem_add)))
# colnames(fem_mm)[which(colnames(fem_mm) == "")] <- fem_add
# 
# male_mm <- cbind(male_mm, matrix(0, nrow = nrow(male_mm), ncol = length(male_add)))
# colnames(male_mm)[which(colnames(male_mm) == "")] <- male_add
# 
# fem_mm <- fem_mm[, order(colnames(fem_mm))]
# male_mm <- male_mm[, order(colnames(male_mm))]

# only look at diseases common to both

dis_all <- intersect(colnames(fem_mm), colnames(male_mm))

fem_mm <- fem_mm[, which(colnames(fem_mm) %in% dis_all)]
male_mm <- male_mm[, which(colnames(male_mm) %in% dis_all)]

# combine

mm <- rbind(fem_mm, male_mm)

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


# 6 for all

mm_ncp <- 6

# run MCA again

mm_mca <- MCA(mm, ncp = mm_ncp)

# get coordinates

mm_coor <- mm_mca$ind$coord

# cluster

k_max <- 15

registerDoParallel(n_cpu)

mm_scores <- foreach(k = 2:k_max, .export = c("mm_coor"), .combine = "rbind") %dopar% {
  set.seed(42)
  cl <- kmeans(mm_coor, k, iter.max = 30, nstart = 100)
  ch <- fpc::calinhara(mm_coor, cl$cluster)
  return(c(rss = cl$tot.withinss, ch = ch))
}

plot <- data.frame(rss = mm_scores[, 1], k = 2:k_max) %>%
  ggplot(aes(x = k, y = rss)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, k_max, 2)) +
  xlab("Number of clusters") +
  ylab("Total within cluster sum of squares")
ggsave("mm_kmeans_k_rss.png", plot, units = "px", width = 1800, height = 1600)

# mm_ch <- foreach(k = 2:k_max, .export = c("mm_coor"), .combine = "c") %dopar% {
#   set.seed(42)
#   return(fpc::calinhara(mm_coor, kmeans(mm_coor, k, iter.max = 30, nstart = 500)$cluster))
# }

plot <- data.frame(ch = mm_scores[, 2], k = 2:k_max) %>%
  ggplot(aes(x = k, y = ch)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, k_max, 2)) +
  xlab("Number of clusters") +
  ylab("Calinski-Harabasz score")
ggsave("mm_kmeans_k_ch.png", plot, units = "px", width = 1800, height = 1600)
