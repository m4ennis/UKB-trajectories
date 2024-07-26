# run Jaccard statistical tests of coincidence in parallel

library(tidyverse)
library(jaccard)
library(parallel)
library(foreach)
library(doParallel)

# Variables

# female prevalence filtered multimorbidity matrix
fem_mm_file = "ukb_fem_multimorbidity_ICD10_physio_block_1pct_table.tsv"
# male prevalence filtered multimorbidity matrix
male_mm_file = "ukb_male_multimorbidity_ICD10_physio_block_1pct_table.tsv"

# female Jaccard statistical test table
fem_stat_file = "ukb_fem_multimorbidity_ICD10_physio_block_1pct_stats.tsv"
# male Jaccard statistical test table
male_stat_file = "ukb_male_multimorbidity_ICD10_physio_block_1pct_stats.tsv"

# number of cores to use in parallel
n_cores <- 32

# set seed for replicability
set.seed(42)

# Load data

fem_mm <- as.matrix(read_tsv(fem_mm_file))

male_mm <- as.matrix(read_tsv(male_mm_file))

# morbidity filter

fem_morb <- rowSums(fem_mm)
male_morb <- rowSums(male_mm)

fem_mm <- fem_mm[which(fem_morb >= 2),]
male_mm <- male_mm[which(male_morb >= 2),]

# run statistical tests

cl <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cl)

# females

fem_stats <- foreach(i = 1:ncol(fem_mm), .combine = "rbind", .packages = c("jaccard")) %dopar% {
  stats_i <- data.frame(matrix(nrow = 0, ncol = 5))
  for(j in i:ncol(fem_mm)) {
    if(j <= i) {
      next
    } else {
      set.seed(42)
      stats_ij <- jaccard.test(fem_mm[,i], fem_mm[,j], method = "bootstrap", verbose = F)
      stats_ij <- data.frame(dis_a = colnames(fem_mm)[i], dis_b = colnames(fem_mm)[j], expected = stats_ij$expectation, stat = stats_ij$statistics, p_value = stats_ij$pvalue)
      stats_i <- rbind(stats_i, stats_ij)
    }
  }
  return(stats_i)
}

# males

male_stats <- foreach(i = 1:ncol(male_mm), .combine = "rbind", .packages = c("jaccard")) %dopar% {
  stats_i <- data.frame(matrix(nrow = 0, ncol = 5))
  for(j in i:ncol(male_mm)) {
    if(j == i) {
      next
    } else {
      set.seed(42)
      stats_ij <- jaccard.test(male_mm[,i], male_mm[,j], method = "bootstrap", verbose = F)
      stats_ij <- data.frame(dis_a = colnames(male_mm)[i], dis_b = colnames(male_mm)[j], expected = stats_ij$expectation, stat = stats_ij$statistics, p_value = stats_ij$pvalue)
      stats_i <- rbind(stats_i, stats_ij)
    }
  }
  return(stats_i)
}

# write out

write_tsv(fem_stats, fem_stat_file)
write_tsv(male_stats, male_stat_file)


