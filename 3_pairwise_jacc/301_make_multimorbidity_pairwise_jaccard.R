# creates pairwise jaccard coincidence matrix

# Load libraries

library(tidyverse)
library(jaccard)

# Variables

# female prevalence filtered multimorbidity matrix
fem_mm_file = "ukb_fem_multimorbidity_ICD10_physio_block_1pct_table.tsv"

# male prevalence filtered multimorbidity matrix
male_mm_file = "ukb_male_multimorbidity_ICD10_physio_block_1pct_table.tsv"

# output files

# female pairwise Jaccard coincidence matrix
fem_mm_jacc_file = "ukb_fem_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_table.tsv"

# male pairwise Jaccard coincidence matrix
male_mm_jacc_file = "ukb_male_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_table.tsv"

# Functions

pair_col_jacc <- function(x) {
  pairwise_matrix <- matrix(nrow = ncol(x), ncol = ncol(x))
  colnames(pairwise_matrix) <- colnames(x)
  rownames(pairwise_matrix) <- colnames(x)
  for(i in 1:ncol(x)) {
    for(j in 1:ncol(x)) {
      jacc_ij <- jaccard(x[,i], x[,j])
      pairwise_matrix[i,j] <- jacc_ij
    }
    print(i)
  }
  return(pairwise_matrix)
}

# Load data

fem_mm <- as.matrix(read_tsv(fem_mm_file))

male_mm <- as.matrix(read_tsv(male_mm_file))

# morbidity filter

fem_morb <- rowSums(fem_mm)
male_morb <- rowSums(male_mm)

fem_mm <- fem_mm[which(fem_morb >= 2),]
male_mm <- male_mm[which(male_morb >= 2),]

# make pairwise column jaccard comparison matrix of multimorbidity matrix

fem_jacc <- pair_col_jacc(fem_mm)

male_jacc <- pair_col_jacc(male_mm)

# write out pairwise jaccard matrix

write_tsv(as.data.frame(fem_jacc), fem_mm_jacc_file)

write_tsv(as.data.frame(male_jacc), male_mm_jacc_file)

