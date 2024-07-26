# filter for physiological and prevalence >= 1%

print("Loading libraries...")

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(data.table)

# Variables

# female multimorbidity matrix
fem_mm_file = "ukb_fem_multimorbidity_ICD10_block_table.tsv"

# male multimorbidity matrix
male_mm_file = "ukb_male_multimorbidity_ICD10_block_table.tsv"

# output files

# female diseases passing prevalence filter
fem_disease_use_file = "fem_physio_disease_use_1pct.tsv"

# male diseases passing prevalence filter
male_disease_use_file = "male_physio_disease_use_1pct.tsv"

# female prevalence filtered multimorbidity matrix 
fem_mm_filt_file = "ukb_fem_multimorbidity_ICD10_physio_block_1pct_table.tsv"

# male prevalence filtered multimorbidity matrix 

male_mm_filt_file = "ukb_male_multimorbidity_ICD10_physio_block_1pct_table.tsv"

# read in multimorbidity matrix

fem_mm <- as.data.frame(fread(fem_mm_file))

male_mm <- as.data.frame(fread(male_mm_file))

# get ICD10 block prevalences

fem_prevs <- colSums(fem_mm)/nrow(fem_mm)

male_prevs <- colSums(male_mm)/nrow(male_mm)

# get blocks >= 1% prevalence

fem_prev_gt <- fem_prevs >= 0.01

male_prev_gt <- male_prevs >= 0.01

# get physiological blocks

fem_physio <- !str_starts(names(fem_prevs), "O|P|Q|R|S|T|U|V|W|X|Y|Z")

male_physio <- !str_starts(names(fem_prevs), "O|P|Q|R|S|T|U|V|W|X|Y|Z")

# get filter ind

fem_filt_ind <- which(fem_prev_gt & fem_physio)

male_filt_ind <- which(male_prev_gt & male_physio)

# get disease use

fem_disease_use <- names(fem_filt_ind)

male_disease_use <- names(male_filt_ind)

# filter multimorbidity matrices

fem_mm_filt <- fem_mm[, fem_filt_ind]

male_mm_filt <- male_mm[, male_filt_ind]

# write out diseases to use

write_tsv(as.data.frame(fem_disease_use), fem_disease_use_file)

write_tsv(as.data.frame(male_disease_use), male_disease_use_file)

# write out physiological and prevalence filtered multimorbidity matrices

write_tsv(fem_mm_filt, fem_mm_filt_file)

write_tsv(male_mm_filt, male_mm_filt_file)








