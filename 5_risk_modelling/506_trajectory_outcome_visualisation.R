# visualisation of trajectories and trajectory outcomes

library(tidyverse)
library(data.table)
library(magrittr)
library(data.tree)
library(igraph)
library(circlize)
library(ggrepel)
library(DiagrammeRsvg)
library(DiagrammeR)
library(UpSetR)
library(ggforce)
library(rpart)
library(rpart.plot)
library(ComplexHeatmap)

# variables

# female significant built trajectories with info
fem_sig_file = "female_significant_trajectory_info.tsv"
# male significant built trajectories with info
male_sig_file = "male_significant_trajectory_info.tsv"

# female trajectory mortality rate coefficients
fem_mort_file = "female_trajectory_risk_of_death_coefficients_new.tsv"

# female prior morbidity mortality rate coefficients
fem_morb_mort_file = "female_prior_morbidity_risk_of_death_coefficients_new.tsv"

# male trajectory mortality rate coefficients
male_mort_file = "male_trajectory_risk_of_death_coefficients_new.tsv"

# male prior morbidity mortality rate coefficients
male_morb_mort_file = "male_prior_morbidity_risk_of_death_coefficients_new.tsv"

# female trajectory hospitalisation rate coefficients
fem_hosp_file = "female_trajectory_risk_of_hospitalisation_coefficients_new.tsv"

# female prior morbidity hospitalisation rate coefficients
fem_morb_hosp_file = "female_prior_morbidity_risk_of_hospitalisation_coefficients_new.tsv"

# male trajectory hospitalisation rate coefficients
male_hosp_file = "male_trajectory_risk_of_hospitalisation_coefficients_new.tsv"

# male prior morbidity hospitalisation rate coefficients
male_morb_hosp_file = "male_prior_morbidity_risk_of_hospitalisation_coefficients_new.tsv"

# female jaccard statistics

fem_stats_file = "ukb_fem_multimorbidity_ICD10_physio_block_1pct_stats.tsv"

# male jaccard statistics

male_stats_file = "ukb_male_multimorbidity_ICD10_physio_block_1pct_stats.tsv"

# ICD10 blocks dictionary
icd10_blocks_file = "icd10_blocks.tsv"

# functions

set_mt_col <- function(node, map) {
  if(is_empty(node$mt_p_val < 0.05 & node$std_mt > 0)) {
    return(paste0("/", non_sig_col, "/", fem_ns_cols))
  } else if(node$mt_p_val < 0.05) {
    return(paste0("/", map$scheme[which(map$std_mt == node$std_mt)], "/", map$col[which(map$std_mt == node$std_mt)]))
  } else {
    return(paste0("/", non_sig_col, "/", fem_ns_cols))
  }
}

set_node_col <- function(node) {
  return(block_col %>% filter(block == node$name) %>% pull(col) %>% unlist())
}

risk_class <- Vectorize(function(estimate, qnt_90, qnt_80, qnt_60, qnt_40, qnt_20) {
  if(estimate > qnt_90) {
    return(6)
  } else if(estimate > qnt_80) {
    return(5)
  } else if(estimate > qnt_60) {
    return(4)
  } else if(estimate > qnt_40) {
    return(3)
  } else if(estimate > qnt_20) {
    return(2)
  } else {
    return(1)
  }
})

# read in 

fem_sig <- fread(fem_sig_file)
male_sig <- fread(male_sig_file)

# read in mortality risk

fem_mort <- fread(fem_mort_file)
male_mort <- fread(male_mort_file)

fem_morb_mort <- fread(fem_morb_mort_file)
male_morb_mort <- fread(male_morb_mort_file)

# read in hospitalisation risk

fem_hosp <- fread(fem_hosp_file)
male_hosp <- fread(male_hosp_file)

fem_morb_hosp <- fread(fem_morb_hosp_file)
male_morb_hosp <- fread(male_morb_hosp_file)

# filter for common outputs

fem_cmn <- intersect(fem_mort$traj, fem_hosp$traj)

male_cmn <- intersect(male_mort$traj, male_hosp$traj)

# calculate confidence intervals

fem_mort %<>%
  as.tibble() %>%
  mutate(upp_ci = estimate + 1.96*std.error,
         low_ci = estimate - 1.96*std.error)

fem_hosp %<>%
  as.tibble() %>%
  mutate(upp_ci = estimate + 1.96*std.error,
         low_ci = estimate - 1.96*std.error)

male_mort %<>%
  as.tibble() %>%
  mutate(upp_ci = estimate + 1.96*std.error,
         low_ci = estimate - 1.96*std.error)

male_hosp %<>%
  as.tibble() %>%
  mutate(upp_ci = estimate + 1.96*std.error,
         low_ci = estimate - 1.96*std.error)


# convert to more interpretable effect

fem_mort %<>%
  mutate(estimate = 1/exp(estimate),
         upp_ci = 1/exp(upp_ci),
         low_ci = 1/exp(low_ci))
  
male_mort %<>%
  mutate(estimate = 1/exp(estimate),
         upp_ci = 1/exp(upp_ci),
         low_ci = 1/exp(low_ci))

fem_hosp %<>%
  mutate(estimate = 1/exp(estimate),
         upp_ci = 1/exp(upp_ci),
         low_ci = 1/exp(low_ci))

male_hosp %<>%
  mutate(estimate = 1/exp(estimate),
         upp_ci = 1/exp(upp_ci),
         low_ci = 1/exp(low_ci))

# read in ICD10 blocks

icd10_blocks <- fread(icd10_blocks_file)

# plot distributions

plot <- fem_morb_mort %>%
  filter(term == "morb") %>%
  mutate(estimate = 1/exp(estimate)) %>%
  ggplot(aes(x = estimate)) +
  stat_ecdf() +
  theme_bw() +
  xlab("Estimated increase in 1-year mortality\nrate per prior accrued diagnosis") +
  ylab("eCDF") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(1, 1.6, 0.1)) +
  geom_text(aes(label = paste0("Median = ", round(median(estimate), 2)), x = quantile(estimate, 0.85), y = 0.5), size = 8) +
  theme(text = element_text(size = 15))
ggsave(filename = "fem_morb_mort_risk_distribution.png", plot = plot,
       units = "px", width = 1600, height = 1400)

plot <- male_morb_mort %>%
  filter(term == "morb") %>%
  mutate(estimate = 1/exp(estimate)) %>%
  ggplot(aes(x = estimate)) +
  stat_ecdf() +
  theme_bw() +
  xlab("Estimated increase in 1-year mortality\nrate per prior accrued diagnosis") +
  ylab("eCDF") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(1, 1.6, 0.1)) +
  geom_text(aes(label = paste0("Median = ", round(median(estimate), 2)), x = quantile(estimate, 0.95), y = 0.5), size = 8) +
  theme(text = element_text(size = 15)) +
  coord_cartesian(xlim = c(1.1, 1.625))
ggsave(filename = "male_morb_mort_risk_distribution.png", plot = plot,
       units = "px", width = 1600, height = 1400)


plot <- fem_morb_hosp %>%
  filter(term == "morb") %>%
  mutate(estimate = 1/exp(estimate)) %>%
  ggplot(aes(x = estimate)) +
  stat_ecdf() +
  theme_bw() +
  xlab("Estimated increase in 1-year hospitalisation\nrate per prior accrued diagnosis") +
  ylab("eCDF") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  geom_text(aes(label = paste0("Median = ", round(median(estimate), 2)), x = quantile(estimate, 0.05), y = 0.5), size = 8) +
  theme(text = element_text(size = 15))
ggsave(filename = "fem_morb_hosp_risk_distribution.png", plot = plot,
       units = "px", width = 1600, height = 1400)

plot <- male_morb_hosp %>%
  filter(term == "morb") %>%
  mutate(estimate = 1/exp(estimate)) %>%
  ggplot(aes(x = estimate)) +
  stat_ecdf() +
  theme_bw() +
  xlab("Estimated increase in 1-year hospitalisation\nrate per prior accrued diagnosis") +
  ylab("eCDF") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(1, 1.25, 0.05)) +
  geom_text(aes(label = paste0("Median = ", round(median(estimate), 2)), x = quantile(estimate, 0.05), y = 0.5), size = 8) +
  theme(text = element_text(size = 15))
ggsave(filename = "male_morb_hosp_risk_distribution.png", plot = plot,
       units = "px", width = 1600, height = 1400)

# effect of general multimorbidity burden on mortality

fem_tmp <- fem_morb_mort %>%
  as_tibble() %>%
  mutate(upp_ci = estimate + 1.96*std.error,
         low_ci = estimate - 1.96*std.error) %>%
  filter(term == "morb") %>%
  mutate(estimate = 1/exp(estimate),
         upp_ci = 1/exp(upp_ci),
         low_ci = 1/exp(low_ci)) %>%
  mutate(sex = "female")

male_tmp <- male_morb_mort %>%
  as_tibble() %>%
  mutate(upp_ci = estimate + 1.96*std.error,
         low_ci = estimate - 1.96*std.error) %>%
  filter(term == "morb") %>%
  mutate(estimate = 1/exp(estimate),
         upp_ci = 1/exp(upp_ci),
         low_ci = 1/exp(low_ci)) %>%
  mutate(sex = "male")

# females

fem_tmp %<>%
  select(disease, estimate, low_ci, upp_ci, sex)

plot <- fem_tmp %>%
  mutate(disease = factor(disease, levels = disease[order(estimate)])) %>%
  ggplot(aes(y = disease, x = estimate, xmin = low_ci, xmax = upp_ci)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_errorbarh() +
  xlab("Estimated increase in 1-year mortality\nrate per prior accrued diagnosis") +
  ylab("Presenting diagnosis") +
  coord_cartesian(xlim = c(1, 1.7)) +
  geom_vline(xintercept = 1, lty = 2) +
  scale_x_continuous(breaks = seq(1, 1.9, 0.1))
ggsave(filename = "fem_morb_mort_risk_presenting_disease.png", plot = plot,
       units = "px", width = 1200, height = 2400)

# males

male_tmp %<>%
  select(disease, estimate, low_ci, upp_ci, sex)

plot <- male_tmp %>%
  mutate(disease = factor(disease, levels = disease[order(estimate)])) %>%
  ggplot(aes(y = disease, x = estimate, xmin = low_ci, xmax = upp_ci)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_errorbarh() +
  xlab("Estimated increase in 1-year mortality\nrate per prior accrued diagnosis") +
  ylab("Presenting diagnosis") +
  coord_cartesian(xlim = c(1, 1.65)) +
  geom_vline(xintercept = 1, lty = 2) +
  scale_x_continuous(breaks = seq(1, 1.6, 0.1))
ggsave(filename = "male_morb_mort_risk_presenting_disease.png", plot = plot,
       units = "px", width = 1200, height = 2400)

# female vs male largest differences

fem_tmp2 <- fem_tmp %>%
  select(-sex, -estimate)
colnames(fem_tmp2) <- c("disease", "fem_low", "fem_upp")

male_tmp2 <- male_tmp %>%
  select(-sex, -estimate)
colnames(male_tmp2) <- c("disease", "male_low", "male_upp")

fm_tmp <- fem_tmp2 %>%
  full_join(male_tmp2, "disease")

non_ol <- fm_tmp$disease[which(fm_tmp$fem_low < fm_tmp$male_upp | fm_tmp$male_low < fm_tmp$fem_upp)]

tmp_diff <- rbind(fem_tmp,
                  male_tmp) %>%
  pivot_wider(disease, names_from = "sex", values_from = "estimate") %>%
  mutate(diff = female - male) %>%
  filter(!is.na(diff)) %>%
  mutate(disease = ifelse(disease %in% non_ol, paste0("*** ", disease), disease)) %>%
  mutate(disease = factor(disease, levels = disease[order(diff)]))

plot <- rbind(fem_tmp,
      male_tmp) %>%
  mutate(disease = ifelse(disease %in% non_ol, paste0("*** ", disease), disease)) %>%
  mutate(disease = factor(disease, levels = levels(tmp_diff$disease))) %>%
  filter(!is.na(disease)) %>%
  ggplot(aes(y = disease, x = estimate, xmin = low_ci, xmax = upp_ci, col = factor(sex))) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_errorbarh() +
  xlab("Estimated increase in 1-year mortality\nrate per prior accrued diagnosis") +
  ylab("Presenting diagnosis (Ordered by difference in effect between females and males)") +
  coord_cartesian(xlim = c(1, 1.65)) +
  geom_vline(xintercept = 1, lty = 2) +
  scale_x_continuous(breaks = seq(1, 1.6, 0.1)) +
  labs(col = "Sex") +
  scale_colour_discrete(labels = c(female = "Female", male = "Male"))
ggsave(filename = "fem_vs_male_morb_mort_risk_presenting_disease.png", plot = plot,
       units = "px", width = 1200, height = 2400)

# effect of general multimorbidity on hospitalisation

fem_tmp <- fem_morb_hosp %>%
  as_tibble() %>%
  mutate(upp_ci = estimate + 1.96*std.error,
         low_ci = estimate - 1.96*std.error) %>%
  filter(term == "morb") %>%
  mutate(estimate = 1/exp(estimate),
         upp_ci = 1/exp(upp_ci),
         low_ci = 1/exp(low_ci)) %>%
  mutate(sex = "female")

male_tmp <- male_morb_hosp %>%
  as_tibble() %>%
  mutate(upp_ci = estimate + 1.96*std.error,
         low_ci = estimate - 1.96*std.error) %>%
  filter(term == "morb") %>%
  mutate(estimate = 1/exp(estimate),
         upp_ci = 1/exp(upp_ci),
         low_ci = 1/exp(low_ci)) %>%
  mutate(sex = "male")

# females

fem_tmp %<>%
  select(disease, estimate, low_ci, upp_ci, sex)

plot <- fem_tmp %>%
  mutate(disease = factor(disease, levels = disease[order(estimate)])) %>%
  ggplot(aes(y = disease, x = estimate, xmin = low_ci, xmax = upp_ci)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_errorbarh() +
  xlab("Estimated increase in 1-year hospitalisation\nrate per prior accrued diagnosis") +
  ylab("Presenting diagnosis") +
  coord_cartesian(xlim = c(1, 1.25)) +
  geom_vline(xintercept = 1, lty = 2) +
  scale_x_continuous(breaks = seq(1, 1.6, 0.05))
ggsave(filename = "fem_morb_hosp_risk_presenting_disease.png", plot = plot,
       units = "px", width = 1200, height = 2400)

# males

male_tmp %<>%
  select(disease, estimate, low_ci, upp_ci, sex)

plot <- male_tmp %>%
  mutate(disease = factor(disease, levels = disease[order(estimate)])) %>%
  ggplot(aes(y = disease, x = estimate, xmin = low_ci, xmax = upp_ci)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_errorbarh() +
  xlab("Estimated increase in 1-year hospitalisation\nrate per prior accrued diagnosis") +
  ylab("Presenting diagnosis") +
  coord_cartesian(xlim = c(1, 1.275)) +
  geom_vline(xintercept = 1, lty = 2) +
  scale_x_continuous(breaks = seq(1, 1.6, 0.05))
ggsave(filename = "male_morb_hosp_risk_presenting_disease.png", plot = plot,
       units = "px", width = 1200, height = 2400)

# female vs male largest differences

fem_tmp2 <- fem_tmp %>%
  select(-sex, -estimate)
colnames(fem_tmp2) <- c("disease", "fem_low", "fem_upp")

male_tmp2 <- male_tmp %>%
  select(-sex, -estimate)
colnames(male_tmp2) <- c("disease", "male_low", "male_upp")

fm_tmp <- fem_tmp2 %>%
  full_join(male_tmp2, "disease")

non_ol <- fm_tmp$disease[which(fm_tmp$fem_low < fm_tmp$male_upp | fm_tmp$male_low < fm_tmp$fem_upp)]

tmp_diff <- rbind(fem_tmp,
                  male_tmp) %>%
  pivot_wider(disease, names_from = "sex", values_from = "estimate") %>%
  mutate(diff = female - male) %>%
  filter(!is.na(diff)) %>%
  mutate(disease = ifelse(disease %in% non_ol, paste0("*** ", disease), disease)) %>%
  mutate(disease = factor(disease, levels = disease[order(diff)]))

plot <- rbind(fem_tmp,
      male_tmp) %>%
  mutate(disease = ifelse(disease %in% non_ol, paste0("*** ", disease), disease)) %>%
  mutate(disease = factor(disease, levels = levels(tmp_diff$disease))) %>%
  filter(!is.na(disease)) %>%
  ggplot(aes(y = disease, x = estimate, xmin = low_ci, xmax = upp_ci, col = factor(sex))) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_errorbarh() +
  xlab("Estimated increase in 1-year hospitalisation\nrate per prior accrued diagnosis") +
  ylab("Presenting diagnosis (Ordered by difference in effect between females and males)") +
  coord_cartesian(xlim = c(1, 1.275)) +
  geom_vline(xintercept = 1, lty = 2) +
  scale_x_continuous(breaks = seq(1, 1.6, 0.05)) +
  labs(col = "Sex") +
  scale_colour_discrete(labels = c(female = "Female", male = "Male"))
ggsave(filename = "fem_vs_male_morb_hosp_risk_presenting_disease.png", plot = plot,
       units = "px", width = 1200, height = 2400)

# where are single disease histories effective on outcomes?

plot <- fem_mort %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(hist_l == 1) %>%
  filter(str_detect(term, "hist")) %>%
  mutate(estimate = ifelse(adj_p_val < 0.05, estimate, 1)) %>%
  separate(traj, into = c("hist1", "pres"), sep = " -> ") %>%
  mutate(chapter = icd10_blocks$chapter[match(hist1, icd10_blocks$block)]) %>%
  ggplot(aes(x = pres, y = log10(estimate), col = factor(ifelse(chapter == "C00-D49", "Cancer", "Non-cancer")))) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(col = "Diagnosis\nhistory type\n") +
  xlab("Presenting diagnosis") +
  ylab("Log10(estimated increase\nin 1-year mortality rate with\nsingle diagnosis history)") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_colour_viridis_d(option = "H", begin = 0.35, end = 1)
ggsave(filename = "fem_traj_mort_risk_presenting_disease.png", plot = plot,
       units = "px", width = 2800, height = 800)

fem_tmp <- fem_mort %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(hist_l == 1) %>%
  filter(str_detect(term, "hist")) %>%
  mutate(estimate = ifelse(adj_p_val < 0.05, estimate, 1)) %>%
  separate(traj, into = c("hist1", "pres"), sep = " -> ") %>%
  mutate(chapter = icd10_blocks$chapter[match(hist1, icd10_blocks$block)]) %>%
  filter(chapter != "C00-D49") %>%
  filter(estimate > 1) %>%
  select(pres, hist1, estimate, upp_ci, low_ci, upp_ci, adj_p_val, pres_size, hist_size, traj_size) %>%
  arrange(pres, hist1)
write_tsv(fem_tmp, "fem_traj_mort_risk_non_cancer_single_disease_histories.tsv")
  

plot <- male_mort %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(hist_l == 1) %>%
  filter(str_detect(term, "hist")) %>%
  mutate(estimate = ifelse(adj_p_val < 0.05, estimate, 1)) %>%
  separate(traj, into = c("hist1", "pres"), sep = " -> ") %>%
  mutate(chapter = icd10_blocks$chapter[match(hist1, icd10_blocks$block)]) %>%
  ggplot(aes(x = pres, y = log10(estimate), col = factor(ifelse(chapter == "C00-D49", "Cancer", "Non-cancer")))) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(col = "Diagnosis\nhistory type\n") +
  xlab("Presenting diagnosis") +
  ylab("Log10(estimated increase\nin 1-year mortality rate with\nsingle diagnosis history)") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_colour_viridis_d(option = "H", begin = 0.35, end = 1)
ggsave(filename = "male_traj_mort_risk_presenting_disease.png", plot = plot,
       units = "px", width = 2800, height = 800)

male_tmp <- male_mort %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(hist_l == 1) %>%
  filter(str_detect(term, "hist")) %>%
  mutate(estimate = ifelse(adj_p_val < 0.05, estimate, 1)) %>%
  separate(traj, into = c("hist1", "pres"), sep = " -> ") %>%
  mutate(chapter = icd10_blocks$chapter[match(hist1, icd10_blocks$block)]) %>%
  filter(chapter != "C00-D49") %>%
  filter(estimate > 1) %>%
  select(pres, hist1, estimate, upp_ci, low_ci, upp_ci, adj_p_val, pres_size, hist_size, traj_size) %>%
  arrange(pres, hist1)
write_tsv(male_tmp, "male_traj_mort_risk_non_cancer_single_disease_histories.tsv")


plot <- fem_hosp %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(hist_l == 1) %>%
  filter(str_detect(term, "hist")) %>%
  mutate(estimate = ifelse(adj_p_val < 0.05, estimate, 1)) %>%
  separate(traj, into = c("hist1", "pres"), sep = " -> ") %>%
  mutate(chapter = icd10_blocks$chapter[match(hist1, icd10_blocks$block)]) %>%
  ggplot(aes(x = pres, y = log10(estimate), col = factor(ifelse(chapter == "C00-D49", "Cancer", "Non-cancer")))) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(col = "Diagnosis\nhistory type\n") +
  xlab("Presenting diagnosis") +
  ylab("Log10(estimated increase\nin 1-year hospitalisation rate\nwith single diagnosis history)") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_colour_viridis_d(option = "H", begin = 0.35, end = 1)
ggsave(filename = "fem_traj_hosp_risk_presenting_disease.png", plot = plot,
       units = "px", width = 2800, height = 800)

fem_tmp <- fem_hosp %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(hist_l == 1) %>%
  filter(str_detect(term, "hist")) %>%
  mutate(estimate = ifelse(adj_p_val < 0.05, estimate, 1)) %>%
  separate(traj, into = c("hist1", "pres"), sep = " -> ") %>%
  mutate(chapter = icd10_blocks$chapter[match(hist1, icd10_blocks$block)]) %>%
  filter(chapter != "C00-D49") %>%
  filter(estimate > 1) %>%
  select(pres, hist1, estimate, upp_ci, low_ci, upp_ci, adj_p_val, pres_size, hist_size, traj_size) %>% 
  arrange(pres, hist1)
write_tsv(fem_tmp, "fem_traj_hosp_risk_non_cancer_single_disease_histories.tsv")


plot <- male_hosp %>% 
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(hist_l == 1) %>%
  filter(str_detect(term, "hist")) %>%
  mutate(estimate = ifelse(adj_p_val < 0.05, estimate, 1)) %>%
  separate(traj, into = c("hist1", "pres"), sep = " -> ") %>%
  mutate(chapter = icd10_blocks$chapter[match(hist1, icd10_blocks$block)]) %>%
  ggplot(aes(x = pres, y = log10(estimate), col = factor(ifelse(chapter == "C00-D49", "Cancer", "Non-cancer")))) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(col = "Diagnosis\nhistory type\n") +
  xlab("Presenting diagnosis") +
  ylab("Log10(estimated increase\nin 1-year hospitalisation rate\nwith single diagnosis history)") +
  geom_hline(yintercept = 0, lty = 2) +
  scale_colour_viridis_d(option = "H", begin = 0.35, end = 1)
ggsave(filename = "male_traj_hosp_risk_presenting_disease.png", plot = plot,
       units = "px", width = 2800, height = 800)

male_tmp <- male_hosp %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(hist_l == 1) %>%
  filter(str_detect(term, "hist")) %>%
  mutate(estimate = ifelse(adj_p_val < 0.05, estimate, 1)) %>%
  separate(traj, into = c("hist1", "pres"), sep = " -> ") %>%
  mutate(chapter = icd10_blocks$chapter[match(hist1, icd10_blocks$block)]) %>%
  filter(chapter != "C00-D49") %>%
  filter(estimate > 1) %>%
  select(pres, hist1, estimate, upp_ci, low_ci, upp_ci, adj_p_val, pres_size, hist_size, traj_size) %>% 
  arrange(pres, hist1)
write_tsv(male_tmp, "male_traj_hosp_risk_non_cancer_single_disease_histories.tsv")

# are there significant interaction effects?

# mortality

fem_mort_inter <- fem_mort %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(str_detect(term, ":")) %>%
  filter(adj_p_val < 0.05) %>%
  pull(traj)

# calculate resultant effect

fem_mort_inter_est <- fem_mort %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(traj %in% fem_mort_inter) %>%
  filter(str_detect(term, "hist")) %>%
  group_by(traj) %>%
  mutate(type = ifelse(estimate[term == "hist11:hist21"] > 1, "Positively interactive", "Negatively interactive")) %>%
  summarise(estimate = reduce(estimate, `*`),
            type = type,
            pres_size = pres_size,
            hist_size = hist_size,
            traj_size = traj_size) %>%
  unique.data.frame()

fem_mort_inter_p <- fem_mort %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(traj %in% fem_mort_inter) %>%
  filter(str_detect(term, "hist")) %>%
  pivot_wider(traj, names_from = term, values_from = adj_p_val)

fem_mort_inter_est %<>%
  full_join(fem_mort_inter_p, "traj")

# males

male_mort_inter <- male_mort %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(str_detect(term, ":")) %>%
  filter(adj_p_val < 0.05) %>%
  pull(traj)

male_mort_inter_est <- male_mort %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(traj %in% male_mort_inter) %>%
  filter(str_detect(term, "hist")) %>%
  group_by(traj) %>%
  mutate(type = ifelse(estimate[term == "hist11:hist21"] > 1, "Positively interactive", "Negatively interactive")) %>%
  summarise(estimate = reduce(estimate, `*`),
            type = type,
            pres_size = pres_size,
            hist_size = hist_size,
            traj_size = traj_size) %>%
  unique.data.frame()

male_mort_inter_p <- male_mort %>%
  mutate(hist_l = str_count(traj, "->")) %>%
  filter(traj %in% male_mort_inter) %>%
  filter(str_detect(term, "hist")) %>%
  pivot_wider(traj, names_from = term, values_from = adj_p_val)

male_mort_inter_est %<>%
  full_join(male_mort_inter_p, "traj")

# hospitalisation

# are there significant additive effects of histories?

# mortality

fem_mort %<>%
  mutate(hist_len = str_count(traj, "->"))

fem_mort_add_est <- fem_mort %>%
  filter(hist_len > 1) %>%
  filter(!str_detect(term, ":")) %>%
  filter(str_detect(term, "hist")) %>%
  filter(!traj %in% fem_mort_inter) %>%
  group_by(traj) %>%
  filter(all(adj_p_val < 0.05)) %>%
  mutate(type =  ifelse(all(estimate > 1) | all(estimate < 1), "Concordantly additive", "Discordantly additive")) %>%
  summarise(estimate = reduce(estimate, `*`),
            type = type,
            pres_size = pres_size,
            hist_size = hist_size,
            traj_size = traj_size) %>%
  unique.data.frame()

fem_mort_add_p <- fem_mort %>%
  filter(hist_len == 2) %>%
  # filter(!str_detect(term, ":")) %>%
  filter(str_detect(term, "hist")) %>%
  filter(!traj %in% fem_mort_inter) %>%
  pivot_wider(traj, names_from = term, values_from = adj_p_val) %>%
  filter(hist11 < 0.05 & hist21 < 0.05)

fem_mort_add_est %<>%
  full_join(fem_mort_add_p, "traj")
  
fem_mort_multi <- rbind(fem_mort_add_est,
                        fem_mort_inter_est) %>%
  # mutate(type = ifelse(traj %in% fem_mort_inter, "Interactive", "Additive")) %>%
  separate(traj, into = c("hist2", "hist1", "pres"), sep = " -> ") %>%
  select(pres, hist1, hist2, estimate, type, hist11, hist21, `hist11:hist21`, pres_size, hist_size, traj_size)

plot <- fem_mort_multi %>%
  arrange(desc(estimate)) %>%
  mutate(chapt1 = icd10_blocks$chapter[match(hist1, icd10_blocks$block)],
         chapt2 = icd10_blocks$chapter[match(hist2, icd10_blocks$block)]) %>%
  mutate(chapt = ifelse(chapt1 == "C00-D49" & chapt2 == "C00-D49", "Cancer + cancer", "Non-cancer + non-cancer")) %>%
  mutate(chapt = ifelse(chapt1 == "C00-D49" & chapt2 != "C00-D49", "Cancer + non-cancer", chapt)) %>%
  mutate(chapt = ifelse(chapt1 != "C00-D49" & chapt2 == "C00-D49", "Cancer + non-cancer", chapt)) %>%
  ggplot(aes(x = pres, y = log10(estimate), col = factor(chapt), shape = type)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Log10(estimated increase\nin 1-year mortality rate with\nmulti-diagnosis history)") +
  xlab("Presenting diagnosis") +
  scale_colour_viridis_d(option = "H", begin = 0.35, end = 1) +
  theme_bw() +
  labs(shape = "Effect type of\nmulti-diagnosis\nhistory",
       col = "Multi-diagnosis\nhistory type") +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 14))
ggsave(filename = "fem_traj_mort_multi_disease_history_presenting_disease.png", plot = plot,
              units = "px", width = 1600, height = 1200)


# write_tsv(fem_tmp, "fem_traj_mort_risk_disease_histories_additive_effects.tsv")

male_mort %<>%
  mutate(hist_len = str_count(traj, "->"))


male_mort_add_est <- male_mort %>%
  filter(hist_len > 1) %>%
  filter(!str_detect(term, ":")) %>%
  filter(str_detect(term, "hist")) %>%
  filter(!traj %in% male_mort_inter) %>%
  group_by(traj) %>%
  filter(all(adj_p_val < 0.05)) %>%
  mutate(type =  ifelse(all(estimate > 1) | all(estimate < 1), "Concordantly additive", "Discordantly additive")) %>%
  summarise(estimate = reduce(estimate, `*`),
            type = type,
            pres_size = pres_size,
            hist_size = hist_size,
            traj_size = traj_size) %>%
  unique.data.frame()

male_mort_add_p <- male_mort %>%
  filter(hist_len == 2) %>%
  # filter(!str_detect(term, ":")) %>%
  filter(str_detect(term, "hist")) %>%
  filter(!traj %in% male_mort_inter) %>%
  pivot_wider(traj, names_from = term, values_from = adj_p_val) %>%
  filter(hist11 < 0.05 & hist21 < 0.05)

male_mort_add_est %<>%
  full_join(male_mort_add_p, "traj")

male_mort_multi <- rbind(male_mort_add_est,
                        male_mort_inter_est) %>%
  # mutate(type = ifelse(traj %in% male_mort_inter, "Interactive", "Additive")) %>%
  separate(traj, into = c("hist2", "hist1", "pres"), sep = " -> ") %>%
  select(pres, hist1, hist2, estimate, type, hist11, hist21, `hist11:hist21`, pres_size, hist_size, traj_size)

plot <- male_mort_multi %>%
  arrange(desc(estimate)) %>%
  mutate(chapt1 = icd10_blocks$chapter[match(hist1, icd10_blocks$block)],
         chapt2 = icd10_blocks$chapter[match(hist2, icd10_blocks$block)]) %>%
  mutate(chapt = ifelse(chapt1 == "C00-D49" & chapt2 == "C00-D49", "Cancer + cancer", "Non-cancer + non-cancer")) %>%
  mutate(chapt = ifelse(chapt1 == "C00-D49" & chapt2 != "C00-D49", "Cancer + non-cancer", chapt)) %>%
  mutate(chapt = ifelse(chapt1 != "C00-D49" & chapt2 == "C00-D49", "Cancer + non-cancer", chapt)) %>%
  ggplot(aes(x = pres, y = log10(estimate), col = factor(chapt), shape = type)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Log10(estimated increase\nin 1-year mortality rate with\nmulti-diagnosis history)") +
  xlab("Presenting diagnosis") +
  scale_colour_viridis_d(option = "H", begin = 0.35, end = 1) +
  theme_bw() +
  labs(shape = "Effect type of\nmulti-diagnosis\nhistory",
       col = "Multi-diagnosis\nhistory type") +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 14))
ggsave(filename = "male_traj_mort_multi_disease_history_presenting_disease.png", plot = plot,
       units = "px", width = 1800, height = 1200)

# hospitalisation

fem_hosp %<>%
  mutate(hist_len = str_count(traj, "->"))

fem_hosp_add_est <- fem_hosp %>%
  filter(hist_len > 1) %>%
  filter(!str_detect(term, ":")) %>%
  filter(str_detect(term, "hist")) %>%
  # filter(!traj %in% fem_hosp_inter) %>%
  group_by(traj) %>%
  filter(all(adj_p_val < 0.05)) %>%
  mutate(type =  ifelse(all(estimate > 1) | all(estimate < 1), "Concordantly additive", "Discordantly additive")) %>%
  summarise(estimate = reduce(estimate, `*`),
            type = type,
            pres_size = pres_size,
            hist_size = hist_size,
            traj_size = traj_size) %>%
  unique.data.frame()

fem_hosp_add_p <- fem_hosp %>%
  filter(hist_len == 2) %>%
  # filter(!str_detect(term, ":")) %>%
  filter(str_detect(term, "hist")) %>%
  pivot_wider(traj, names_from = term, values_from = adj_p_val) %>%
  filter(hist11 < 0.05 & hist21 < 0.05)

fem_hosp_add_est %<>%
  full_join(fem_hosp_add_p, "traj")

fem_hosp_multi <- fem_hosp_add_est %>%
  # mutate(type = "Additive") %>%
  separate(traj, into = c("hist2", "hist1", "pres"), sep = " -> ") %>%
  select(pres, hist1, hist2, estimate, type, hist11, hist21, `hist11:hist21`, pres_size, hist_size, traj_size)

plot <- fem_hosp_multi %>%
  arrange(desc(estimate)) %>%
  mutate(chapt1 = icd10_blocks$chapter[match(hist1, icd10_blocks$block)],
         chapt2 = icd10_blocks$chapter[match(hist2, icd10_blocks$block)]) %>%
  mutate(chapt = ifelse(chapt1 == "C00-D49" & chapt2 == "C00-D49", "Cancer + cancer", "Non-cancer + non-cancer")) %>%
  mutate(chapt = ifelse(chapt1 == "C00-D49" & chapt2 != "C00-D49", "Cancer + non-cancer", chapt)) %>%
  mutate(chapt = ifelse(chapt1 != "C00-D49" & chapt2 == "C00-D49", "Cancer + non-cancer", chapt)) %>%
  ggplot(aes(x = pres, y = log10(estimate), col = factor(chapt ,levels = c("Cancer + cancer", "Cancer + non-cancer", "Non-cancer + non-cancer")), shape = factor(type, levels = c("Concordantly additive", "Discordantly additive", "Negatively interactive")))) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Log10(estimated increase\nin 1-year hospitalisation rate with\nmulti-diagnosis history)") +
  xlab("Presenting diagnosis") +
  scale_shape(drop = F) +
  scale_colour_viridis_d(option = "H", begin = 0.35, end = 1, drop = FALSE) +
  theme_bw() +
  labs(shape = "Effect type of\nmulti-diagnosis\nhistory",
       col = "Multi-diagnosis\nhistory type") +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 14))
ggsave(filename = "fem_traj_hosp_multi_disease_history_presenting_disease.png", plot = plot,
       units = "px", width = 1600, height = 1200)


# write_tsv(fem_tmp, "fem_traj_hosp_risk_disease_histories_additive_effects.tsv")

male_hosp %<>%
  mutate(hist_len = str_count(traj, "->"))

male_hosp_add_est <- male_hosp %>%
  filter(hist_len > 1) %>%
  filter(!str_detect(term, ":")) %>%
  filter(str_detect(term, "hist")) %>%
  # filter(!traj %in% male_hosp_inter) %>%
  group_by(traj) %>%
  filter(all(adj_p_val < 0.05)) %>%
  mutate(type =  ifelse(all(estimate > 1) | all(estimate < 1), "Concordantly additive", "Discordantly additive")) %>%
  summarise(estimate = reduce(estimate, `*`),
            type = type,
            pres_size = pres_size,
            hist_size = hist_size,
            traj_size = traj_size) %>%
  unique.data.frame()

male_hosp_add_p <- male_hosp %>%
  filter(hist_len == 2) %>%
  # filter(!str_detect(term, ":")) %>%
  filter(str_detect(term, "hist")) %>%
  pivot_wider(traj, names_from = term, values_from = adj_p_val) %>%
  filter(hist11 < 0.05 & hist21 < 0.05)

male_hosp_add_est %<>%
  full_join(male_hosp_add_p, "traj")

male_hosp_multi <- male_hosp_add_est %>%
  # mutate(type = "Additive") %>%
  separate(traj, into = c("hist2", "hist1", "pres"), sep = " -> ") %>%
  select(pres, hist1, hist2, estimate, type, hist11, hist21, `hist11:hist21`, pres_size, hist_size, traj_size)

plot <- male_hosp_multi %>%
  arrange(desc(estimate)) %>%
  mutate(chapt1 = icd10_blocks$chapter[match(hist1, icd10_blocks$block)],
         chapt2 = icd10_blocks$chapter[match(hist2, icd10_blocks$block)]) %>%
  mutate(chapt = ifelse(chapt1 == "C00-D49" & chapt2 == "C00-D49", "Cancer + cancer", "Non-cancer + non-cancer")) %>%
  mutate(chapt = ifelse(chapt1 == "C00-D49" & chapt2 != "C00-D49", "Cancer + non-cancer", chapt)) %>%
  mutate(chapt = ifelse(chapt1 != "C00-D49" & chapt2 == "C00-D49", "Cancer + non-cancer", chapt)) %>%
  ggplot(aes(x = pres, y = log10(estimate), col = factor(chapt ,levels = c("Cancer + cancer", "Cancer + non-cancer", "Non-cancer + non-cancer")), shape = factor(type, levels = c("Concordantly additive", "Discordantly additive", "Negatively interactive")))) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Log10(estimated increase\nin 1-year hospitalisation rate with\nmulti-diagnosis history)") +
  xlab("Presenting diagnosis") +
  scale_shape(drop = F) +
  scale_colour_viridis_d(option = "H", begin = 0.35, end = 1, drop = FALSE) +
  theme_bw() +
  labs(shape = "Effect type of\nmulti-diagnosis\nhistory",
       col = "Multi-diagnosis\nhistory type") +
  theme(axis.text.x = element_text(angle = 90),
        text = element_text(size = 14))
ggsave(filename = "male_traj_hosp_multi_disease_history_presenting_disease.png", plot = plot,
       units = "px", width = 1800, height = 1200)

fem_mort_multi %<>%
  arrange(pres)

fem_hosp_multi %<>%
  arrange(pres)

male_mort_multi %<>%
  arrange(pres)

male_hosp_multi %<>%
  arrange(pres)

write_tsv(fem_mort_multi, "fem_traj_mort_risk_multi_disease_histories.tsv")
write_tsv(fem_hosp_multi, "fem_traj_hosp_risk_multi_disease_histories.tsv")

write_tsv(male_mort_multi, "male_traj_mort_risk_multi_disease_histories.tsv")
write_tsv(male_hosp_multi, "male_traj_hosp_risk_multi_disease_histories.tsv")