# get morbidity distributions of each sex

library(tidyverse)

# female prevalence filtered multimorbidity matrix
fem_mm_file = "ukb_fem_multimorbidity_ICD10_physio_block_1pct_table.tsv"

# male prevalence filtered multimorbidity matrix
male_mm_file = "ukb_male_multimorbidity_ICD10_physio_block_1pct_table.tsv"

# read in mm

fem_mm <- read_tsv(fem_mm_file)
male_mm <- read_tsv(male_mm_file)

# get morbidity

fem_morb <- rowSums(fem_mm)
male_morb <- rowSums(male_mm)

# plot

fem_col <- RColorBrewer::brewer.pal(11, name = "RdYlBu")[2]
male_col <- RColorBrewer::brewer.pal(11, name = "RdYlBu")[10]

plot <- data.frame(morb = c(fem_morb, male_morb), 
           sex = c(rep("fem", length(fem_morb)), rep("male", length(male_morb)))) %>%
  ggplot(aes(x = morb, col = sex)) +
  stat_ecdf() +
  theme_bw() +
  scale_colour_manual(values = c(fem_col, male_col), labels = c("Female", "Male")) +
  labs(col = "Sex") +
  xlab("Length of patient diagnoses record\n(Number of ICD10 blocks)") +
  ylab("Cumulative proportion") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(xlim = c(0, 17)) +
  theme(text = element_text(size = 25))
ggsave("fm_block_morb_distrib.png", plot = plot,
       units = "px", width = 2000, height = 1800)




