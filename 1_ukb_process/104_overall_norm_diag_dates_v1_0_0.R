# normalise diagnosis dates

rm(list = ls())

library(tidyverse)

# Variables

# input

# female diagnosis dates matrix
fem_traj_dates_file = "fem_diag_block_1pct_dates.tsv"

# male diagnosis dates matrix
male_traj_dates_file = "male_diag_block_1pct_dates.tsv"

# output

# female diagnosis time from first diagnosis dates
fem_norm_dates_file = "fem_diag_block_1pct_norm_dates.tsv"
# male diagnosis time from first diagnosis dates
male_norm_dates_file = "male_diag_block_1pct_norm_dates.tsv"

# female quartile thresholds of distribution of time from first diagnosis
fem_interval_file = "fem_1pct_block_quartile_interval_norm_dates.tsv"
# male quartile thresholds of distribution of time from first diagnosis
male_interval_file = "male_1pct_block_quartile_interval_norm_dates.tsv"

# read in dates

fem_dates <- read.table(fem_traj_dates_file, sep = "\t")
fem_dates <- as.matrix(fem_dates)

male_dates <- read.table(male_traj_dates_file, sep = "\t")
male_dates <- as.matrix(male_dates)

# convert to relative time differences

fem_dates_norm <- fem_dates

for(i in 1:nrow(fem_dates)) {
  dates_i <- as.Date(fem_dates[i,])
  dates_i_n <- dates_i - dates_i[1]
  fem_dates_norm[i,] <- dates_i_n

}

male_dates_norm <- male_dates

for(i in 1:nrow(male_dates)) {
  dates_i <- as.Date(male_dates[i,])
  dates_i_n <- dates_i - dates_i[1]
  male_dates_norm[i,] <- dates_i_n
}

fem_dates_norm <- apply(fem_dates_norm, 2, as.numeric)
male_dates_norm <- apply(male_dates_norm, 2, as.numeric)

fem_days <- as.vector(fem_dates_norm)

male_days <- as.vector(male_dates_norm)

# remove NAs

fem_na_ind <- which(is.na(fem_days))
fem_days <- fem_days[-fem_na_ind]

fem_days_df <- as.data.frame(fem_days)
fem_diag_qnts <- sapply(c(seq(0.25, 0.75, 0.25), 0.9), function(x) quantile(fem_days, x))

male_na_ind <- which(is.na(male_days))
male_days <- male_days[-male_na_ind]

male_days_df <- as.data.frame(male_days)
male_diag_qnts <- sapply(c(seq(0.25, 0.75, 0.25), 0.9), function(x) quantile(male_days, x))

days_df <- data.frame(days = c(fem_days, male_days)) %>%
  mutate(sex = c(rep("fem", nrow(fem_days_df)), rep("male", nrow(male_days_df))))

# write out normalised diagnosis intervals

write_tsv(as.data.frame(fem_diag_qnts), fem_interval_file)

write_tsv(as.data.frame(male_diag_qnts), male_interval_file)

# plot diagnosis ecdf

fem_col <- RColorBrewer::brewer.pal(11, name = "RdYlBu")[2]
male_col <- RColorBrewer::brewer.pal(11, name = "RdYlBu")[10]

plot <- days_df %>%
  mutate(days = days/365) %>%
  ggplot(aes(x = days, col = sex)) +
  stat_ecdf() +
  scale_colour_manual(values = c(fem_col, male_col), labels = c("Female", "Male")) +
  geom_vline(xintercept = fem_diag_qnts/365, col = fem_col, lty = 2) +
  geom_vline(xintercept = male_diag_qnts/365, col = male_col, lty = 2) +
  theme_bw() +
  labs(col = "Sex") +
  xlab("Time from first diagnosis (years)") +
  ylab("Cumulative proportion") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(0, 22, 4)) +
  theme(text = element_text(size = 25))
ggsave("fm_days_from_first_diag_distrib_quartiles.png",  plot = plot,
       units = "px",
       width = 2000, height = 1800)

# write out normalised dates

write_tsv(as.data.frame(fem_dates_norm), fem_norm_dates_file)

write_tsv(as.data.frame(male_dates_norm), male_norm_dates_file)
