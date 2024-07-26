# distribution value scatterplots

rm(list = ls())

print("Loading libraries...")

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(magrittr)
library(ggplot2)
library(sf)
library(sp)  
library(ggrepel)

# Variables

# ICD10 blocks dictionary
icd10_blocks_file = "icd10_blocks.tsv"

# functions

get_distrib <- function(type, data_type, gender) {
  # read in disease trajectoires
  
  filename <- paste0(gender, "_", data_type, "_", type)
  
  trajs <- read_tsv(traj_file, col_types = cols(.default = "c"))
  
  colnames(trajs) <- paste0("diag_", 1:ncol(trajs))
  
  # read in diagnosis dates
  
  dates <- read.table(traj_dates_file, sep = "\t", header = T)
  
  # read in blocks file
  
  blocks <- read_tsv(icd10_blocks_file)
  
  # remove all NAs
  
  all_na_inds <- which(apply(trajs, 1, function(x) sum(is.na(x))) == ncol(trajs))
  
  trajs <- trajs[-all_na_inds,]
  dates <- dates[-all_na_inds,]
  
  # pivot table
  
  trajs <- as.data.frame(trajs)
  
  colnames(trajs) <- paste0(1:ncol(trajs))
  
  trajs_diags <- trajs %>%
    pivot_longer(everything(), names_to = "diag_order", values_to = "diag")
  
  colnames(dates) <- paste0(1:ncol(dates))
  
  dates_diags <- dates %>%
    pivot_longer(everything(), names_to = "diag_order", values_to = "date")
  
  # convert to named 
  
  if(data_type == "icd10_block") {
    trajs_diags <- trajs_diags %>% mutate(name = blocks$name[match(diag, blocks$block)])
  } else if(data_type == "icd10_chapter") {
    trajs_diags <- trajs_diags %>% mutate(name = blocks$chapter_name[match(diag, blocks$chapter)])
  }
  
  # join with dates
  
  if(type == "day_interval") {
    trajs_diags <- cbind(trajs_diags, dates_diags)
    trajs_diags <- trajs_diags[, -4]
  }
  
  # remove NAs
  
  trajs_diags <- trajs_diags %>%
    filter(!is.na(diag))
  
  # plot distribution
  
  dis_tots <- trajs_diags %>%
    group_by(name) %>%
    count()
  
  trajs_diags <- trajs_diags %>%
    full_join(dis_tots, by = "name")
  
  if(type == "diag_order") {
    colnames(trajs_diags)[4] <- "total"
  } else if(type == "day_interval") {
    colnames(trajs_diags)[5] <- "total"
  }
  
  if(type == "diag_order") {
    trajs_diag_cnts <- trajs_diags %>%
      group_by(diag_order) %>%
      count(name)
  } else if(type == "day_interval") {
    
    # trajs_diags$interval <- cut(trajs_diags$date, breaks = seq(-1, 7801, 200), labels = FALSE)
    # 
    # trajs_diag_cnts <- trajs_diags %>%
    #   group_by(interval) %>%
    #   count(name)
    trajs_diags %<>%
      mutate(date = date/365) %>%
      group_by(diag) %>%
      summarise(pct_25 = quantile(date, 0.25), pct_50 = median(date), pct_75 = quantile(date, 0.75))
  }
  
  trajs_diags <- distinct(trajs_diags)
  
  # create abbreviation
  
  if(gender == "f") {
    abbrev <- c("D60-D64" = "Aplas.Oth.Anaem.",
                "M15-M19" = "Arthrosis",
                "B95-B98" = "Bact.Vir.Oth.Infect.",
                "D10-D36" = "Beni.Neop.",
                "I60-I69" = "Cerebro.Dis.",
                "J40-J47" = "Chron.Low.Resp.",
                "M40-M43" = "Deform.Dorso.",
                "E10-E14" = "Diab.Mell.",
                "I70-I79" = "Dis.Arter.Capill.",
                "K70-K77" = "Dis.Liver",
                "K20-K31" = "Dis.Oeso.Stom.Duod.",
                "K00-K14" = "Dis.Oral.Saliv.Jaws",
                "K65-K67" = "Dis.Peritoneum",
                "I80-I89" = "Dis.Vein.Lymph.",
                "M80-M85" = "Dis.Bone.Struc",
                "N60-N64" = "Dis.Breast",
                "H30-H36" = "Dis.Chor.Retina",
                "H00-H06" = "Dis.Eyelid.Lacr.Orb",
                "K80-K87" = "Dis.Gallb.Bili.Panc.",
                "H25-H28" = "Dis.Lens",
                "L60-L75" = "Dis.Skin.App.",
                "M65-M68" = "Dis.Synov.Tendon",
                "E00-E07" = "Dis.Thyroid.Gla.",
                "G40-G47" = "Episod.Paroxy.Dis.",
                "H40-H42" = "Glaucoma",
                "K40-K46" = "Hernia",
                "I10-I15" = "Hypertensive.Dis.",
                "D00-D09" = "In.Situ.Neopl.",
                "L00-L08" = "Infect.Skin.Subcut.Tiss.",
                "N70-N77" = "Inflamm.Dis.Fem.Pelvic.Org.",
                "M05-M14" = "Inflamm.Polyarthrop.",
                "J09-J18" = "Influ.Pneumon.",
                "A00-A09" = "Intest.Infect.Dis.",
                "I20-I25" = "Isch.Heart.Dis.",
                "C50-C50" = "Malig.Neopl.Breast",
                "C15-C26" = "Malig.Neopl.Digest.",
                "C51-C58" = "Malig.Neopl.Fem.Genit.",
                "C76-C80" = "Malig.Neopl.Ill-Def.Second.",
                "C43-C44" = "Melano.Oth.Malig.Neopl.Skin",
                "F10-F19" = "Mental.Behav.Dis.Psycho.Subst.",
                "E70-E90" = "Metabolic.Dis.",
                "F30-F39" = "Mood [Affective] Dis.",
                "D37-D48" = "Neopl.Uncert.Behav.",
                "G50-G59" = "Nerve.Plexus.Dis.",
                "F40-F48" = "Neurot.Stress.Somato.Dis.",
                "K50-K52" = "Noninfect.Enterit.Colit.",
                "N80-N98" = "Noninflamm.Dis.Fem.Genit.",
                "D50-D53" = "Nutrition.Anaem.",
                "E65-E68" = "Obesity.Oth.Hyperalim.",
                "J20-J22" = "Oth.Acute.Low.Resp.Infect.",
                "I95-I99" = "Oth.Dis.Circulat.Sys.",
                "A30-A49" = "Oth.Bacter.Dis.",
                "D70-D77" = "Oth.Dis.Blood.Org.",
                "K55-K64" = "Oth.Dis.Intestin.",
                "J90-J94" = "Oth.Dis.Pleura",
                "K90-K93" = "Oth.Dis.Digestive.Sys.",
                "J95-J99" = "Oth.Dis.Respir.Sys.",
                "J30-J39" = "Oth.Dis.Upper.Respir.",
                "N30-N39" = "Oth.Dis.Urinary.Sys.",
                "H90-H95" = "Oth.Dis.Ear",
                "G90-G99" = "Oth.Dis.Nervous.Sys.",
                "L80-L99" = "Oth.Dis.Skin.Subcutan.",
                "M50-M54" = "Oth.Dorsopathies",
                "I30-I52" = "Oth.Heart.Dis.",
                "M20-M25" = "Oth.Joint.Dis.",
                "M70-M79" = "Oth.Soft.Tiss.Dis.",
                "I26-I28" = "Pulmon.Heart.Dis.Pulmon.Circ.",
                "N17-N19" = "Renal Failure",
                "N10-N16" = "Renal.Tubulo-Interst.Dis.",
                "M45-M49" = "Spondylopathies",
                "N20-N23" = "Urolithiasis",
                "H53-H54" = "Vis.Disturb.Blind.")
  } else if(gender == "m") {
    abbrev <- c("D60-D64" = "Aplas.Oth.Anaem.",
                "M15-M19" = "Arthrosis",
                "B95-B98" = "Bact.Vir.Oth.Infect.",
                "D10-D36" = "Beni.Neop.",
                "I60-I69" = "Cerebro.Dis.",
                "J40-J47" = "Chron.Low.Resp.",
                "E10-E14" = "Diab.Mell.",
                "I70-I79" = "Dis.Arter.Capill.",
                "K70-K77" = "Dis.Liver",
                "N40-N51" = "Dis.Male.Genit.",
                "K20-K31" = "Dis.Oeso.Stom.Duod.",
                "K00-K14" = "Dis.Oral.Saliv.Jaws",
                "K65-K67" = "Dis.Peritoneum",
                "I80-I89" = "Dis.Vein.Lymph.",
                "M80-M85" = "Dis.Bone.Struc",
                "H30-H36" = "Dis.Chor.Retina",
                "H00-H06" = "Dis.Eyelid.Lacr.Orb",
                "K80-K87" = "Dis.Gallb.Bili.Panc.",
                "H25-H28" = "Dis.Lens",
                "L60-L75" = "Dis.Skin.App.",
                "M65-M68" = "Dis.Synov.Tendon",
                "E00-E07" = "Dis.Thyroid.Gla.",
                "G40-G47" = "Episod.Paroxy.Dis.",
                "H40-H42" = "Glaucoma",
                "K40-K46" = "Hernia",
                "I10-I15" = "Hypertensive.Dis.",
                "L00-L08" = "Infect.Skin.Subcut.Tiss.",
                "M05-M14" = "Inflamm.Polyarthrop.",
                "J09-J18" = "Influ.Pneumon.",
                "A00-A09" = "Intest.Infect.Dis.",
                "I20-I25" = "Isch.Heart.Dis.",
                "C15-C26" = "Malig.Neopl.Digest.",
                "C76-C80" = "Malig.Neopl.Ill-Def.Second.",
                "C60-C63" = "Malig.Neopl.Male.Genit.",
                "C64-C68" = "Malig.Neopl.Urinary.",
                "C81-C96" = "Malig.Neopl.Lymph.Haemato.",
                "C43-C44" = "Melano.Oth.Malig.Neopl.Skin",
                "F10-F19" = "Mental.Behav.Dis.Psycho.Subst.",
                "E70-E90" = "Metabolic.Dis.",
                "F30-F39" = "Mood [Affective] Dis.",
                "D37-D48" = "Neopl.Uncert.Behav.",
                "G50-G59" = "Nerve.Plexus.Dis.",
                "F40-F48" = "Neurot.Stress.Somato.Dis.",
                "K50-K52" = "Noninfect.Enterit.Colit.",
                "D50-D53" = "Nutrition.Anaem.",
                "E65-E68" = "Obesity.Oth.Hyperalim.",
                "J20-J22" = "Oth.Acute.Low.Resp.Infect.",
                "I95-I99" = "Oth.Dis.Circulat.Sys.",
                "A30-A49" = "Oth.Bacter.Dis.",
                "D70-D77" = "Oth.Dis.Blood.Org.",
                "K55-K64" = "Oth.Dis.Intestin.",
                "J90-J94" = "Oth.Dis.Pleura",
                "K90-K93" = "Oth.Dis.Digestive.Sys.",
                "J95-J99" = "Oth.Dis.Respir.Sys.",
                "J30-J39" = "Oth.Dis.Upper.Respir.",
                "N30-N39" = "Oth.Dis.Urinary.Sys.",
                "H90-H95" = "Oth.Dis.Ear",
                "N25-N29" = "Oth.Dis.Kidney.Uret.",
                "G90-G99" = "Oth.Dis.Nervous.Sys.",
                "L80-L99" = "Oth.Dis.Skin.Subcutan.",
                "M50-M54" = "Oth.Dorsopathies",
                "I30-I52" = "Oth.Heart.Dis.",
                "M20-M25" = "Oth.Joint.Dis.",
                "M70-M79" = "Oth.Soft.Tiss.Dis.",
                "I26-I28" = "Pulmon.Heart.Dis.Pulmon.Circ.",
                "N17-N19" = "Renal Failure",
                "N10-N16" = "Renal.Tubulo-Interst.Dis.",
                "M45-M49" = "Spondylopathies",
                "N20-N23" = "Urolithiasis",
                "H53-H54" = "Vis.Disturb.Blind.")
  }
  
  # get chapter groups
  
  if(type == "diag_order") {
    trajs_diags$diag_order <- as.integer(trajs_diags$diag_order)
    trajs_diag_cnts$diag_order <- as.integer(trajs_diag_cnts$diag_order)
    trajs_cumfreq <- trajs_diags %>%
      inner_join(trajs_diag_cnts, by = c("name", "diag_order")) %>%
      group_by(name) %>%
      arrange(name, diag_order) %>%
      summarise(diag_order = diag_order, cumfreq = cumsum(n)/total)
    block_cumfreq_qnts <- trajs_cumfreq %>%
      group_by(name) %>%
      summarise(pct_25 = min(which(cumfreq >= 0.25)), pct_50 = min(which(cumfreq >= 0.5)), pct_75 = min(which(cumfreq >= 0.75)))
    
    if(data_type == "icd10_block") {
      chapt_grps <- blocks$chapter[match(block_cumfreq_qnts$name, blocks$name)]
      chapt_grps <- factor(chapt_grps, levels = sort(unique(chapt_grps)))
      block_cumfreq_qnts$name <- abbrev
    } else if(data_type == "icd10_chapter") {
      block_cumfreq_qnts$name <- str_to_title(block_cumfreq_qnts$name)
    }
    
    block_order <- block_cumfreq_qnts$name[order(block_cumfreq_qnts$pct_50)]
    
    block_cumfreq_qnts$name <- factor(block_cumfreq_qnts$name, levels = block_order)
    
  } else if(type == "day_interval") {
    # trajs_diags$interval <- as.integer(trajs_diags$interval)
    # trajs_diag_cnts$interval <- as.integer(trajs_diag_cnts$interval)
    # trajs_cumfreq <- trajs_diags %>%
    #   inner_join(trajs_diag_cnts, by = c("name", "interval")) %>%
    #   arrange(name, interval) %>%
    #   select(name, interval, n, total) %>%
    #   distinct() %>%
    #   group_by(name) %>%
    #   summarise(interval = interval, cumfreq =  cumsum(n)/total)
    # block_cumfreq_qnts <- trajs_cumfreq %>%
    #   group_by(name) %>%
    #   summarise(pct_25 = min(which(cumfreq >= 0.25)), pct_50 = min(which(cumfreq >= 0.5)), pct_75 = min(which(cumfreq >= 0.75)))
    # 
    
    if(data_type == "icd10_block") {
      # chapt_grps <- blocks$chapter[match(block_cumfreq_qnts$name, blocks$name)]
      # block_cumfreq_qnts$name <- abbrev
      trajs_diags$name <- abbrev[match(trajs_diags$diag, names(abbrev))]
      block_cumfreq_qnts <- trajs_diags
      block_cumfreq_qnts %<>% select(name, pct_25, pct_50, pct_75)
    } else if(data_type == "icd10_chapter") {
      # block_cumfreq_qnts$name <- str_to_title(block_cumfreq_qnts$name)
    }
    
    block_order <- block_cumfreq_qnts$name[order(block_cumfreq_qnts$pct_50)]
    block_cumfreq_qnts$name <- factor(block_cumfreq_qnts$name, levels = block_order)
    
  }
  if(data_type == "icd10_block") {
    # block_cumfreq_qnts$chapter <- factor(chapt_grps, levels = sort(unique(chapt_grps)))
  }
  
  block_cumfreq_qnts$sex <- gender
  block_cumfreq_qnts$data_type <- data_type
  block_cumfreq_qnts$distrib_type <- type
  
  return(block_cumfreq_qnts)
}


add_nl <- function(x) {
  x_sub <- substring(x, 50)
  spc_loc <- str_locate(x_sub, " ")
  spc_loc <- spc_loc[, 2]
  x_sub2 <- paste0(substring(x, 1, 49), substring(x_sub, 1, spc_loc), "/n", substring(x_sub, spc_loc+1), collapse = "")
  return(x_sub2)
}

get_bf_distrib <- function(df, coef, x_var, y_var) {
  # get x values
  x_vals <- df %>% pull(x_var)
  # get y coordinates of best fit line
  df_line_y <- x_vals*coef
  # choose arbitrary point to calculate orthonormal basis vector
  df_pt_y <- df_line_y[50]
  df_pt_x <- x_vals[50]
  # calculate orthonormal basis vector
  df_orth_norm <- (1/sqrt(df_pt_x*df_pt_x + df_pt_y*df_pt_y))*c(df_pt_x, df_pt_y)
  # get coordinates of points on best fit line
  df_line_pt <- vector(length = nrow(df))
  for(i in 1:length(df_line_pt)) {
    df_line_pt[i] <- c(unlist(df[i, x_var]), unlist(df[i, y_var])) %*% df_orth_norm
  }
  return(df_line_pt)
}

get_uncentre <- function(mean_x, mean_y, coef) {
  return(mean_y - (coef*mean_x))
}

# diagnosis order, male block

# male trajectory matrix
traj_file = "male_disease_block_1pct_traj.tsv"
# male time to diagnosis
traj_dates_file = "male_diag_block_1pct_norm_dates.tsv"

male_block_diag_ord <- get_distrib("diag_order", "icd10_block", "m")

# diagnosis order, female block

# female trajectory matrix
traj_file = "fem_disease_block_1pct_traj.tsv"
# female time to diagnosis
traj_dates_file = "fem_diag_block_1pct_norm_dates.tsv"

fem_block_diag_ord <- get_distrib("diag_order", "icd10_block", "f")

# time to diagnosis, male block

# male trajectory matrix
traj_file = "male_disease_block_1pct_traj.tsv"
# male time to diagnosis matrix
traj_dates_file = "male_diag_block_1pct_norm_dates.tsv"

male_block_day_interv <- get_distrib("day_interval", "icd10_block", "m")

# time to diagnosis, female block

# female trajectory matrix
traj_file = "fem_disease_block_1pct_traj.tsv"
# female time to diagnosis matrix
traj_dates_file = "fem_diag_block_1pct_norm_dates.tsv"

fem_block_day_interv <- get_distrib("day_interval", "icd10_block", "f")

# combine

fm_distrib <- rbind(male_block_diag_ord,
                    male_block_day_interv,
                    fem_block_diag_ord,
                    fem_block_day_interv)

# plot of male median time to diagnosis against median diagnosis order

# get pre-centred means
male_means <- fm_distrib %>%
  select(-pct_25, -pct_75) %>%
  pivot_wider(names_from = distrib_type, values_from = pct_50) %>%
  filter(sex == "m") %>%
  summarise(mean_diag_order = mean(diag_order), mean_day_interval = mean(day_interval)) %>%
  as.matrix()

# create temporary centred data
male_tmp <- fm_distrib %>%
  select(-pct_25, -pct_75) %>%
  pivot_wider(names_from = distrib_type, values_from = pct_50) %>%
  filter(sex == "m") %>%
  mutate(diag_order = scale(diag_order, scale = F, center = T), day_interval = scale(day_interval, scale = F, center = T))

# get line of best fit
male_bf <- lm(day_interval ~ diag_order, data = male_tmp)
# get best fit coefficients
male_bf_coef <- coef(male_bf)
# get perpendicular line coefficients
male_perp_coef <- -(1/male_bf_coef[2])

male_line_pt <- get_bf_distrib(male_tmp, male_bf_coef[2], "diag_order", "day_interval")
male_perp_pt <- get_bf_distrib(male_tmp, male_perp_coef, "diag_order", "day_interval")

# calculate SD of line
male_line_sd <- sd(male_line_pt)
male_perp_sd <- sd(male_perp_pt)

# uncentring scalar

uncentre <- get_uncentre(male_means[1], male_means[2], male_bf_coef[2])
uncentre2 <- get_uncentre(male_means[1], male_means[2], male_perp_coef)

# get uncentred coordinates on bf line

male_uncentred_pt <- male_line_pt/cos(atan(male_perp_coef))+uncentre2
male_uncentred_perp <- male_perp_pt/cos(atan(male_bf_coef[2]))+uncentre

plot <- fm_distrib %>%
  select(-pct_25, -pct_75) %>%
  pivot_wider(names_from = distrib_type, values_from = pct_50) %>%
  filter(sex == "m") %>%
  ggplot(aes(x = diag_order, y = day_interval)) +
  geom_point(size = 3) +
  annotate("text", x = 6.5, y = 2.5, label = paste0("Beta = ", round(male_bf_coef[2], 2)), size = 15) +
  geom_abline(slope = male_bf_coef[2], intercept = uncentre) +
  # geom_abline(slope = male_perp_coef, intercept = (male_line_sd/cos(atan(male_perp_coef)))+uncentre2, lty = 2) +
  # geom_abline(slope = male_perp_coef, intercept = -(male_line_sd/cos(atan(male_perp_coef)))+uncentre2, lty = 2) +
  # geom_abline(slope = male_perp_coef, intercept = (2*male_line_sd/cos(atan(male_perp_coef)))+uncentre2, lty = 3) +
  # geom_abline(slope = male_perp_coef, intercept = (-2*male_line_sd/cos(atan(male_perp_coef)))+uncentre2, lty = 3) +
  xlab("Median position of diagnosis\nin male diagnoses records") +
  ylab("Median time to diagnosis\nin male diagnoses records (years)") +
  scale_x_continuous(breaks = seq(2, 9, 3)) +
  scale_y_continuous(breaks = seq(0, 9, 3)) +
  theme_bw() +
  theme(text = element_text(size = 25)) +
  coord_cartesian(xlim = c(2, 9))
ggsave(filename = "male_diag_order_against_day_interval_scatter.png",
       units = "px", width = 2000, height = 2000,
       plot = plot)

# create result df

male_df <- fm_distrib %>%
  select(-pct_25, -pct_75) %>%
  pivot_wider(names_from = distrib_type, values_from = pct_50) %>%
  filter(sex == "m") %>%
  filter(!is.na(diag_order) & !is.na(day_interval)) %>%
  mutate(bf_sd_thresh = "within 1sd",
         perp_sd_thresh = "within 1sd")

male_df[male_line_pt > male_line_sd, "bf_sd_thresh"] <- ">1sd"
male_df[male_line_pt > male_line_sd*2, "bf_sd_thresh"] <- ">2sd"
male_df[male_line_pt < -male_line_sd, "bf_sd_thresh"] <- "<-1sd"
male_df[male_line_pt < -male_line_sd*2, "bf_sd_thresh"] <- "<-2sd"

male_df[male_perp_pt > male_perp_sd, "perp_sd_thresh"] <- ">1sd"
male_df[male_perp_pt > male_perp_sd*2, "perp_sd_thresh"] <- ">2sd"
male_df[male_perp_pt < -male_perp_sd, "perp_sd_thresh"] <- "<-1sd"
male_df[male_perp_pt < -male_perp_sd*2, "perp_sd_thresh"] <- "<-2sd"


# plot of female day interval against diag order

# get pre-centred means
fem_means <- fm_distrib %>%
  select(-pct_25, -pct_75) %>%
  pivot_wider(names_from = distrib_type, values_from = pct_50) %>%
  filter(sex == "f") %>%
  summarise(mean_diag_order = mean(diag_order), mean_day_interval = mean(day_interval)) %>%
  as.matrix()

# create temporary centred data
fem_tmp <- fm_distrib %>%
  select(-pct_25, -pct_75) %>%
  pivot_wider(names_from = distrib_type, values_from = pct_50) %>%
  filter(sex == "f") %>%
  mutate(diag_order = scale(diag_order, scale = F, center = T), day_interval = scale(day_interval, scale = F, center = T))

# get line of best fit
fem_bf <- lm(day_interval ~ diag_order, data = fem_tmp)
# get best fit coefficients
fem_bf_coef <- coef(fem_bf)
# get perpendicular line coefficients
fem_perp_coef <- -(1/fem_bf_coef[2])

fem_line_pt <- get_bf_distrib(fem_tmp, fem_bf_coef[2], "diag_order", "day_interval")
fem_perp_pt <- get_bf_distrib(fem_tmp, fem_perp_coef, "diag_order", "day_interval")

# calculate SD of line
fem_line_sd <- sd(fem_line_pt)
fem_perp_sd <- sd(fem_perp_pt)

# uncentring scalar

uncentre <- get_uncentre(fem_means[1], fem_means[2], fem_bf_coef[2])
uncentre2 <- get_uncentre(fem_means[1], fem_means[2], fem_perp_coef)

# get uncentred coordinates on bf line

fem_uncentred_pt <- fem_line_pt/cos(atan(fem_perp_coef))+uncentre2
fem_uncentred_perp <- fem_perp_pt/cos(atan(fem_bf_coef[2]))+uncentre


plot <- fm_distrib %>%
  select(-pct_25, -pct_75) %>%
  pivot_wider(names_from = distrib_type, values_from = pct_50) %>%
  filter(sex == "f") %>%
  ggplot(aes(x = as.integer(diag_order), y = day_interval)) +
  geom_point(size = 3) +
  annotate("text", x = 7, y = 3, label = paste0("Beta = ", round(fem_bf_coef[2], 2)), size = 15) +
  geom_abline(slope = fem_bf_coef[2], intercept = uncentre) +
  # geom_abline(slope = fem_perp_coef, intercept = (fem_line_sd/cos(atan(fem_perp_coef)))+uncentre2, lty = 2) +
  # geom_abline(slope = fem_perp_coef, intercept = -(fem_line_sd/cos(atan(fem_perp_coef)))+uncentre2, lty = 2) +
  # geom_abline(slope = fem_perp_coef, intercept = (2*fem_line_sd/cos(atan(fem_perp_coef)))+uncentre2, lty = 3) +
  # geom_abline(slope = fem_perp_coef, intercept = (-2*fem_line_sd/cos(atan(fem_perp_coef)))+uncentre2, lty = 3) +
  xlab("Median position of diagnosis\nin female diagnoses records") +
  ylab("Median time to diagnosis\nin female diagnoses records (years)") +
  scale_x_continuous(breaks = seq(2, 9, 3)) +
  scale_y_continuous(breaks = seq(0, 11, 3)) +
  theme_bw() +
  theme(text = element_text(size = 25))
ggsave(filename = "fem_diag_order_against_day_interval_scatter.png",
       units = "px", width = 2000, height = 2000,
       plot = plot)

# create result df

fem_df <- fm_distrib %>%
  select(-pct_25, -pct_75) %>%
  pivot_wider(names_from = distrib_type, values_from = pct_50) %>%
  filter(sex == "f") %>%
  filter(!is.na(diag_order) & !is.na(day_interval)) %>%
  mutate(bf_sd_thresh = "within 1sd",
         perp_sd_thresh = "within 1sd")

fem_df[fem_line_pt > fem_line_sd, "bf_sd_thresh"] <- ">1sd"
fem_df[fem_line_pt > fem_line_sd*2, "bf_sd_thresh"] <- ">2sd"
fem_df[fem_line_pt < -fem_line_sd, "bf_sd_thresh"] <- "<-1sd"
fem_df[fem_line_pt < -fem_line_sd*2, "bf_sd_thresh"] <- "<-2sd"

fem_df[fem_perp_pt > fem_perp_sd, "perp_sd_thresh"] <- ">1sd"
fem_df[fem_perp_pt > fem_perp_sd*2, "perp_sd_thresh"] <- ">2sd"
fem_df[fem_perp_pt < -fem_perp_sd, "perp_sd_thresh"] <- "<-1sd"
fem_df[fem_perp_pt < -fem_perp_sd*2, "perp_sd_thresh"] <- "<-2sd"

# combine

fm_df <- rbind(male_df,
               fem_df)

# read in blocks

blocks <- read_tsv(icd10_blocks_file)

# get codes for abbreviations

fem_abbrev <- c("D60-D64" = "Aplas.Oth.Anaem.",
                "M15-M19" = "Arthrosis",
                "B95-B98" = "Bact.Vir.Oth.Infect.",
                "D10-D36" = "Beni.Neop.",
                "I60-I69" = "Cerebro.Dis.",
                "J40-J47" = "Chron.Low.Resp.",
                "M40-M43" = "Deform.Dorso.",
                "E10-E14" = "Diab.Mell.",
                "I70-I79" = "Dis.Arter.Capill.",
                "K70-K77" = "Dis.Liver",
                "K20-K31" = "Dis.Oeso.Stom.Duod.",
                "K00-K14" = "Dis.Oral.Saliv.Jaws",
                "K65-K67" = "Dis.Peritoneum",
                "I80-I89" = "Dis.Vein.Lymph.",
                "M80-M85" = "Dis.Bone.Struc",
                "N60-N64" = "Dis.Breast",
                "H30-H36" = "Dis.Chor.Retina",
                "H00-H06" = "Dis.Eyelid.Lacr.Orb",
                "K80-K87" = "Dis.Gallb.Bili.Panc.",
                "H25-H28" = "Dis.Lens",
                "H49-H52" = "Dis.Ocul.Binoc.Accom.Refr.",
                "H15-H22" = "Dis.Scle.Corn.Iris.Cili.",
                "L60-L75" = "Dis.Skin.App.",
                "M65-M68" = "Dis.Synov.Tendon",
                "E00-E07" = "Dis.Thyroid.Gla.",
                "G40-G47" = "Episod.Paroxy.Dis.",
                "H40-H42" = "Glaucoma",
                "K40-K46" = "Hernia",
                "I10-I15" = "Hypertensive.Dis.",
                "D00-D09" = "In.Situ.Neopl.",
                "L00-L08" = "Infect.Skin.Subcut.Tiss.",
                "N70-N77" = "Inflamm.Dis.Fem.Pelvic.Org.",
                "M05-M14" = "Inflamm.Polyarthrop.",
                "J09-J18" = "Influ.Pneumon.",
                "A00-A09" = "Intest.Infect.Dis.",
                "I20-I25" = "Isch.Heart.Dis.",
                "C50-C50" = "Malig.Neopl.Breast",
                "C15-C26" = "Malig.Neopl.Digest.",
                "C69-C72" = "Malig.Neopl.Eye.Brain.Nerv.",
                "C51-C58" = "Malig.Neopl.Fem.Genit.",
                "C76-C80" = "Malig.Neopl.Ill-Def.Second.",
                "C43-C44" = "Melano.Oth.Malig.Neopl.Skin",
                "F10-F19" = "Mental.Behav.Dis.Psycho.Subst.",
                "E70-E90" = "Metabolic.Dis.",
                "F30-F39" = "Mood [Affective] Dis.",
                "D37-D48" = "Neopl.Uncert.Behav.",
                "G50-G59" = "Nerve.Plexus.Dis.",
                "F40-F48" = "Neurot.Stress.Somato.Dis.",
                "K50-K52" = "Noninfect.Enterit.Colit.",
                "N80-N98" = "Noninflamm.Dis.Fem.Genit.",
                "D50-D53" = "Nutrition.Anaem.",
                "E65-E68" = "Obesity.Oth.Hyperalim.",
                "J20-J22" = "Oth.Acute.Low.Resp.Infect.",
                "I95-I99" = "Oth.Dis.Circulat.Sys.",
                "A30-A49" = "Oth.Bacter.Dis.",
                "D70-D77" = "Oth.Dis.Blood.Org.",
                "K55-K64" = "Oth.Dis.Intestin.",
                "J90-J94" = "Oth.Dis.Pleura",
                "K90-K93" = "Oth.Dis.Digestive.Sys.",
                "J95-J99" = "Oth.Dis.Respir.Sys.",
                "J30-J39" = "Oth.Dis.Upper.Respir.",
                "N30-N39" = "Oth.Dis.Urinary.Sys.",
                "H90-H95" = "Oth.Dis.Ear",
                "G90-G99" = "Oth.Dis.Nervous.Sys.",
                "L80-L99" = "Oth.Dis.Skin.Subcutan.",
                "M50-M54" = "Oth.Dorsopathies",
                "I30-I52" = "Oth.Heart.Dis.",
                "M20-M25" = "Oth.Joint.Dis.",
                "M86-M90" = "Oth.Osteopathies",
                "M70-M79" = "Oth.Soft.Tiss.Dis.",
                "I26-I28" = "Pulmon.Heart.Dis.Pulmon.Circ.",
                "N17-N19" = "Renal Failure",
                "N10-N16" = "Renal.Tubulo-Interst.Dis.",
                "M45-M49" = "Spondylopathies",
                "N20-N23" = "Urolithiasis",
                "H53-H54" = "Vis.Disturb.Blind.")

male_abbrev <- c("Aplas.Oth.Anaem.",
                 "Arthrosis",
                 "Bact.Vir.Oth.Infect.",
                 "Beni.Neop.",
                 "Cerebro.Dis.",
                 "Chron.Low.Resp.",
                 "Diab.Mell.",
                 "Dis.Arter.Capill.",
                 "Dis.Liver",
                 "Dis.Male.Genit.",
                 "Dis.Oeso.Stom.Duod.",
                 "Dis.Oral.Saliv.Jaws",
                 "Dis.Peritoneum",
                 "Dis.Vein.Lymph.",
                 "Dis.Bone.Struc",
                 "Dis.Chor.Retina",
                 "Dis.Eyelid.Lacr.Orb",
                 "Dis.Gallb.Bili.Panc.",
                 "Dis.Lens",
                 "Dis.Ocul.Binoc.Accom.Refr.",
                 "Dis.Scle.Corn.Iris.Cili.",
                 "Dis.Skin.App.",
                 "Dis.Synov.Tendon",
                 "Dis.Thyroid.Gla.",
                 "Episod.Paroxy.Dis.",
                 "Glaucoma",
                 "Hernia",
                 "Hypertensive.Dis.",
                 "Infect.Skin.Subcut.Tiss.",
                 "Inflamm.Polyarthrop.",
                 "Influ.Pneumon.",
                 "Intest.Infect.Dis.",
                 "Isch.Heart.Dis.",
                 "Malig.Neopl.Digest.",
                 "Malig.Neopl.Eye.Brain.Nerv.",
                 "Malig.Neopl.Ill-Def.Second.",
                 "Malig.Neopl.Male.Genit.",
                 "Malig.Neopl.Urinary.",
                 "Malig.Neopl.Lymph.Haemato.",
                 "Melano.Oth.Malig.Neopl.Skin",
                 "Mental.Behav.Dis.Psycho.Subst.",
                 "Metabolic.Dis.",
                 "Mood [Affective] Dis.",
                 "Neopl.Uncert.Behav.",
                 "Nerve.Plexus.Dis.",
                 "Neurot.Stress.Somato.Dis.",
                 "Noninfect.Enterit.Colit.",
                 "Nutrition.Anaem.",
                 "Obesity.Oth.Hyperalim.",
                 "Oth.Acute.Low.Resp.Infect.",
                 "Oth.Dis.Circulat.Sys.",
                 "Oth.Bacter.Dis.",
                 "Oth.Dis.Blood.Org.",
                 "Oth.Dis.Intestin.",
                 "Oth.Dis.Pleura",
                 "Oth.Dis.Digestive.Sys.",
                 "Oth.Dis.Respir.Sys.",
                 "Oth.Dis.Upper.Respir.",
                 "Oth.Dis.Urinary.Sys.",
                 "Oth.Dis.Ear",
                 "Oth.Dis.Kidney.Uret.",
                 "Oth.Dis.Nervous.Sys.",
                 "Oth.Dis.Skin.Subcutan.",
                 "Oth.Dorsopathies",
                 "Oth.Heart.Dis.",
                 "Oth.Joint.Dis.",
                 "Oth.Osteopathies",
                 "Oth.Soft.Tiss.Dis.",
                 "Pulmon.Heart.Dis.Pulmon.Circ.",
                 "Renal Failure",
                 "Renal.Tubulo-Interst.Dis.",
                 "Spondylopathies",
                 "Urolithiasis",
                 "Vis.Disturb.Blind.")

male_uni <- male_abbrev[which(!male_abbrev %in% fem_abbrev)]
male_abbrev2 <- fem_abbrev[-which(!fem_abbrev %in% male_abbrev)]
names(male_abbrev2) <- names(fem_abbrev)[match(male_abbrev2, fem_abbrev)]
male_abbrev2 <- c(male_abbrev2, male_uni)
names(male_abbrev2)[names(male_abbrev2) == ""] <- c("N40-N51",
                                                    "C60-C63",
                                                    "C64-C68",
                                                    "C81-C96",
                                                    "N25-N29")

fem_abbrev_df <- as.data.frame(fem_abbrev)
fem_code <- names(fem_abbrev)
fem_abbrev_df %<>%
  mutate(block = fem_code)
colnames(fem_abbrev_df) <- c("abbrev", "block")

male_abbrev_df <- as.data.frame(male_abbrev2)
male_code <- names(male_abbrev2)
male_abbrev_df %<>%
  mutate(block = male_code)
colnames(male_abbrev_df) <- c("abbrev", "block")

comb_df <- rbind(fem_abbrev_df,
                 male_abbrev_df)
comb_df <- comb_df[!duplicated(comb_df),]

# clean

fm_df %<>%
  select(-data_type)


fm_df <- cbind(code = comb_df$block[match(fm_df$name, comb_df$abbrev)], fm_df)
colnames(fm_df) <- c("Code", "Name", "Sex", "Diagnosis order median", "Day interval median", "Best fit SD threshold", "Deviation SD threshold")
fem_out <- fm_df %>% filter(Sex == "f") %>% filter(`Best fit SD threshold` != "within 1sd" | `Deviation SD threshold` != "within 1sd") %>% select(-Sex)
male_out <- fm_df %>% filter(Sex == "m") %>% filter(`Best fit SD threshold` != "within 1sd" | `Deviation SD threshold` != "within 1sd") %>% select(-Sex)

write_tsv(fem_out, file = "fem_diag_order_vs_day_interval_sd_thresh_table.tsv")
write_tsv(male_out, file = "male_diag_order_vs_day_interval_sd_thresh_table.tsv")
