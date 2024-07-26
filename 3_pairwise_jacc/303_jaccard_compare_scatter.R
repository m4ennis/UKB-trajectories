# scatterplot comparing male/female jaccard values

library(tidyverse)
library(magrittr)
library(sf)
library(sp)  
library(ggrepel)

# Variables

# female pairwise Jaccard coincidence matrix
fem_mm_jacc_file = "ukb_fem_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_table.tsv"

# male pairwise Jaccard coincidence matrix
male_mm_jacc_file = "ukb_male_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_table.tsv"

# ICD10 block dictionary
icd10_blocks_file = "icd10_blocks.tsv"

# block or chapter
icd10_type = "block"

# functions

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

# Read in multimorbidity pairwise matrix  

fem_jacc <- read_tsv(fem_mm_jacc_file)
fem_jacc <- as.matrix(fem_jacc)

male_jacc <- read_tsv(male_mm_jacc_file)
male_jacc <- as.matrix(male_jacc)

if(icd10_type == "block") {
  # report top sex-specific combinations
  
  dis_use <- intersect(colnames(fem_jacc), colnames(male_jacc))
  
  fem_spec <- fem_jacc
  
  fem_spec %<>%
    as.data.frame() %>%
    mutate(comp1 = colnames(fem_spec)) %>%
    pivot_longer(-comp1, names_to = "comp2", values_to = "fem_jacc")
  
  male_spec <- male_jacc
  
  male_spec %<>%
    as.data.frame() %>%
    mutate(comp1 = colnames(male_spec)) %>%
    pivot_longer(-comp1, names_to = "comp2", values_to = "male_jacc")
  
  
  fem_spec %<>%
    filter(!comp1 == comp2)
  
  tmp <- fem_spec
  tmp %<>%
    arrange(desc(fem_jacc))
  tmp <- tmp[seq(1, nrow(tmp), 2), ]
  
  fem_spec <- tmp
  
  male_spec %<>%
    filter(!comp1 == comp2) 
  
  tmp <- male_spec
  tmp %<>%
    arrange(desc(male_jacc))
  tmp <- tmp[seq(1, nrow(tmp), 2), ]
  male_spec <- tmp
  
  fem_spec %<>%
    unite(col = "comp", -fem_jacc, sep = "+")
  
  male_spec %<>%
    unite(col = "comp", -male_jacc, sep = "+")
  
  plot <- fem_spec %>%
    filter(!comp %in% male_spec$comp) %>%
    mutate(comp = factor(comp, levels = comp[order(fem_jacc)])) %>%
    arrange(desc(fem_jacc)) %>%
    top_n(5, wt = fem_jacc) %>%
    ggplot(aes(x = fem_jacc, y = comp)) +
    geom_col() +
    theme_bw() +
    xlab("Jaccard coefficient") +
    ylab("Female-specific block comparison") +
    scale_x_continuous(breaks = seq(0, 0.4, 0.1)) +
    coord_cartesian(xlim = c(0, 0.4))
  ggsave(filename = "fem_physio_1pct_block_top_5_fem_spec_pairwise_jaccard_barplot.png",
         plot = plot,
         units = "px",
         width = 1800,
         height = 1600)
  
  plot <- male_spec %>%
    filter(!comp %in% fem_spec$comp) %>%
    mutate(comp = factor(comp, levels = comp[order(male_jacc)])) %>%
    arrange(desc(male_jacc)) %>%
    top_n(5, wt = male_jacc) %>%
    ggplot(aes(x = male_jacc, y = comp)) +
    geom_col() +
    theme_bw() +
    xlab("Jaccard coefficient") +
    ylab("Male-specific block comparison") +
    scale_x_continuous(breaks = seq(0, 0.4, 0.1)) +
    coord_cartesian(xlim = c(0, 0.4))
  ggsave(filename = "male_physio_1pct_block_top_5_male_spec_pairwise_jaccard_barplot.png",
         plot = plot,
         units = "px",
         width = 1800,
         height = 1600)
  
  # remove diseases not present in both
  
  fem_jacc <- fem_jacc[which(colnames(fem_jacc) %in% dis_use), which(colnames(fem_jacc) %in% dis_use)]
  male_jacc <- male_jacc[which(colnames(male_jacc) %in% dis_use), which(colnames(male_jacc) %in% dis_use)]
}

# format

fem_jacc %<>%
  as.data.frame() %>%
  mutate(comp1 = colnames(fem_jacc)) %>%
  pivot_longer(-comp1, names_to = "comp2", values_to = "fem_jacc")

male_jacc %<>%
  as.data.frame() %>%
  mutate(comp1 = colnames(male_jacc)) %>%
  pivot_longer(-comp1, names_to = "comp2", values_to = "male_jacc")

fm_jacc <- fem_jacc %>%
  full_join(male_jacc, by = c("comp1", "comp2"))

# remove duplicates

fm_jacc %<>%
  filter(!comp1 == comp2)

tmp <- fm_jacc
tmp %<>%
  arrange(desc(fem_jacc))
tmp <- tmp[seq(1, nrow(tmp), 2), ]

fm_jacc <- tmp

# calculate distance from line

# get pre-centred means
fm_means <- fm_jacc %>%
  summarise(male_mean = mean(male_jacc), fem_mean = mean(fem_jacc)) %>%
  as.matrix()

# create temporary centred data
fm_tmp <- fm_jacc %>%
  mutate(male_jacc = scale(male_jacc, scale = F), fem_jacc = scale(fem_jacc, scale = F))

# get line of best fit
fm_bf <- lm(fem_jacc ~ male_jacc, data = fm_tmp)
# get best fit coefficients
fm_bf_coef <- coef(fm_bf)
# get perpendicular line coefficients
fm_perp_coef <- -(1/fm_bf_coef[2])

# get best fit coordinates
fm_line_pt <- get_bf_distrib(fm_tmp, fm_bf_coef[2], "male_jacc", "fem_jacc")

# get perpendicular coordinates
fm_line_pt2 <- get_bf_distrib(fm_tmp, fm_perp_coef, "male_jacc", "fem_jacc")

# calculate SD of line
fm_line_sd <- sd(fm_line_pt)

# SD of perpendicular distribution
fm_line_sd2 <- sd(fm_line_pt2)

# uncentring scalar

uncentre <- get_uncentre(fm_means[1], fm_means[2], fm_bf_coef[2])
uncentre2 <- get_uncentre(fm_means[1], fm_means[2], fm_perp_coef)

# add bf distance

fm_jacc$bf_dist <- fm_line_pt

# sd_num <- ifelse(icd10_type == "block", -4, -1)

# plot
if(icd10_type == "block") {
  corr <- cor(fm_jacc$male_jacc, fm_jacc$fem_jacc)
  plot <- fm_jacc %>%
    mutate(label = paste0(comp1, "+", comp2)) %>%
    mutate(label = ifelse(fm_jacc$bf_dist < -5.5*fm_line_sd | abs(fm_line_pt2) > 4.5*fm_line_sd2, label, NA)) %>%
    ggplot(aes(x = male_jacc, y = fem_jacc, label = label)) +
    geom_point(alpha = ifelse(fm_jacc$bf_dist < -5.5*fm_line_sd | abs(fm_line_pt2) > 4.5*fm_line_sd2, 1, 0.05)) +
    geom_text_repel(size = 4, force = 10, max.overlaps = 1000, seed = 42) +
    annotate("text", x = 0.1, y = 0.3, label = paste0("Beta = ", round(fm_bf_coef[2], 2), "\n Corr. = ", round(corr, 2)), size = 10) +
    geom_abline(slope = fm_bf_coef[2], intercept = uncentre) +
    # geom_abline(slope = fm_perp_coef, intercept = (2*fm_line_sd/cos(atan(fm_perp_coef)))+uncentre2, lty = 2) +
    # geom_abline(slope = fm_perp_coef, intercept = -(2*fm_line_sd/cos(atan(fm_perp_coef)))+uncentre2, lty = 2) +
    # geom_abline(slope = fm_perp_coef, intercept = (4*fm_line_sd/cos(atan(fm_perp_coef)))+uncentre2, lty = 3) +
    # geom_abline(slope = fm_perp_coef, intercept = (-4*fm_line_sd/cos(atan(fm_perp_coef)))+uncentre2, lty = 3) +
    # geom_abline(slope = fm_bf_coef[2], intercept = (2*fm_line_sd2/cos(atan(fm_bf_coef[2])))+uncentre, lty = 2) +
    # geom_abline(slope = fm_bf_coef[2], intercept = -(2*fm_line_sd2/cos(atan(fm_bf_coef[2])))+uncentre, lty = 2) +
    # geom_abline(slope = fm_bf_coef[2], intercept = (4*fm_line_sd2/cos(atan(fm_bf_coef[2])))+uncentre, lty = 3) +
    # geom_abline(slope = fm_bf_coef[2], intercept = (-4*fm_line_sd2/cos(atan(fm_bf_coef[2])))+uncentre, lty = 3) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    xlab("Male Jaccard of pairwise ICD10 block") +
    ylab("Female Jaccard of pairwise ICD10 block")
  ggsave("fm_block_pairwise_jaccard_scatter.png", plot = plot,
         units = "px", width = 2000, height = 1800)
  
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
  
  male_abbrev <- c("D60-D64" = "Aplas.Oth.Anaem.",
                   "M15-M19" = "Arthrosis",
                   "B95-B98" = "Bact.Vir.Oth.Infect.",
                   "D10-D36" = "Beni.Neop.",
                   "I60-I69" = "Cerebro.Dis.",
                   "J40-J47" = "Chron.Low.Resp.",
                   "E10-E14" = "Diab.Mell.",
                   "I70-I79" = "Dis.Arter.Capill.",
                   "K20-K31" = "Dis.Liver",
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
  
  # create output df
  
  fm_jacc %<>%
    mutate(bf_sd_thresh = "within -4 SD", perp_devi_sd_thresh = "within +-4 SD")
  
  fm_jacc$bf_sd_thresh[fm_jacc$bf_dist < -4*fm_line_sd] <- "<-4SD"
  fm_jacc$perp_devi_sd_thresh[fm_line_pt2 > 4*fm_line_sd2] <- ">4SD"
  fm_jacc$perp_devi_sd_thresh[fm_line_pt2 < -4*fm_line_sd2] <- "<-4SD"
  
  fm_jacc %<>%
    mutate(comp1_name = NA, comp2_name = NA)
  
  fm_jacc$comp1_name <- comb_df$abbrev[match(fm_jacc$comp1, comb_df$block)]
  fm_jacc$comp2_name <- comb_df$abbrev[match(fm_jacc$comp2, comb_df$block)]
  
  fm_jacc %<>%
    select(comp1, comp2, comp1_name, comp2_name, fem_jacc, male_jacc, bf_sd_thresh, perp_devi_sd_thresh)
  
  colnames(fm_jacc) <- c("Comparison 1", "Comparison 2", "Comparison 1 name", "Comparison 2 name", "Female Jaccard", "Male Jaccard",
                         "Best fit SD", "Perpendicular SD")
  
  
  write_tsv(fm_jacc, "fm_block_pairwise_jaccard_table.tsv")
  
  write_tsv(fm_jacc %>% filter(`Best fit SD` != "within -4 SD" | `Perpendicular SD` != "within +-4 SD"), "fm_block_pairwise_jaccard_table_filt.tsv")
  
  
} else {
  corr <- cor(fm_jacc$male_jacc, fm_jacc$fem_jacc)
  plot <- fm_jacc %>%
    mutate(label = paste0(comp1, "+", comp2)) %>%
    mutate(label = ifelse(fm_jacc$bf_dist < -1*fm_line_sd | abs(fm_line_pt2) > fm_line_sd2, label, NA)) %>%
    ggplot(aes(x = male_jacc, y = fem_jacc, label = label)) +
    geom_point(alpha = ifelse(fm_jacc$bf_dist < -1*fm_line_sd | abs(fm_line_pt2) > fm_line_sd2, 1, 0.2)) +
    geom_text_repel(size = 3, force = 1, seed = 42, max.overlaps = 1000) +
    annotate("text", x = 0.1, y = 0.4, label = paste0("Beta = ", round(fm_bf_coef[2], 2), "\n Corr. = ", round(corr, 2))) +
    geom_abline(slope = fm_bf_coef[2], intercept = uncentre) +
    # geom_abline(slope = fm_perp_coef, intercept = (fm_line_sd/cos(atan(fm_perp_coef)))+uncentre2, lty = 2) +
    # geom_abline(slope = fm_perp_coef, intercept = -(fm_line_sd/cos(atan(fm_perp_coef)))+uncentre2, lty = 2) +
    # geom_abline(slope = fm_perp_coef, intercept = (2*fm_line_sd/cos(atan(fm_perp_coef)))+uncentre2, lty = 3) +
    # geom_abline(slope = fm_perp_coef, intercept = (-2*fm_line_sd/cos(atan(fm_perp_coef)))+uncentre2, lty = 3) +
    # geom_abline(slope = fm_bf_coef[2], intercept = (-1*fm_line_sd2/cos(atan(fm_bf_coef[2])))+uncentre, lty = 2) +
    # geom_abline(slope = fm_bf_coef[2], intercept = -(-1*fm_line_sd2/cos(atan(fm_bf_coef[2])))+uncentre, lty = 2) +
    # geom_abline(slope = fm_bf_coef[2], intercept = (-2*fm_line_sd2/cos(atan(fm_bf_coef[2])))+uncentre, lty = 3) +
    # geom_abline(slope = fm_bf_coef[2], intercept = -(-2*fm_line_sd2/cos(atan(fm_bf_coef[2])))+uncentre, lty = 3) +
    theme_bw() +
    xlab("Male Jaccard of pairwise chapter") +
    ylab("Female Jaccard of pairwise chapter")
  ggsave("fm_chapter_pairwise_jaccard_scatter.png", plot = plot,
         units = "px", width = 1800, height = 1600)
  
  # create output df
  
  fm_jacc %<>%
    mutate(bf_sd_thresh = "within -1 SD", perp_devi_sd_thresh = "within +-1 SD")
  
  fm_jacc$bf_sd_thresh[fm_jacc$bf_dist < -1*fm_line_sd] <- "<-1SD"
  fm_jacc$bf_sd_thresh[fm_jacc$bf_dist < -2*fm_line_sd] <- "<-2SD"
  fm_jacc$perp_devi_sd_thresh[fm_line_pt2 > fm_line_sd2] <- ">1SD"
  fm_jacc$perp_devi_sd_thresh[fm_line_pt2 > 2*fm_line_sd2] <- ">2SD"
  fm_jacc$perp_devi_sd_thresh[fm_line_pt2 < -fm_line_sd2] <- "<-1SD"
  fm_jacc$perp_devi_sd_thresh[fm_line_pt2 < -2*fm_line_sd2] <- "<-2SD"
  
  # add in names
  
  abbrev <- c("A00-B99" = "Cert.infect.paras.dis.", 
              "C00-D49" = "Neoplasms", 
              "D50-D89" = "Dis.blood.org.immune.", 
              "E00-E89" = "Endocr.nutrit.metabol.", 
              "F01-F99" = "Mental.behav.dis.", 
              "G00-G99" = "Dis.nervous.sys.", 
              "H00-H59" = "Dis.eye.adnexa", 
              "H60-H95" = "Dis.ear.mast.proc.", 
              "I00-I99" = "Dis.circul.sys.", 
              "J00-J99" = "Dis.resp.sys.", 
              "K00-K95" = "Dis.digest.sys.",
              "L00-L99" = "Dis.skin.subcut.tiss.",
              "M00-M99" = "Dis.musculo.sys.connect.tiss.",
              "N00-N99" = "Dis.genitourinary.sys.")
  
  fm_jacc %<>%
    mutate(comp1_name = NA, comp2_name = NA)
  
  fm_jacc$comp1_name <- unname(abbrev)[match(fm_jacc$comp1, names(abbrev))]
  fm_jacc$comp2_name <- unname(abbrev)[match(fm_jacc$comp2, names(abbrev))]
  
  fm_jacc %<>%
    select(comp1, comp2, comp1_name, comp2_name, fem_jacc, male_jacc, bf_sd_thresh, perp_devi_sd_thresh)
  
  colnames(fm_jacc) <- c("Comparison 1", "Comparison 2", "Comparison 1 name", "Comparison 2 name", "Female Jaccard", "Male Jaccard",
                         "Best fit SD", "Perpendicular SD")
  
  write_tsv(fm_jacc, "fm_chapter_pairwise_jaccard_table.tsv")
  
  write_tsv(fm_jacc %>% filter(`Best fit SD` != "within -1 SD" | `Perpendicular SD` != "within +-1 SD"), "fm_chapter_pairwise_jaccard_table_filt.tsv")
  
}



