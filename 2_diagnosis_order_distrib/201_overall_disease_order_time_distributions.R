# overall icd10 block disease trajectory analysis

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(RColorBrewer)

# Variables

# ICD10 block dictionary
icd10_blocks_file = "icd10_blocks.tsv"

# functions

distrib_plot <- function(type, data_type, gender, max_val) {
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
      summarise(pct_25 = quantile(date, 0.25), pct_50 = median(date), pct_75 = quantile(date, 0.75)) %>%
      filter(!is.na(diag))
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
      block_cumfreq_qnts <- trajs_diags
      chapt_grps <- blocks$chapter[match(block_cumfreq_qnts$diag, blocks$block)]
      # block_cumfreq_qnts$name <- abbrev
      block_cumfreq_qnts$name <- abbrev[match(block_cumfreq_qnts$diag, names(abbrev))]
      block_cumfreq_qnts %<>% select(name, pct_25, pct_50, pct_75)
    } else if(data_type == "icd10_chapter") {
      # block_cumfreq_qnts$name <- str_to_title(block_cumfreq_qnts$name)
    }
    
    block_order <- block_cumfreq_qnts$name[order(block_cumfreq_qnts$pct_50)]
    block_cumfreq_qnts$name <- factor(block_cumfreq_qnts$name, levels = block_order)
    
  }
  
  # create colours
  
  if(data_type == "icd10_block") {
    n <- length(unique(chapt_grps))
  } else if(data_type == "icd10_chapter") {
    n <- nrow(block_cumfreq_qnts)
  }
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col_vector <- col_vector[5:(n+4)]
  
  if(data_type == "icd10_block") {
    block_cumfreq_qnts$chapter <- factor(chapt_grps, levels = sort(unique(chapt_grps)))
    names(col_vector) <- unique(chapt_grps)
  } else if(data_type == "icd10_chapter") {
    names(col_vector) <- unique(block_cumfreq_qnts$name)
  }
  
  #max_val <- max(block_cumfreq_qnts$pct_75)
  
  # split plot
  
  # split <- "lower"
  # 
  # upper <- c("A00-B99", "C00-D49", "E00-E89", "F01-F99", "G00-G99", "H00-H59")
  # lower <- c("H60-H95", "I00-I99", "J00-J99", "K00-K95", "L00-L99", "M00-M99", "N00-N99")
  # 
  # if(split == "upper") {
  #   split_use <- upper
  # } else if(split == "lower") {
  #   split_use <- lower
  # }
  # 
  # if(split == "lower") {
  #   tick_len <- 1
  # } else if(split == "upper") {
  #   tick_len <- 0
  # }
  # 
  # if(data_type == "icd10_block") {
  #   block_cumfreq_qnts %>%
  #     #filter(chapter %in% split_use) %>%
  #     ggplot(aes(x = pct_50, xmin = pct_25, xmax = pct_75, y = name, col = chapter)) +
  #     geom_pointrange(size = 0.75, shape = "diamond",) +
  #     scale_x_continuous(breaks = seq(ifelse(type == "diag_order", 1, 0), max_val, 1), limits = c(0, max_val)) +
  #     scale_colour_manual(values = col_vector, breaks = levels(chapt_grps)) +
  #     scale_y_discrete(guide = guide_axis(check.overlap = T)) +
  #     facet_grid(chapter ~ ., space="free", scales="free", switch="y") +
  #     theme_bw() +
  #     theme(axis.text.x = element_text(size = 11, face = "bold", colour = ifelse(split == "lower", "black", "white")),
  #           axis.text.y = element_text(size = 11, face = "bold"),
  #           axis.title = element_text(size = 12, face = "bold"),
  #           axis.ticks.length.x = unit(tick_len, "mm"),
  #           strip.text.y = element_blank(),
  #           legend.position = "none") +
  #     xlab(ifelse(type == "diag_order", "Diagnosis order", "Day interval/200 days")) +
  #     ylab("ICD10 block") +
  #     labs(color = "ICD10 Chapter")
  # } else if(data_type == "icd10_chapter") {
  #   block_cumfreq_qnts %>%
  #     ggplot(aes(x = pct_50, xmin = pct_25, xmax = pct_75, y = name)) +
  #     geom_pointrange(size = 0.75, shape = "diamond") +
  #     scale_x_continuous(breaks = seq(ifelse(type == "diag_order", 1, 0), max_val, 2), limits = c(0, max_val)) +
  #     scale_colour_manual(values = col_vector) +
  #     theme_bw() +
  #     theme(axis.text = element_text(size = 13)) +
  #     xlab(ifelse(type == "diag_order", "Diagnosis order", "Day interval/200 days")) +
  #     ylab("ICD10 chapter") +
  #     scale_y_discrete(labels = function(x) str_wrap(x, 50))
  # }
  
  # non split
  
  name_ordered <- block_cumfreq_qnts %>%
    arrange(pct_50) %>%
    pull(name)
  
  n_dis <- length(name_ordered)
  cutoffs <- round(sapply(1:3, function(x) (n_dis/4)*x))
  cutoff_names <- sapply(cutoffs, function(i) name_ordered[i])
  
  if(data_type == "icd10_block") {
    png(paste0(filename, "_distrib.png"), units = "px", width = 600, height = 900)
    print(block_cumfreq_qnts %>%
            ggplot(aes(x = pct_50, xmin = pct_25, xmax = pct_75, y = factor(name, levels = name_ordered))) +
            geom_point(shape = "diamond", size = 4) +
            geom_errorbarh(size = 0.75) +
            scale_x_continuous(breaks = seq(ifelse(type == "diag_order", 1, 0), max_val, ifelse(type == "diag_order", 1, 1)), limits = c(0, max_val)) +
            scale_colour_manual(values = col_vector, breaks = levels(chapt_grps)) +
            scale_y_discrete(guide = guide_axis(check.overlap = T)) +
            theme_bw() +
            theme(axis.text.x = element_text(size = 14, face = "bold"),
                  axis.text.y = element_text(size = 14, face = "bold"),
                  axis.title = element_text(size = 15, face = "bold"),
                  strip.text.y = element_blank(),
                  legend.position = "none") +
            xlab(ifelse(type == "diag_order", "Diagnosis order", "Time to diagnosis / Years")) +
            ylab("ICD10 block") +
            labs(color = "ICD10 Chapter") +
            geom_hline(yintercept = cutoff_names, col = "red"))
    dev.off()
  } else if(data_type == "icd10_chapter") {
    png(paste0(filename, "_distrib.png"), units = "px", width = 600, height = 900)
    print(block_cumfreq_qnts %>%
            ggplot(aes(x = pct_50, xmin = pct_25, xmax = pct_75, y = name)) +
            geom_point(shape = "diamond", size = 4) +
            geom_errorbarh(size = 0.75) +
            scale_x_continuous(breaks = seq(ifelse(type == "diag_order", 1, 0), max_val, 2), limits = c(0, max_val)) +
            scale_colour_manual(values = col_vector) +
            theme_bw() +
            theme(axis.text = element_text(size = 15)) +
            xlab(ifelse(type == "diag_order", "Diagnosis order", "Day interval/200 days")) +
            ylab("ICD10 chapter") +
            scale_y_discrete(labels = function(x) str_wrap(x, 50)))
    dev.off()
  }
  
  
  # chapter counts of cutoffs
  
  cutoff_1_chapt_cnts <- block_cumfreq_qnts %>%
    arrange(pct_25, pct_50, pct_75) %>%
    slice(1:cutoffs[1]) %>%
    count(chapter) %>%
    mutate(cutoff = 1)
  
  cutoff_2_chapt_cnts <- block_cumfreq_qnts %>%
    arrange(pct_25, pct_50, pct_75) %>%
    slice((cutoffs[1]+1):cutoffs[2]) %>%
    count(chapter) %>%
    mutate(cutoff = 2)
  
  cutoff_3_chapt_cnts <- block_cumfreq_qnts %>%
    arrange(pct_25, pct_50, pct_75) %>%
    slice((cutoffs[2]+1):cutoffs[3]) %>%
    count(chapter) %>%
    mutate(cutoff = 3)
  
  cutoff_4_chapt_cnts <- block_cumfreq_qnts %>%
    arrange(pct_25, pct_50, pct_75) %>%
    slice((cutoffs[3]+1):nrow(block_cumfreq_qnts)) %>%
    count(chapter) %>%
    mutate(cutoff = 4)
  
  # plot
  
  # cutoff_1_chapt_cnts %>%
  #   ggplot(aes(y = chapter, x = n)) +
  #   geom_bar(stat = "identity") +
  #   theme_bw() +
  #   ylab("ICD10 Chapter") +
  #   xlab("Count")
  #   
  # cutoff_2_chapt_cnts %>%
  #   ggplot(aes(y = chapter, x = n)) +
  #   geom_bar(stat = "identity") +
  #   theme_bw() +
  #   ylab("ICD10 Chapter") +
  #   xlab("Count")
  # 
  # cutoff_3_chapt_cnts %>%
  #   ggplot(aes(y = chapter, x = n)) +
  #   geom_bar(stat = "identity") +
  #   theme_bw() +
  #   ylab("ICD10 Chapter") +
  #   xlab("Count")
  # 
  # cutoff_4_chapt_cnts %>%
  #   ggplot(aes(y = chapter, x = n)) +
  #   geom_bar(stat = "identity") +
  #   theme_bw() +
  #   ylab("ICD10 Chapter") +
  #   xlab("Count")
  
  # combined plot
  
  cutoff_chapt_cnts <- rbind(cutoff_1_chapt_cnts,
                             cutoff_2_chapt_cnts,
                             cutoff_3_chapt_cnts,
                             cutoff_4_chapt_cnts)
  
  lab_names <- c("A00-B99" = "Certain infectious and parasitic diseases",
                 "C00-D49" = "Neoplasms",
                 "D50-D89" = "Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism",
                 "E00-E89" = "Endocrine, nutritional and metabolic diseases",
                 "F01-F99" = "Mental and behavioural disorders",
                 "G00-G99" = "Diseases of the nervous system",
                 "H00-H59" = "Diseases of the eye and adnexa",
                 "H60-H95" = "Diseases of the ear and mastoid process",
                 "I00-I99" = "Diseases of the circulatory system",
                 "J00-J99" = "Diseases of the respiratory system",
                 "K00-K95" = "Diseases of the digestive system",
                 "L00-L99" = "Diseases of the skin and subcutaneous tissue",
                 "M00-M99" = "Diseases of the musculoskeletal system and connective tissue",
                 "N00-N99" = "Diseases of the genitourinary system")
  
  cutoff_chapt_cnts$name <- lab_names[match(cutoff_chapt_cnts$chapter, names(lab_names))]
  cutoff_chapt_cnts$name <- factor(cutoff_chapt_cnts$name, levels = lab_names)
  
  png(paste0(filename, "_chapt_cnts.png"), units = "px", width = 600, height = 900)
  print(cutoff_chapt_cnts %>%
          ggplot(aes(y = name, x = n, fill = factor(cutoff, levels = 4:1))) +
          geom_bar(stat = "identity") +
          theme_bw() +
          ylab("ICD10 Chapter") +
          xlab("Number of Blocks") +
          labs(fill = "Diagnosis order\nquartile") +
          scale_x_continuous(breaks = seq(0, 10, 1)) +
          scale_y_discrete(labels = function(x) str_wrap(str_to_title(x), 40)) +
          scale_fill_brewer(type =  "div", palette = 5) +
          theme(legend.position = "none", 
                axis.title.x = element_text(size = 15, face = "bold"), 
                axis.title.y = element_text(size = 15, face = "bold"),
                axis.text.x = element_text(size = 14, face = "bold"),
                axis.text.y = element_text(size = 14, face = "bold")))
  dev.off()
}


add_nl <- function(x) {
  x_sub <- substring(x, 50)
  spc_loc <- str_locate(x_sub, " ")
  spc_loc <- spc_loc[, 2]
  x_sub2 <- paste0(substring(x, 1, 49), substring(x_sub, 1, spc_loc), "/n", substring(x_sub, spc_loc+1), collapse = "")
  return(x_sub2)
}

# make plots

# type = "diag_order" or "day_interval"

# data_type = "icd10_block" or "icd10_chapter"

# gender =  "m" or "f"

# max_val = max x value of block
# 36 for day interval block
# 11 for diag order for block
# 5 for diag order for chapter
# 16 for day interval chapter

# diagnosis order, male block

# male trajectory matrix
traj_file = "male_disease_block_1pct_traj.tsv"
# male relative time to diagnosis matrix
traj_dates_file = "male_diag_block_1pct_norm_dates.tsv"

distrib_plot("diag_order", "icd10_block", "m", 13)

# diagnosis order, female block

# female trajectory matrix
traj_file = "fem_disease_block_1pct_traj.tsv"

# female time to diagnosis matrix
traj_dates_file = "fem_diag_block_1pct_norm_dates.tsv"

distrib_plot("diag_order", "icd10_block", "f", 13)

# time to diagnosis, male block

# male trajectory matrix
traj_file = "male_disease_block_1pct_traj.tsv"
# male time to diagnosis matrix
traj_dates_file = "male_diag_block_1pct_norm_dates.tsv"

distrib_plot("day_interval", "icd10_block", "m", 14)

# time to diagnosis, female block

# female trajectory matrix
traj_file = "fem_disease_block_1pct_traj.tsv"
# female time to diagnosis matrix
traj_dates_file = "fem_diag_block_1pct_norm_dates.tsv"

distrib_plot("day_interval", "icd10_block", "f", 15)


