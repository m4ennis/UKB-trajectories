# pairwise Jaccard coincidence heatmap

library(tidyverse)
library(ComplexHeatmap)
library(cluster)
library(magrittr)

# Variables

# female Jaccard coincidence matrix
fem_mm_jacc_file = "ukb_fem_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_table.tsv"

# male Jaccard coincidence matrix
male_mm_jacc_file = "ukb_male_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_table.tsv"

# female prevalence filtered multimorbidity matrix
fem_mm_file = "ukb_fem_multimorbidity_ICD10_physio_block_1pct_table.tsv"

# male prevalence filtered multimorbidity matrix
male_mm_file = "ukb_male_multimorbidity_ICD10_physio_block_1pct_table.tsv"

# ICD10 blocks dictionary
icd10_blocks_file = "icd10_blocks.tsv"

# block or chapter
icd10_type = "block"

# group blocks into chapters?

blocks_grouped <- T


# Read in multimorbidity pairwise matrix  

fem_jacc <- read_tsv(fem_mm_jacc_file)
fem_jacc <- as.matrix(fem_jacc)

male_jacc <- read_tsv(male_mm_jacc_file)
male_jacc <- as.matrix(male_jacc)

# read in multimorbidity matrix

fem_mm <- read_tsv(fem_mm_file)
fem_mm <- as.matrix(fem_mm)

male_mm <- read_tsv(male_mm_file)
male_mm <- as.matrix(male_mm)

# get sizes

fem_sizes <- signif(colSums(fem_mm)/nrow(fem_mm), 2)*100
male_sizes <- signif(colSums(male_mm)/nrow(male_mm), 2)*100

# read in icd10 blocks table

icd10_blocks <- read_tsv(icd10_blocks_file)

if(icd10_type == "block") {

        # convert block codes to names
        
        colnames(fem_jacc) <- icd10_blocks$name[match(colnames(fem_jacc), icd10_blocks$block)]
        
        colnames(male_jacc) <- icd10_blocks$name[match(colnames(male_jacc), icd10_blocks$block)]
        
        # set row names
        
        rownames(fem_jacc) <- colnames(fem_jacc)
        
        rownames(male_jacc) <- colnames(male_jacc)
        
        # add NAs
        
        diag(fem_jacc) <- NA
        diag(male_jacc) <- NA
        
        # colours
        
        col_fun <- circlize::colorRamp2(c(0, 0.1, 0.2), c("lightblue", "yellow1", "red"))
        
        # female block version
        
        if(blocks_grouped) {
          
          fem_row_anno <- rowAnnotation(row_bar = anno_barplot(fem_sizes, axis_param = list(side = "top", gp = gpar(fontsize = 10))), show_annotation_name = F)
          fem_col_anno <- columnAnnotation(col_bar = anno_barplot(fem_sizes, axis_param = list(gp = gpar(fontsize = 10))), show_annotation_name = F)
          
          png("ukb_fem_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_heatmap.png", width = 800, height = 700)
          
          print(Heatmap(fem_jacc,
                        bottom_annotation = fem_col_anno,
                        left_annotation = fem_row_anno,
                        na_col = "white",
                        show_row_names = F,
                        show_column_names = F,
                        row_gap = unit(0.1, "mm"),
                        col = col_fun,
                        column_split = icd10_blocks$chapter[match(colnames(fem_jacc), icd10_blocks$name)],
                        row_split = icd10_blocks$chapter[match(rownames(fem_jacc), icd10_blocks$name)],
                        row_title_gp = gpar(fontsize = 18, fill = "white"),
                        column_title_gp = gpar(fontsize = 18, fill = "white"),
                        row_title_rot = 0,
                        column_title_rot = 90,
                        column_title_side = "bottom",
                        border = "white",
                        column_gap = unit(0.1, "mm"),
                        column_dend_height = unit(50, "mm"),
                        cluster_row_slices = T,
                        cluster_rows = T,
                        cluster_columns = T,
                        cluster_column_slices = T,
                        show_row_dend = F,
                        heatmap_legend_param = list(
                          title = "  Jaccard\ncoefficient\n",
                          title_gp = gpar(fontsize = 22, fontface = "bold"),
                          grid_height = unit(12, "mm"),
                          grid_width = unit(16, "mm"),
                          labels = c("0", "0.05", "0.1", "0.15", ">0.2"),
                          border = T,
                          title_position = "topleft",
                          labels_gp = gpar(fontsize = 18, fontface = "bold")
                        )))
          
        dev.off()
          
          
          # male block version
          
          male_row_anno <- rowAnnotation(row_bar = anno_barplot(male_sizes, axis_param = list(side = "top", gp = gpar(fontsize = 10))), show_annotation_name = F)
          male_col_anno <- columnAnnotation(col_bar = anno_barplot(male_sizes, axis_param = list(gp = gpar(fontsize = 10))), show_annotation_name = F)
          
          png("ukb_male_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_heatmap.png", width = 800, height = 700)
          
          print(Heatmap(male_jacc,
                        bottom_annotation = male_col_anno,
                        left_annotation = male_row_anno,
                        na_col = "white",
                        show_row_names = F,
                        show_column_names = F,
                        row_gap = unit(0.1, "mm"),
                        col = col_fun,
                        column_split = icd10_blocks$chapter[match(colnames(male_jacc), icd10_blocks$name)],
                        row_split = icd10_blocks$chapter[match(rownames(male_jacc), icd10_blocks$name)],
                        row_title_gp = gpar(fontsize = 18, fill = "white"),
                        column_title_gp = gpar(fontsize = 18, fill = "white"),
                        row_title_rot = 0,
                        column_title_rot = 90,
                        column_title_side = "bottom",
                        border = "white",
                        column_gap = unit(0.1, "mm"),
                        column_dend_height = unit(50, "mm"),
                        cluster_row_slices = T,
                        cluster_rows = T,
                        cluster_columns = T,
                        cluster_column_slices = T,
                        show_row_dend = F,
                        heatmap_legend_param = list(
                          title = "  Jaccard\ncoefficient\n",
                          title_gp = gpar(fontsize = 22, fontface = "bold"),
                          grid_height = unit(12, "mm"),
                          grid_width = unit(16, "mm"),
                          border = T,
                          labels = c("0", "0.05", "0.1", "0.15", ">0.2"),
                          title_position = "topleft",
                          labels_gp = gpar(fontsize = 18, fontface = "bold")
                        )))
          
          dev.off()
          
          fem_jacc_df <- fem_jacc %>%
            as.data.frame() %>%
            mutate(comp1 = rownames(fem_jacc)) %>%
            pivot_longer(-comp1, names_to = "comp2", values_to = "jacc") %>%
            mutate(sex = "female")
          
          male_jacc_df <- male_jacc %>%
            as.data.frame() %>%
            mutate(comp1 = rownames(male_jacc)) %>%
            pivot_longer(-comp1, names_to = "comp2", values_to = "jacc") %>%
            mutate(sex = "male")
          
          fm_jacc_df <- rbind(fem_jacc_df,
                              male_jacc_df)
  
          write_tsv(fm_jacc_df, file = "ukb_fm_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_table.tsv")
          
        } else {
          
          fem_model <- hclust(dist(fem_jacc))
          male_model <- hclust(dist(male_jacc))
          
          fem_sil_k <- sapply(2:15, function(k) {
            cut_k <- cutree(fem_model, k = k)
            sil_k <- summary(silhouette(cut_k, dist = dist(fem_jacc)))$avg.width
            return(sil_k)
          })
          
          fem_k <- (3:15)[which.max(fem_sil_k[-1])]
          fem_clus <- cutree(fem_model, k = fem_k)
          fem_clus <- data.frame(block = names(fem_clus), sex = "female", clus = as.vector(fem_clus))
          
          male_sil_k <- sapply(2:15, function(k) {
            cut_k <- cutree(male_model, k = k)
            sil_k <- summary(silhouette(cut_k, dist = dist(male_jacc)))$avg.width
            return(sil_k)
          })
          
          male_k <- (3:15)[which.max(male_sil_k[-1])]
          male_clus <- cutree(male_model, k = male_k)
          male_clus <- data.frame(block = names(male_clus), sex = "male", clus = as.vector(male_clus))
          
          fm_clus <- rbind(fem_clus,
                           male_clus)
          
          fm_clus %<>% pivot_wider(names_from = "sex", values_from = "clus", names_prefix = "clus_")
          
          fem_jacc_df <- fem_jacc %>%
            as.data.frame() %>%
            mutate(comp1 = rownames(fem_jacc)) %>%
            pivot_longer(-comp1, names_to = "comp2", values_to = "jacc") %>%
            mutate(sex = "female")
          
          male_jacc_df <- male_jacc %>%
            as.data.frame() %>%
            mutate(comp1 = rownames(male_jacc)) %>%
            pivot_longer(-comp1, names_to = "comp2", values_to = "jacc") %>%
            mutate(sex = "male")
          
          fm_jacc_df <- rbind(fem_jacc_df,
                              male_jacc_df)
          
          write_tsv(fm_clus, file = "ukb_fm_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_cluster_table.tsv")
          
          write_tsv(fm_jacc_df, file = "ukb_fm_mm_physio_1pct_block_morb_gt_1_pairwise_jaccard_table.tsv")
          
          png("ukb_fem_mm_physio_1pct_block_morb_gt_1_pairwise_ungrouped_jaccard_heatmap.png", width = 800, height = 700)
          
          fem_row_anno <- rowAnnotation(row_bar = anno_barplot(fem_sizes, axis_param = list(side = "top")), show_annotation_name = F)
          fem_col_anno <- columnAnnotation(col_bar = anno_barplot(fem_sizes), show_annotation_name = F)
          
          print(Heatmap(fem_jacc,
                        bottom_annotation = fem_col_anno,
                        left_annotation = fem_row_anno,
                        na_col = "white",
                        show_row_names = F,
                        show_column_names = F,
                        row_gap = unit(0.1, "mm"),
                        col = col_fun,
                        column_split = fem_clus$clus,
                        show_row_dend = F,
                        show_column_dend = F,
                        row_split = fem_clus$clus,
                        row_title_gp = gpar(fontsize = 18, fill = "white"),
                        column_title_gp = gpar(fontsize = 18, fill = "white"),
                        row_title_rot = 0,
                        column_title_rot = 90,
                        column_title_side = "bottom",
                        border = "black",
                        column_gap = unit(0.1, "mm"),
                        column_dend_height = unit(50, "mm"),
                        cluster_row_slices = F,
                        cluster_rows = T,
                        cluster_columns = T,
                        cluster_column_slices = F,
                        heatmap_legend_param = list(
                          title = "  Jaccard\ncoefficient\n",
                          title_gp = gpar(fontsize = 22, fontface = "bold"),
                          grid_height = unit(12, "mm"),
                          grid_width = unit(16, "mm"),
                          border = T,
                          title_position = "topleft",
                          labels_gp = gpar(fontsize = 18, fontface = "bold")
                        )))
          
          dev.off()
          
          
          # male block version
          
          male_row_anno <- rowAnnotation(row_bar = anno_barplot(male_sizes, axis_param = list(side = "top")), show_annotation_name = F)
          male_col_anno <- columnAnnotation(col_bar = anno_barplot(male_sizes), show_annotation_name = F)
          
          png("ukb_male_mm_physio_1pct_block_morb_gt_1_pairwise_ungrouped_jaccard_heatmap.png", width = 800, height = 700)
          
          print(Heatmap(male_jacc,
                        bottom_annotation = male_col_anno,
                        left_annotation = male_row_anno,
                        na_col = "white",
                        show_row_names = F,
                        show_column_names = F,
                        row_gap = unit(0.1, "mm"),
                        col = col_fun,
                        column_split = male_clus$clus,
                        row_split = male_clus$clus,
                        show_row_dend = F,
                        show_column_dend = F,
                        row_dend_reorder = F,
                        column_dend_reorder = F,
                        row_title_gp = gpar(fontsize = 18, fill = "white"),
                        column_title_gp = gpar(fontsize = 18, fill = "white"),
                        row_title_rot = 0,
                        column_title_rot = 90,
                        column_title_side = "bottom",
                        border = "black",
                        column_gap = unit(0.1, "mm"),
                        column_dend_height = unit(50, "mm"),
                        cluster_row_slices = F,
                        cluster_rows = T,
                        cluster_columns = T,
                        cluster_column_slices = F,
                        heatmap_legend_param = list(
                          title = "  Jaccard\ncoefficient\n",
                          title_gp = gpar(fontsize = 22, fontface = "bold"),
                          grid_height = unit(12, "mm"),
                          grid_width = unit(16, "mm"),
                          border = T,
                          title_position = "topleft",
                          labels_gp = gpar(fontsize = 18, fontface = "bold")
                        )))
          
          dev.off()
        }
        
        
} else if(icd10_type == "chapter") {
        
        # convert chapter codes to names
        
        #colnames(fem_jacc) <- icd10_blocks$chapter_name[match(colnames(fem_jacc), icd10_blocks$chapter)]
        
        #colnames(male_jacc) <- icd10_blocks$chapter_name[match(colnames(male_jacc), icd10_blocks$chapter)]
        
        # set row names
        
        rownames(fem_jacc) <- colnames(fem_jacc)
        
        rownames(male_jacc) <- colnames(male_jacc)
        
        # add NAs
        
        diag(fem_jacc) <- NA
        diag(male_jacc) <- NA
        
        # colours
        
        col_fun <- circlize::colorRamp2(c(0, 0.25, 0.5), c("lightblue", "yellow1", "red"))
        
        # female chapter version
        
        fem_model <- hclust(dist(fem_jacc))
        male_model <- hclust(dist(male_jacc))

        fem_sil_k <- sapply(2:10, function(k) {
          cut_k <- cutree(fem_model, k = k)
          sil_k <- summary(silhouette(cut_k, dist = dist(fem_jacc)))$avg.width
          return(sil_k)
        })
        
        fem_k <- (2:10)[which.max(fem_sil_k)]
        fem_clus <- cutree(fem_model, k = fem_k)
        fem_clus <- data.frame(block = names(fem_clus), sex = "female", clus = as.vector(fem_clus))
        
        male_sil_k <- sapply(2:10, function(k) {
          cut_k <- cutree(male_model, k = k)
          sil_k <- summary(silhouette(cut_k, dist = dist(male_jacc)))$avg.width
          return(sil_k)
        })
        
        male_k <- (2:10)[which.max(male_sil_k)]
        male_clus <- cutree(male_model, k = male_k)
        male_clus <- data.frame(block = names(male_clus), sex = "male", clus = as.vector(male_clus))
        
        fm_clus <- rbind(fem_clus,
                         male_clus)
        
        fm_clus %<>% pivot_wider(names_from = "sex", values_from = "clus", names_prefix = "clus_")
        
        fem_jacc_df <- fem_jacc %>%
          as.data.frame() %>%
          mutate(comp1 = rownames(fem_jacc)) %>%
          pivot_longer(-comp1, names_to = "comp2", values_to = "jacc") %>%
          mutate(sex = "female")
        
        male_jacc_df <- male_jacc %>%
          as.data.frame() %>%
          mutate(comp1 = rownames(male_jacc)) %>%
          pivot_longer(-comp1, names_to = "comp2", values_to = "jacc") %>%
          mutate(sex = "male")
        
        fm_jacc_df <- rbind(fem_jacc_df,
                            male_jacc_df)
        
        write_tsv(fm_clus, file = "ukb_fm_mm_physio_chapter_morb_gt_1_pairwise_jaccard_cluster_table.tsv")
        
        write_tsv(fm_jacc_df, file = "ukb_fm_mm_physio_chapter_morb_gt_1_pairwise_jaccard_table.tsv")
        
        fem_row_anno <- rowAnnotation(row_bar = anno_barplot(fem_sizes, axis_param = list(side = "top", gp = gpar(fontsize = 10), at = c(0, 20, 40))), show_annotation_name = F)
        fem_col_anno <- columnAnnotation(col_bar = anno_barplot(fem_sizes, axis_param = list(gp = gpar(fontsize = 10), at = c(0, 20, 40))), show_annotation_name = F)
        
        png("ukb_fem_mm_physio_chapter_morb_gt_1_pairwise_jaccard_heatmap.png", width = 800, height = 700)
        
        print(Heatmap(fem_jacc,
                      bottom_annotation = fem_col_anno,
                      left_annotation = fem_row_anno,
                      na_col = "white",
                      show_row_names = T,
                      show_column_names = T,
                      # column_split = fem_k,
                      # row_split = fem_k,
                      row_names_side = "left",
                      row_names_gp = gpar(fill = "white", fontsize = 18),
                      column_names_gp = gpar(fill = "white", fontsize = 18),
                      row_gap = unit(0.1, "mm"),
                      col = col_fun,
                      row_title_rot = 0,
                      column_title_rot = 90,
                      column_title_side = "bottom",
                      border = "black",
                      column_gap = unit(0.1, "mm"),
                      column_dend_height = unit(50, "mm"),
                      cluster_rows = T,
                      cluster_columns = T,
                      show_row_dend = F,
                      heatmap_legend_param = list(
                        title = "  Jaccard\ncoefficient\n",
                        title_gp = gpar(fontsize = 22, fontface = "bold"),
                        grid_height = unit(12, "mm"),
                        grid_width = unit(16, "mm"),
                        border = T,
                        title_position = "topleft",
                        labels_gp = gpar(fontsize = 18, fontface = "bold")
                )))
        
        dev.off()
        
        
        # male chapter version
        
        male_row_anno <- rowAnnotation(row_bar = anno_barplot(male_sizes, axis_param = list(side = "top", gp = gpar(fontsize = 10), at = c(0, 20, 40))), show_annotation_name = F)
        male_col_anno <- columnAnnotation(col_bar = anno_barplot(male_sizes, axis_param = list(gp = gpar(fontsize = 10), at = c(0, 20, 40))), show_annotation_name = F)
        
        png("ukb_male_mm_physio_chapter_morb_gt_1_pairwise_jaccard_heatmap.png", width = 800, height = 700)
        
        print(Heatmap(male_jacc,
                      bottom_annotation = male_col_anno,
                      left_annotation = male_row_anno,
                      na_col = "white",
                      show_row_names = T,
                      show_column_names = T,
                      row_names_side = "left",
                      # column_split = male_k,
                      # row_split = male_k,
                      row_names_gp = gpar(fill = "white", fontsize = 18),
                      column_names_gp = gpar(fill = "white", fontsize = 18),
                      row_gap = unit(0.1, "mm"),
                      col = col_fun,
                      row_title_rot = 0,
                      column_title_rot = 90,
                      column_title_side = "bottom",
                      border = "black",
                      column_gap = unit(0.1, "mm"),
                      column_dend_height = unit(50, "mm"),
                      cluster_rows = T,
                      cluster_columns = T,
                      show_row_dend = F,
                      heatmap_legend_param = list(
                        title = "  Jaccard\ncoefficient\n",
                        title_gp = gpar(fontsize = 22, fontface = "bold"),
                        grid_height = unit(12, "mm"),
                        grid_width = unit(16, "mm"),
                        border = T,
                        title_position = "topleft",
                        labels_gp = gpar(fontsize = 18, fontface = "bold")
                )))
        
        dev.off()
        
        
}

