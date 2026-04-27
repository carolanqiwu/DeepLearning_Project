library(ggplot2)
.libPaths('/data1/lesliec/carolw/my_envs/archr_2/lib/R/library')

library(ArchR)
#library(tidyr)
library(dplyr)
library(patchwork)
library(cowplot)
library(BSgenome.Hsapiens.UCSC.hg38)

genome_assembly <- "hg38"
addArchRGenome(genome_assembly)
addArchRThreads(threads = 4)
output_path <- "/data1/lesliec/carolw/projects/chromafold/preprocess"
setwd(output_path)
# ============================================================
# Get all ArchR proj
# ============================================================
brca_dir <- "/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA-chan/1_3k/archr_data"
gbm_dir  <- "/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/1_3k/archr_data"
brca_projects <- list.dirs(brca_dir, recursive = FALSE)
gbm_projects  <- list.dirs(gbm_dir, recursive = FALSE)

all_projects <- c(brca_projects, gbm_projects)
proj_list <- lapply(all_projects, loadArchRProject)
names(proj_list) <- basename(all_projects)

brca_base <- basename(brca_projects)
gbm_base  <- basename(gbm_projects)
brca_map <- c(
  "BRCA_08499A64_3FD8_4E62_AF08_3C66AF93CAE7_X009_S01_B1_T1" = 'LumA2',
  "BRCA_14AD76EE_12F9_40B3_8DCD_4A256E02CF8D_X003_S02_B1_T1" = 'LumB2',
  "BRCA_5C54B79C_DA02_4B22_9FC2_3D61BFFC5559_X011_S02_B1_T1" = 'LumA1',
  "BRCA_7C6A3AE4_E2EA_42B3_B3F1_81C19E6F2170_X005_S02_B1_T1" = 'Basal1',
  "BRCA_8D1E6006_85CB_484A_8B5C_30766D90137B_X003_S03_B1_T2" = 'LumB1',
  "BRCA_94AF19F0_1F2A_41EC_8CB6_96C76227811F_X013_S01_B1_T1" = 'Her2',
  "BRCA_C147AAD5_A8F1_41D5_8709_21820BE50902_X008_S02_B1_T1" = 'Basal2',
  "BRCA_C9C8D426_A3FD_4455_89A9_768BC01D66A9_X009_S02_B1_T1" = 'LumB3',
  "BRCA_CB96A542_7AC1_4FEC_A5D2_458D8EEDC6C4_X013_S06_B1_T1" = 'Basal3',
  "BRCA_DD69EDE9_142D_46E2_AA06_58D07D3230FB_X014_S08_B1_T1" = 'Basal4'
)

gbm_map <- c(
  "GBMx_09C0DCE7_D669_4D28_980D_BF71179116A4_X005_S04_B1_T1" = 5,
  "GBMx_6BEE2CB6_9AFD_42A6_9C26_9C4428FBABFA_X004_S04_B1_T1" = 4,
  "GBMx_9976F952_23A5_431A_A431_01E544324A26_X010_S05_B1_T1" = 6,
  "GBMx_A90B18B6_6056_46D1_BF4E_70710236B8DD_X006_S04_B1_T1" = 7
)
brca_names <- paste0("BRCA_", brca_map[brca_base])
gbm_names  <- paste0("GBM_", gbm_map[gbm_base])
names(proj_list) <- c(brca_names, gbm_names)

sample_info <- data.frame(
  sample = names(proj_list),
  cancer = ifelse(grepl("BRCA", names(proj_list)), "BRCA", "GBM")
)

# ============================================================
# Extract QC
# ============================================================
qc_list <- lapply(names(proj_list), function(samp){
  
  proj <- proj_list[[samp]]
  df <- as.data.frame(proj@cellColData)
  
  df <- df %>%
    dplyr::select(TSSEnrichment, nFrags) %>%
    mutate(
      sample = samp,
      cancer = ifelse(grepl("BRCA", samp), "BRCA", "GBM")
    )
  
  return(df)
})

qc_df <- dplyr::bind_rows(qc_list)
qc_df$log10_nFrags <- log10(qc_df$nFrags)

# ============================================================
# Plot TSS & nFrags
# ============================================================
cols <- c(
  BRCA = "#E64B35",
  GBM = "#4DBBD5"
)

qc_df$sample <- factor(
  qc_df$sample,
  levels = qc_df %>%
    dplyr::distinct(sample, cancer) %>%
    arrange(cancer, sample) %>%
    dplyr::pull(sample)
)

p1 <- ggplot(qc_df, aes(x = sample, y = TSSEnrichment, fill = cancer)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(y = "TSS Enrichment", x = "")

p2 <- ggplot(qc_df, aes(x = sample, y = log10_nFrags, fill = cancer)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(y = "log10(nFrags)", x = "")

combined_plot <- plot_grid(p1, p2, ncol = 1)

ggsave(
  "qc_tss_fragments.pdf",
  plot = print(combined_plot),
  width = 10,
  height = 5
)