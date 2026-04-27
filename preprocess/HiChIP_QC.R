.libPaths(c(
  "/data1/lesliec/carolw/my_envs/r442_seurat/lib/R/library",
  .libPaths()
))

library(dplyr)
library(readr)
library(tidyverse)

output_path <- "/data1/lesliec/carolw/projects/chromafold/preprocess"
setwd(output_path)

# ============================================================
# Get GI List 
# ============================================================
gbm_dir  <- "/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/hichip"
brca_dir <- "/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA-chan/hichip"
gbm_files <- list.files(gbm_dir, pattern = "*_norm.txt.gz$", 
                        recursive = TRUE, full.names = TRUE)
brca_files <- list.files(brca_dir, pattern = "*_norm.txt.gz$", 
                         recursive = TRUE, full.names = TRUE)
# GBM
gbm_list <- lapply(gbm_files, read_tsv)
names(gbm_list) <- basename(gbm_files) %>% gsub("_norm.txt.gz","",.)

# BRCA
brca_list <- lapply(brca_files, read_tsv)
names(brca_list) <- basename(brca_files) %>% gsub("_norm.txt.gz","",.)
gi_list <- c(gbm_list, brca_list)

brca_map <- c(
  "BRCA-08499A64" = 'LumA2',
  "BRCA-14AD76EE" = 'LumB2',
  "BRCA-5C54B79C" = 'LumA1',
  "BRCA-7C6A3AE4" = 'Basal1',
  "BRCA-8D1E6006" = 'LumB1',
  "BRCA-94AF19F0" = 'Her2',
  "BRCA-C147AAD5" = 'Basal2',
  "BRCA-C9C8D426" = 'LumB3',
  "BRCA-CB96A542" = 'Basal3',
  "BRCA-DD69EDE9" = 'Basal4'
)

gbm_map <- c(
  "GBMx-09C0DCE7" = 5,
  "GBMx-6BEE2CB6" = 4,
  "GBMx-9976F952" = 6,
  "GBMx-A90B18B6" = 7
)


sample_ids <- sub("^(GBMx|BRCA)-([0-9A-F]+).*", "\\1-\\2", names(gi_list))

new_names <- sapply(sample_ids, function(id){
  if(grepl("^BRCA", id)){
    paste0("BRCA_", brca_map[id])
  } else if(grepl("^GBM", id) | grepl("^GBMx", id)){
    paste0("GBM_", gbm_map[id])
  } else {
    id
  }
})

names(gi_list) <- new_names

# ============================================================
# Sig interaction distirbution per sample
# ============================================================
sig_counts <- sapply(gi_list, function(gi){
  nrow(gi)  # HiC-DC+ output usually already filtered for significance
})

dist_list <- lapply(gi_list, function(gi){
  # Keep only cis interactions
  gi <- gi[gi$chrI == gi$chrJ, ]
  gi$dist <- abs((gi$startI + gi$endI)/2 - (gi$startJ + gi$endJ)/2)
  gi$dist
})

# Convert to long format
dist_df <- do.call(rbind, lapply(names(dist_list), function(samp){
  data.frame(sample = samp, distance = dist_list[[samp]])
}))

all_distances <- unlist(dist_list)
universal_quantiles <- quantile(all_distances, probs = c(0.25, 0.5, 0.75))
universal_quantiles

p<-ggplot(dist_df, aes(x = distance/1000, color = sample)) +
  geom_density() +
  scale_x_log10() +
  labs(x = "Interaction distance (kb, log10 scale)", y = "Density") +
  theme_classic()

uq_kb <- universal_quantiles / 1000

ggsave(
  "Interacion_dist_distribution_per_sample.pdf",
  plot = print(p + 
                 geom_vline(xintercept = uq_kb, linetype = "dashed", color = "black") +
                 annotate("text", x = uq_kb, y = 0, label = c("25%", "50%", "75%"),
                          angle = 90, vjust = -0.5, hjust = 0, size = 3)),
  width = 10,
  height = 5
)


# ============================================================
# Define distance threshold
# ============================================================

short_range <- 170000   # <170kb
mid_range   <- 450000   # <450kb
long_range  <- 980000   # <980 kb
#ultra_long  <- 2e6   # <2Mb

range_counts <- lapply(dist_list, function(d){
  data.frame(
    short   = sum(d < short_range),
    mid     = sum(d >= short_range & d < mid_range),
    long    = sum(d >= mid_range & d < long_range),
    ultra   = sum(d >= long_range)
  )
})

range_df <- do.call(rbind, lapply(names(range_counts), function(samp){
  df <- range_counts[[samp]]
  df$sample <- samp
  reshape2::melt(df, id.vars = "sample", variable.name = "range", value.name = "count")
}))

p<- ggplot(range_df, aes(x = sample, y = count, fill = range)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(y = "Number of interactions", fill = "Distance range")

ggsave(
  "Interacion_dist_TH_per_sample.pdf",
  plot = print(p),
  width = 10,
  height = 5
)

# ============================================================
# Valid Pairs
# ============================================================
df <- read_tsv("/data1/lesliec/carolw/data_raw/TCGA_HiChIP_allValidPairs/used/valid_pairs_summary.tsv")

df <- df %>%
  mutate(sample_id = str_extract(sample, "^[A-Za-z]+-[A-Za-z0-9]+"))

df <- df %>%
  mutate(
    sample_label = case_when(
      sample_id %in% names(brca_map) ~ paste0("BRCA_", brca_map[sample_id]),
      sample_id %in% names(gbm_map) ~ paste0("GBM_", gbm_map[sample_id]),
      TRUE ~ sample_id
    )
  )

df$sample_label <- factor(df$sample_label,
                          levels = sort(unique(df$sample_label)))
df <- df %>%
  mutate(cancer = str_extract(sample_id, "^[A-Za-z]+"))

ggplot(df, aes(x = sample_label, y = valid_pairs/1e6,
               fill = cancer)) +
  geom_bar(stat = "identity") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Sample",
    y = "Valid pairs (millions)",
    title = "HiChIP sequencing depth per sample"
  ) 

ggplot(df, aes(x = cancer, y = valid_pairs/1e6, fill=cancer)) +
  geom_boxplot() +
  geom_jitter(width=0.15, alpha=0.7) +
  labs(
    x="Caner Type",
    y="Valid pairs (millions)",
    title="HiChIP depth by cancer type"
  ) +
  theme_classic()