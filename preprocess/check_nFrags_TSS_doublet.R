.libPaths('/data1/lesliec/carolw/my_envs/archr_2/lib/R/library')
library(ArchR)
library(tidyr)
library(dplyr)
library(ggplot2)

archr_path <- "/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA-chan/archr_data_subsampled"
out_path <- "/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA-chan/merged_UMAP"
setwd(out_path)

genome_assembly <- "hg38"
addArchRGenome(genome_assembly)
addArchRThreads(threads = 4)

proj <- loadArchRProject(
  path = archr_path,
  showLogo = FALSE
)

table <- getCellColData(proj, select = c("Sample", "nFrags", "TSSEnrichment")) %>%
  as.data.frame() %>%
  group_by(Sample) %>%
  summarise(
    n_cells = n(),
    total_fragments = sum(nFrags),
    median_frags = median(nFrags),
    median_TSS = median(TSSEnrichment)
  )

getCellColData(proj, select = c("Sample", "nFrags", "TSSEnrichment")) %>%
  as.data.frame() %>%
  group_by(Sample) %>%
  summarise(
    n_cells = n(),
    total_fragments = sum(nFrags),
    median_frags = median(nFrags),
    median_TSS = median(TSSEnrichment)
  )%>%
  write.csv("sample_fragment_TSS_summary.csv", row.names = FALSE)