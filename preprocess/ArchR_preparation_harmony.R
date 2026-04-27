#!/usr/bin/env Rscript

# Run ArchR to generate metacell information for running ChromaFold

library(ArchR)
library(tidyr)
library(ggplot2)
library(harmony)

source("utils_archr.R")
args <- commandArgs(trailingOnly = TRUE)

#####################################
#      Step 0. Initial set-ups      #
#####################################

#1. Initial set-up
set.seed(1234)
archr_path <- args[1]
frag_dir <- args[2] #"/data1/lesliec/carolw/data_raw/scATAC-chan/Cancer_scATACseq_data"
genome_assembly <- args[3] #"hg38"
subsample <- as.logical(args[4])
cohort <- args[5]
cell_number <- NULL
if (subsample && length(args) >= 6) {
  cell_number <- as.integer(args[6])
}
setwd(archr_path)

archr_path <- "/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/archr_data_allcells"
frag_dir <- "/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/fragments_sorted"
genome_assembly <- "hg38"
subsample <- "FALSE"
proj <- loadArchRProject(
  path = archr_path,
  showLogo = FALSE
)
setwd(archr_path)

#2. Hyper-parameters
# addArchRGenome("mm10")  # nolint
addArchRGenome(genome_assembly)
addArchRThreads(threads = 1)
tile_mat_params <- list()
tile_mat_params$tileSize <- 500

#########################################
#      Step 1. Load fragmens files      #
#########################################

#1. Load fragments files
inputFiles <- getArchRInputFiles(
  frag_dir = frag_dir,
  cohort = cohort
)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  TileMatParams = tile_mat_params,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

print("Finished creating ArrowFiles")

#2. Create archr project
proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    # outputDirectory = paste0(archr_path, "/archr_out/"), #nolint
    outputDirectory = paste0(archr_path, "/"),
    copyArrows = FALSE
)

print("Finished creating Archr project")

# subsampling of cells
if (subsample) {
  cells.keep <- unlist(
    lapply(
      split(rownames(proj@cellColData), proj$Sample),
      function(cells) {
        if (length(cells) > cell_number) {
          sample(cells, cell_number)
        } else {
          cells
        }
      }
    )
  )

  sub_out <- file.path(dirname(archr_path), "archr_data_subsampled")
  proj <- subsetArchRProject(
    proj,
    cells = cells.keep,
    dropCells = TRUE,
    outputDirectory = sub_out
  )
  setwd(sub_out)
}

sub_out <- file.path(dirname(archr_path), "archr_data_subsampled")
output_path <- if (subsample) sub_out else archr_path
#3. Run lsi
proj <- addIterativeLSI(
  ArchRProj = proj, useMatrix = "TileMatrix",
  name = "IterativeLSI", iterations = 2,
  clusterParams = list(
    resolution = c(0.2), sampleCells = 10000, n.start = 10
  ),
  varFeatures = 25000, dimsToUse = 1:30
)

print("Finished running lsi")

#3.1. Harmony 
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)
print("Finished harmony batch correction")

#3.2. Clustering
proj <- addClusters(
    input = proj, reducedDims = "Harmony",
    method = "Seurat", name = "Clusters",
    resolution = 2, force = TRUE
)

print("Finished clustering")

#3.2. tsne
# proj <- addTSNE(
#     ArchRProj = proj, reducedDims = "IterativeLSI",
#     name = "TSNE", perplexity = 30
# )
proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "Harmony",
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5
)

####################################
#      Step 2. Save all files      #
####################################

#1. Save
final_bc <- rownames(proj@cellColData)
write.csv(final_bc, file = paste0(output_path,  "/archr_filtered_barcode.csv"))
lsi <- getReducedDims(proj, reducedDims = "Harmony")
write.csv(lsi, file = paste0(output_path, "/archr_filtered_lsi.csv"))
cat("After subsetting:", nrow(proj@cellColData), "\n")
cat("Final barcodes:", length(final_bc), "\n")
cat("LSI rows:", nrow(lsi), "\n")
saveArchRProject(ArchRProj = proj, outputDirectory = output_path, load = FALSE)

p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP",
  reducedDims = "Harmony"
)

outfile <- file.path(
  output_path,
  "UMAP_Sample_Harmony.png"
)

ggsave(
  filename = outfile,
  plot = p,
  width = 6,
  height = 5,
  dpi = 300
)

#2. save fragment
cell_names <- rownames(proj@cellColData)

frags <- getFragmentsFromProject(
  ArchRProj = proj,
  cellNames = cell_names
)

frags_unlisted <- unlist(frags)

frags_df <- data.frame(
  seqnames = seqnames(frags_unlisted),
  start = start(frags_unlisted),
  end = end(frags_unlisted),
  RG = frags_unlisted$RG,
  score = 1  # standard fragment file format uses 1
)

write.table(
  frags_df,
  file = gzfile(file.path(output_path, "merged_fragments.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

