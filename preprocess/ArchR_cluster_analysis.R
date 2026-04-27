#!/usr/bin/env Rscript

.libPaths('/data1/lesliec/carolw/my_envs/archr_2/lib/R/library')
library(ArchR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readxl)
library(BSgenome.Hsapiens.UCSC.hg38)

args <- commandArgs(trailing = TRUE)

#####################################
#      Step 0. Initial set-ups      #
#####################################

#1. Initial set-up
archr_path <- "/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/archr_data_allcells"
out_path <- "/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/archr_data_allcells"
genome_assembly <- "hg38"

#archr_path <- args[1]
#out_path <- args[2]
#genome_assembly <- args[3]

if(!dir.exists(out_path)){
  dir.create(out_path, recursive = TRUE)
}
setwd(out_path)

#2. Hyper-parameters
addArchRGenome(genome_assembly)
addArchRThreads(threads = 4)

#########################################
#      Step 1. Load ArchR Proj      #
#########################################

proj <- loadArchRProject(
  path = archr_path,
  showLogo = FALSE
)

# sanity check
table(proj$Clusters)
clusters <- data.frame(table(proj$Clusters))
write.csv(clusters, "clusters_count.csv", row.names = FALSE)
#########################################
#      Step 2. Add GeneScoreMatrix      #
#########################################

if (!"GeneScoreMatrix" %in% getAvailableMatrices(proj)) {
  proj <- addGeneScoreMatrix(
    proj,
    force = TRUE
  )
}

#########################################
#      Step 3. Marker Discovery      #
#########################################
# ArchR github issue #124
# a named list of marker results, one per cluster
markersGS <- lapply(unique(proj$Clusters), function(cl) {
  getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    useGroups = cl
  )
})
names(markersGS) <- unique(proj$Clusters)

# Extract significant markers
markerList <- lapply(markersGS, function(x) {
  getMarkers(x, cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
})

saveRDS(markersGS, file = file.path(out_path, "GeneScore_markers.rds"))
markersGS<-readRDS("GeneScore_markers.rds")
#########################################
#      Step 4. Gene Modules      #
#########################################
# Example marker modules 
df <- read.csv("/data1/lesliec/carolw/data_raw/chromafold_metadata/canonical_marker.csv", stringsAsFactors = FALSE)
#df <- read.csv("/data1/lesliec/carolw/data_raw/chromafold_metadata/BRCA_subtype_markers.csv", stringsAsFactors = FALSE)
#df <- read.csv("/data1/lesliec/carolw/data_raw/chromafold_metadata/genes_by_GBMsubtype.csv")

gene_modules <- lapply(df, function(x) {
  x <- unique(na.omit(x))
  x[x != ""]
})

# Keep only genes present in matrix
genes <- getFeatures(proj, useMatrix = "GeneScoreMatrix")

gene_modules <- lapply(gene_modules, function(x) {
  x[x %in% genes]
})


proj <- addModuleScore(
  proj,
  useMatrix = "GeneScoreMatrix",
  name = "Subtype_Module",
  features = gene_modules
)

proj <- addImputeWeights(proj, reducedDims="Harmony")


#########################################
#      Step 6. Visualization      #
#########################################
# ---- UMAP colored by clusters ----
proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
# Save as PNG
png(file.path(out_path, "UMAP_clusters.png"), width = 1200, height = 1000, res = 150) 

plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "Clusters"
)

# Close the device
dev.off()

# ---- UMAP colored by GeneScore ----

my_pal <- colorRampPalette(c("blue", "white", "red"))(100)
genes_of_interest <- c(
  "TFAP2C","TFAP2A"
)

#"ESR1" , "PGR", "ERBB2", # malignant BRCA
#"SOX2", "EGFR",       # malignant GBM
for (g in genes_of_interest) {
  png(
    filename = file.path(out_path, paste0("UMAP_GeneScore_", g, ".png")),
    width = 2000,
    height = 1600,
    res = 300
  )
  
  p<-plotEmbedding(
    ArchRProj = proj,
    colorBy = "GeneScoreMatrix",
    name = g,
    imputeWeights = getImputeWeights(proj),
    pal = my_pal
  )
  
  print(p)
  
  dev.off()
}

# ---- UMAP colored by module ----
#names(proj@cellColData)

png(file.path(out_path, "UMAP_modules_NL.png"), width = 1200, height = 1000, res = 150) 

p<-plotEmbedding(proj,
                    embedding = "UMAP",
                    colorBy = "cellColData",
                    name="GBM_subtype_Module.NL",
                    imputeWeights = getImputeWeights(proj),
              pal = my_pal)

print(p)
dev.off()

png(file.path(out_path, "UMAP_modules_immune.png"), width = 1200, height = 1000, res = 150) 

p<-plotEmbedding(proj,
              embedding = "UMAP",
              colorBy = "cellColData",
              name="nonTumor_Module.Immune",
              imputeWeights = getImputeWeights(proj),
              pal = my_pal)

print(p)
dev.off()

#########################################
#      Step 5. add peak set for immune  #
#########################################
# Create a new column with cell type labels
proj$CellType <- proj$Clusters  # copy existing clusters

# Define threshold
min_cells <- 50

# Get small cluster names
smallClusters <- names(table(proj$CellType)[table(proj$CellType) < min_cells])
smallClusters  # check which ones will be merged

# Create new column merging small clusters into "Other"
proj$CellType_merged <- as.character(proj$CellType)
proj$CellType_merged[proj$CellType %in% smallClusters] <- "Other"

# Verify
table(proj$CellType_merged)

proj <- addGroupCoverages(
  ArchRProj = proj,
  groupBy = "CellType_merged",
  minCells = 50,
  maxCells = 500,
  threads = 4,
  force = TRUE
)


pathToMacs2 <- '/data1/lesliec/carolw/my_envs/archr_2/bin/macs2'
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "CellType_merged",
  pathToMacs2 = pathToMacs2
)

proj <- addPeakMatrix(proj)


#########################################
#      Step 6. ChromVAR  #
#########################################

if("Motif" %ni% names(proj@peakAnnotation)){
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
}

# id bg peaks
tryCatch(
  proj <- addBgdPeaks(proj),
  error = function(e) {
    traceback()
    message(e)
  }
)


proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

markerTFs <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "MotifMatrix",
  groupBy = "CellType_merged",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

saveArchRProject(ArchRProj = proj)

motifs <- c("SPI1", "IRF8", "RUNX1","bcl11b")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs  # verify before plotting

p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "MotifMatrix",
  name = sort(markerMotifs),
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj),
  pal = my_pal
)


png(file.path(out_path, "UMAP_TFmotifs_immune.png"), 
    width = 2400, height = 1600, res = 150)
print(wrap_plots(p, ncol = 2))
dev.off()

