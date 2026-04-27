.libPaths('/data1/lesliec/carolw/my_envs/archr_2/lib/R/library')
library(ArchR)
library(readxl)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)


genome_assembly <- "hg38"
addArchRGenome(genome_assembly)
addArchRThreads(threads = 4)

#archr_path <- "/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/archr_data_subsampled"
archr_path <- "/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/archr_data_allcells"
proj <- loadArchRProject(
  path = archr_path,
  showLogo = FALSE
)

#####################################
#            immune cells           #
#####################################
immune_cells <- proj$Clusters == "C2" #GBM

immune_cell_names <- proj$cellNames[immune_cells]
#immune_cell_names <- rownames(proj@cellColData)[proj$Clusters %in% c("C1","C2","C3","C4")] #BRCA

nonimmune_cell_names <- setdiff(proj$cellNames, immune_cell_names)

# skip
proj_immune <- subsetArchRProject(
  ArchRProj = proj,
  cells = immune_cell_names,
  outputDirectory = "/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/archr_data_allcells_immune",
  dropCells = TRUE,
  force=TRUE
)

final_bc <- rownames(proj_immune@cellColData)
write.csv(final_bc, file = paste0("/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/archr_data_subsampled_immune",  "/archr_filtered_barcode.csv"))

lsi <- getReducedDims(proj_immune, reducedDims = "Harmony")
write.csv(lsi, file = paste0("/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/archr_data_subsampled_immune", "/archr_filtered_lsi.csv"))

#save fragment
cell_names <- rownames(proj_immune@cellColData)

frags <- getFragmentsFromProject(
  ArchRProj = proj_immune,
  cellNames = immune_cell_names
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
  file = gzfile(file.path("/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/archr_data_subsampled_immune", "immuneCells_fragments.tsv.gz")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)



#####################################
#          nonimmune cells          #
#####################################

proj_nonimmune <- subsetArchRProject(
  ArchRProj = proj,
  cells = nonimmune_cell_names,
  outputDirectory = "/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA-chan/archr_data_subsampled_nonimmune",
  dropCells = TRUE
)

# resume from cluster analysis
proj_nonimmune <- loadArchRProject(
  path = archr_path,
  showLogo = FALSE
)
# Extract the module scores
scores <- proj_nonimmune@cellColData[, c(
  "GBM_subtype_Module.Proneural",
  "GBM_subtype_Module.Neural",
  "GBM_subtype_Module.Classical",
  "GBM_subtype_Module.Mesenchymal"
)]

scores <- proj_nonimmune@cellColData[, c(
  "BRCA_subtype_Module.Luminal",
  "BRCA_subtype_Module.HER2",
  "BRCA_subtype_Module.Basal"
)]
# Assign each cell to the subtype with the highest score
cell_subtype <- apply(scores, 1, function(x) {
  names(x)[which.max(x)]
})

# Clean the names (optional: remove "GBM_subtype_" prefix)
cell_subtype <- sub("GBM_subtype_Module\\.", "", cell_subtype)
cell_subtype <- sub("BRCA_subtype_Module\\.", "", cell_subtype)

# Add to cellColData
proj_nonimmune$Subtype <- cell_subtype
table(proj_nonimmune$Subtype)
plotEmbedding(
  ArchRProj = proj_nonimmune,
  colorBy = "cellColData",
  name = "Subtype",      # Column containing your subtype labels
  embedding = "UMAP"     # Or "TSNE" if you prefer
)

setwd("/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA-chan/archr_data_subsampled_nonimmune")
# Export LSI, barcode, fragments for individual subtypes
subtypes <- c("Proneural", "Neural", "Classical", "Mesenchymal")
subtypes <- c("Luminal","HER2","Basal")
  
for(st in subtypes){
  
  # Cells of this subtype
  cells <- rownames(proj_nonimmune@cellColData)[proj_nonimmune$Subtype == st]
  
  # Export barcodes
  write.csv(cells, file = paste0(st, "_barcodes.csv"), row.names = FALSE)
  
  # Export LSI/Harmony embeddings
  lsi_subset <- getReducedDims(proj_nonimmune, "Harmony")[cells, ]
  write.csv(lsi_subset, file = paste0(st, "_lsi.csv"))
  
  # Export fragment file for these cells
  frags <- getFragmentsFromProject(
    ArchRProj = proj_nonimmune,
    cellNames = cells
  )
  
  frags_unlisted <- unlist(frags)
  
  frags_df <- data.frame(
    seqnames = seqnames(frags_unlisted),
    start = start(frags_unlisted),
    end = end(frags_unlisted),
    RG = frags_unlisted$RG,
    score = 1
  )
  
  write.table(
    frags_df,
    file = gzfile(paste0("NonImmune_", st, "_fragments.tsv.gz")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
}

#####################################
#         subset subtype to 2k          #
#####################################
scores <- proj_nonimmune@cellColData[, c(
  "GBM_subtype_Module.Proneural",
  "GBM_subtype_Module.Neural",
  "GBM_subtype_Module.Classical",
  "GBM_subtype_Module.Mesenchymal"
)]

cell_subtype <- apply(scores, 1, function(x) {
  names(x)[which.max(x)]
})

# Clean the names (optional: remove "GBM_subtype_" prefix)
cell_subtype <- sub("GBM_subtype_Module\\.", "", cell_subtype)
#cell_subtype <- sub("BRCA_subtype_Module\\.", "", cell_subtype)

# Add to cellColData
proj_nonimmune$Subtype <- cell_subtype

setwd("/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/archr_data_subsampled_nonimmune")
# Export LSI, barcode, fragments for individual subtypes
subtypes <- c("Proneural", "Neural", "Classical", "Mesenchymal")

set.seed(42)

for(st in c("Classical", "Proneural")){
  
  cells <- rownames(proj_nonimmune@cellColData)[proj_nonimmune$Subtype == st]
  
  # Subsample to 2k if more than 2000 cells
  if(length(cells) > 2000){
    cells <- sample(cells, 2000)
  }
  
  # Export barcodes
  write.csv(cells, file = paste0(st, "_barcodes_2k.csv"), row.names = FALSE)
  
  # Export LSI/Harmony embeddings
  lsi_subset <- getReducedDims(proj_nonimmune, "Harmony")[cells, ]
  write.csv(lsi_subset, file = paste0(st, "_lsi_2k.csv"))
  
  # Export fragment file for these cells
  frags <- getFragmentsFromProject(
    ArchRProj = proj_nonimmune,
    cellNames = cells
  )
  
  frags_unlisted <- unlist(frags)
  
  frags_df <- data.frame(
    seqnames = seqnames(frags_unlisted),
    start = start(frags_unlisted),
    end = end(frags_unlisted),
    RG = frags_unlisted$RG,
    score = 1
  )
  
  write.table(
    frags_df,
    file = gzfile(paste0("NonImmune_",st,"_2k_fragments.tsv.gz")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
}

#####################################
#       subset w/o proj dir          #
#####################################
cellsUse <- nonimmune_cell_names

scores <- proj@cellColData[cellsUse, c(
  "Subtype_Module.MES",
  "Subtype_Module.CL",
  "Subtype_Module.NL",
  "Subtype_Module.PN"
)] %>% as.data.frame()

cell_subtype <- apply(scores, 1, function(x) {
  names(x)[which.max(x)]
})

# clean names
cell_subtype <- sub("Subtype_Module\\.", "", cell_subtype)
names(cell_subtype) <- rownames(scores)

# Build a full-length vector (NA for everyone)
all_cells <- rownames(proj@cellColData)
subtype_full <- rep(NA_character_, length(all_cells))
names(subtype_full) <- all_cells

# Fill in only the cells you care about
subtype_full[cellsUse] <- cell_subtype

# Assign the whole column at once — no Rle indexing issues
proj$Subtype <- subtype_full

plotEmbedding(
  ArchRProj = proj,
  cells = cellsUse,
  colorBy = "cellColData",
  name = "Subtype",
  embedding = "UMAP"
)