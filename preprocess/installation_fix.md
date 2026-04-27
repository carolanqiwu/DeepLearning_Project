For archR: /data1/lesliec/erinc/env_ymls/r_env.yml
Then in R: devtools::install_github("GreenleafLab/ArchR", ref="v1.0.2", repos = BiocManager::repositories())
Install ArchR Worked!
packageVersion("SeuratObject")                     
[1] ‘4.1.3’
Then issue with Seurat dependency is LSI
https://github.com/satijalab/seurat/issues/9169 : do not update anything they recommend
remotes::install_version("spatstat.geom", version = "3.2-1")
remotes::install_version("spatstat.utils", version = "3.1-0")
remotes::install_version("spatstat.data", version = "3.0-1")
remotes::install_version("spatstat.random", version = "3.1-5")
remotes::install_version("spatstat.sparse", version = "3.0-2")
remotes::install_version("spatstat.core", version = "2.4-4")

remotes::install_version(
  package = "Seurat",
  version = "4.1.1",
  repos = "https://cran.r-project.org"
)
Install Seurat worked!

For ArchR subsampling: 
# Optional subsampling of cells
if (cell_number!=0  && nCells(proj) > cell_number) {
        message("Subsampling to ", cell_number, " cells")
        cells <- sample(1:nCells(proj), cell_number, replace = FALSE)
        proj <- proj[cells, ]
} else {
        message("No subsampling — only ", nCells(proj), " cells")
}

