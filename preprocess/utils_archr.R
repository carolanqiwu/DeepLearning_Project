getArchRInputFiles <- function(
  frag_dir,
  cohort
) {

  files <- list.files(
    path = frag_dir,
    pattern = paste0("^", cohort, "_.*\\.tsv\\.gz$"),
    full.names = TRUE
  )

  sample_names <- sub(
    paste0("^(", cohort, "_[^_]+).*"),
    "\\1",
    basename(files)
  )

  names(files) <- sample_names
  return(files)
}

