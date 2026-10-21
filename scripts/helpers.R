# Helper functions shared by the SoMoSeq analysis notebooks

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Read config/paths.example.yaml

load_paths <- function(config_dir = file.path("..", "config")) {
  f <- file.path(config_dir, "paths.example.yaml")
  if (!file.exists(f)) f <- file.path(config_dir, "paths_example.yaml")
  p <- yaml::read_yaml(f)
  base <- normalizePath(config_dir)
  lapply(p, function(x) {
    if (is.character(x) && length(x) == 1L && nzchar(x) && !grepl("^(/|~)", x))
      file.path(base, x) else x
  })
}

# write tsv files
save_tsv <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(df, file = path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "NA")
}

# create counts and logcounts with unique gene-symbol row names

symbol_matrices <- function(sce, symbol_field = "symbol") {
  syms <- as.character(SummarizedExperiment::rowData(sce)[[symbol_field]])
  bad <- is.na(syms) | !nzchar(trimws(syms))
  syms[bad] <- rownames(sce)[bad]
  syms <- make.unique(syms)
  counts_sym <- SummarizedExperiment::assay(sce, "counts")
  log_sym    <- SummarizedExperiment::assay(sce, "logcounts")
  rownames(counts_sym) <- syms
  rownames(log_sym)    <- syms
  list(counts = counts_sym, logcounts = log_sym)
}
