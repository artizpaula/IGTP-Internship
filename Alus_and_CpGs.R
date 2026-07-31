# Week 3
# Counts CpGs and unique Alu elements

# 1. Reading + format auto-detection

# Read one raw methylation table, auto-detecting the standard vs run134 layout
read_methylation_table <- function(filepath) {

  first_line <- readLines(filepath, n = 1L, warn = FALSE)
  if (length(first_line) == 0L || !nzchar(first_line)) {
    return(data.table(chr = character(0), pos = integer(0),
                       CpG = character(0), alu_region = character(0)))
  }
  first_fields <- trimws(strsplit(first_line, "\t", fixed = TRUE)[[1]])
  has_header <- tolower(first_fields[1]) == "chr"

  if (has_header) {
    # Identify columns by name, never by fixed position.
    header_lower <- tolower(first_fields)
    idx_chr <- match("chr", header_lower)
    idx_end <- match("end", header_lower)
    idx_cpg <- match("cpg", header_lower)
    idx_alu <- match("alu_region", header_lower)
    if (anyNA(c(idx_chr, idx_end, idx_cpg, idx_alu))) {
      stop("Could not find required column(s) chr/end/CpG/alu_region in header of: ", filepath)
    }
    dt <- fread(filepath, sep = "\t", header = TRUE,
                select = c(idx_chr, idx_end, idx_cpg, idx_alu),
                na.strings = c("", "NA", "."), showProgress = FALSE)
    setnames(dt, c("chr", "pos", "CpG", "alu_region"))

  } else {
    # No header row results in fall back to the documented column layouts by position
    ncol_detected <- length(first_fields)
    if (ncol_detected == 8L) {
      # chr start end CpG meth_pct meth unmeth alu_region
      idx_chr <- 1L; idx_end <- 3L; idx_cpg <- 4L; idx_alu <- 8L
    } else if (ncol_detected == 10L) {
      # chr start end CpG strand meth_pct meth unmeth compartment alu_region, (as in run134)
      idx_chr <- 1L; idx_end <- 3L; idx_cpg <- 4L; idx_alu <- 10L
    } else {
      # Fallback for unexpected layouts
      idx_chr <- 1L; idx_end <- 3L; idx_cpg <- 4L; idx_alu <- ncol_detected
      warning("Unexpected column count (", ncol_detected, ") with no header in ", filepath)
    }
    dt <- fread(filepath, sep = "\t", header = FALSE,select = c(idx_chr, idx_end, idx_cpg, idx_alu),
                na.strings = c("", "NA", "."), showProgress = FALSE)
    setnames(dt, c("chr", "pos", "CpG", "alu_region"))
  }

  # Ensure consistent column types
  dt[,chr:=as.character(chr)]
  dt[, pos := suppressWarnings(as.integer(pos))]
  dt[,CpG := as.character(CpG)]
  dt[,alu_region := as.character(alu_region)]
  dt[!is.na(pos)]
}

# 2. Binning + per-sample summary

# Assign each row to a 1 Mb bin based in position
assign_bins <- function(dt, bin_size = 1e6) {
  bin_size <- as.integer(bin_size)
  # Compute genomic bin boundaries
  dt[, bin_start := (pos - 1L) %/% bin_size * bin_size + 1L]
  dt[, bin_end   := bin_start + bin_size - 1L]
  dt
}

# Summarize one sample's already-read methylation table into per-bin counts
summarize_bins <- function(dt, bin_size = 1e6) {
  if (nrow(dt) == 0L) {
    return(data.table(chr = character(0), bin_start = integer(0), bin_end = integer(0), n_CpGs = integer(0), n_Alus = integer(0)))}
  dt <- assign_bins(dt, bin_size)
  
  # counting CpGs and unique Alu elements per bin
  out <- dt[, .(
    n_CpGs = sum(!is.na(CpG) & nzchar(CpG)),
    n_Alus = uniqueN(alu_region[!is.na(alu_region) & nzchar(alu_region)])
  ), by = .(chr, bin_start, bin_end)]
  setorder(out, chr, bin_start)
  out[]
}

# End-to-end per-sample pipeline
bin_methylation_file <- function(filepath, bin_size = 1e6) {
  dt <- read_methylation_table(filepath)
  summarize_bins(dt, bin_size = bin_size)
}

# 3. Sample discovery: loose files and/or samples bundled inside .tar.gz archives

#Turn a raw file/entry name into a short, readable sample ID.
derive_sample_id <- function(entry_path) {
  fname <- basename(entry_path)
  fname <- sub("\\.(tsv|txt)$", "", fname, ignore.case = TRUE)
  fname <- sub("\\.RGok.*$", "", fname) # strip long aligner/QC suffix chain
  fname <- sub("_collapsed$", "", fname)
  fname <- sub("_sorted$", "", fname)
  fname
}

# Process every sample file bundled inside one .tar.gz archive
process_tar_archive <- function(tar_path, bin_size = 1e6, exdir = tempdir(), cleanup = TRUE) {
  entries <- utils::untar(tar_path, list = TRUE)
  entries <- entries[grepl("\\.(tsv|txt)$", entries, ignore.case = TRUE)]
  results <- list()
  for (entry in entries) {
    sample_id <- derive_sample_id(entry)
    message("  -> ", sample_id)
    utils::untar(tar_path, files = entry, exdir = exdir)
    extracted_path <- file.path(exdir, entry)
    results[[sample_id]] <- bin_methylation_file(extracted_path, bin_size = bin_size)
    if (cleanup) unlink(extracted_path)
  }
  results
}

# Scan a directory for raw methylation input: .tar.gz archives
list_raw_methylation_files <- function(input_dir, tar_pattern = "\\.tar\\.gz$|_tar\\.gz$|\\.tgz$", flat_pattern = "\\.tsv$|\\.txt$") {
  tar_files  <- list.files(input_dir, pattern = tar_pattern, full.names = TRUE, ignore.case = TRUE)
  flat_files <- list.files(input_dir, pattern = flat_pattern, full.names = TRUE, ignore.case = TRUE)
  list(tar_files = tar_files, flat_files = flat_files)
}

# Main entry point: process every raw sample found in input_dir
 
process_all_raw_methylation <- function(input_dir, output_dir = NULL, bin_size = 1e6,
                                         cleanup = TRUE) {
  found <- list_raw_methylation_files(input_dir)
  all_results <- list()

  # Process compressed files
  for (tar_path in found$tar_files) {
    message("Processing archive: ", basename(tar_path))
    res <- process_tar_archive(tar_path, bin_size = bin_size, cleanup = cleanup)
    all_results[names(res)] <- res
  }

  # Process standalone files
  for (flat_path in found$flat_files) {
    sample_id <- derive_sample_id(flat_path)
    message("Processing file: ", sample_id)
    all_results[[sample_id]] <- bin_methylation_file(flat_path, bin_size = bin_size)
  }

  # Save results if an output directory was provided
  if (!is.null(output_dir) && length(all_results) > 0) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    for (sample_id in names(all_results)) {
      safe_name <- gsub("[^A-Za-z0-9_.-]", "_", sample_id)
      out_path <- file.path(output_dir, paste0(safe_name, "_cpg_alu_bins.csv"))
      fwrite(all_results[[sample_id]], out_path)
    }
  }

  all_results
}

# 4. Aggregating per-sample tables into one bin-level annotation

aggregate_bin_annotation <- function(sample_tables, fun = c("max", "sum", "mean", "median")) {
  fun <- match.arg(fun)
  sample_tables <- sample_tables[vapply(sample_tables, nrow, integer(1)) > 0] # ignore empty samples
  if (length(sample_tables) == 0) {
    return(data.table(chr = character(0), bin_start = integer(0), bin_end = integer(0),
                       n_cpg = integer(0), n_alu = integer(0)))
  }
  combined <- rbindlist(sample_tables, idcol = "sample_id")
  agg_fun <- switch(fun,
    max= function(x) max(x, na.rm = TRUE),
    sum= function(x) sum(x, na.rm = TRUE),
    mean= function(x) mean(x, na.rm = TRUE),
    median= function(x) stats::median(x, na.rm = TRUE)
  )
  # Aggregate counts across samples
  out <- combined[, .(
    n_cpg = as.integer(round(agg_fun(n_CpGs))),
    n_alu = as.integer(round(agg_fun(n_Alus)))), 
    by = .(chr, bin_start, bin_end)]
  setorder(out, chr, bin_start)
  out[ ]
}
