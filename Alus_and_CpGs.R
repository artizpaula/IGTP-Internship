# Week 3
# Counts CpGs and unique Alu elements

# 1. Reading + format auto-detection

# reads one raw methylation table and figures out on its own whether it's the standard layout or the run134 layout
read_methylation_table <- function(filepath) {
  
  first_line <- readLines(filepath, n = 1L, warn = FALSE)
  if (length(first_line) == 0L || !nzchar(first_line)) {
    return(data.table(chr = character(0), pos = integer(0), CpG = character(0), alu_region = character(0), meth = numeric(0))) # empty file -> empty table but with the right columns
  }
  first_fields <- trimws(strsplit(first_line, "\t", fixed = TRUE)[[1]]) # split first line into columns
  has_header <- tolower(first_fields[1]) == "chr" # if the first entry says "chr" we assume there's a header row
  
  if (has_header) {
    # match columns by name instead of position, safer this way
    header_lower <- tolower(first_fields) # lowercase the header so matching isn't case sensitive
    idx_chr <- match("chr", header_lower) # which column has the chromosome
    idx_end <- match("end", header_lower) # which column has the end position
    idx_cpg <- match("cpg", header_lower) # which column has the CpG values
    idx_alu <- match("alu_region", header_lower) # which column has the alu info
    idx_meth <- match("meth", header_lower)
    if (is.na(idx_meth)) idx_meth <- match("meth_pct", header_lower)
    if (anyNA(c(idx_chr, idx_end, idx_cpg, idx_alu))) {
      stop("Could not find required column(s) chr/end/CpG/alu_region in header of: ", filepath)
    }
    select_idx <- c(idx_chr, idx_end, idx_cpg, idx_alu, if (!is.na(idx_meth)) idx_meth)
    col_names  <- c("chr", "pos", "CpG", "alu_region", if (!is.na(idx_meth)) "meth")
    dt <- fread(filepath, sep = "\t", header = TRUE, select = select_idx,na.strings = c("", "NA", "."), showProgress = FALSE)
    setnames(dt, col_names) # rename to something consistent, makes the rest of the pipeline simpler
    if (is.na(idx_meth)) dt[, meth := NA_real_]
    
  } else {
    # no header row, so fall back to the documented column layouts by position
    ncol_detected <- length(first_fields)
    if (ncol_detected == 8L) {
      # chr start end CpG meth_pct meth unmeth alu_region
      idx_chr <- 1L; idx_end <- 3L; idx_cpg <- 4L; idx_meth <- 6L; idx_alu <- 8L
    } else if (ncol_detected == 10L) {
      # chr start end CpG strand meth_pct meth unmeth compartment alu_region (run134 layout)
      idx_chr <- 1L; idx_end <- 3L; idx_cpg <- 4L; idx_meth <- 7L; idx_alu <- 10L
    } else {
      idx_chr <- 1L; idx_end <- 3L; idx_cpg <- 4L; idx_meth <- NA_integer_; idx_alu <- ncol_detected
      warning("Unexpected column count (", ncol_detected, ") with no header in ", filepath) # fallback for weird layouts
    }
    select_idx <- c(idx_chr, idx_end, idx_cpg, idx_alu, if (!is.na(idx_meth)) idx_meth)
    col_names  <- c("chr", "pos", "CpG", "alu_region", if (!is.na(idx_meth)) "meth")
    dt <- fread(filepath, sep = "\t", header = FALSE, select = select_idx, na.strings = c("", "NA", "."), showProgress = FALSE)
    setnames(dt, col_names)
    if (is.na(idx_meth)) dt[, meth := NA_real_]
  }
  
  # make sure the column types are consistent no matter which branch above ran
  dt[,chr:=as.character(chr)] # "1" -> "chr1"
  dt[, pos := suppressWarnings(as.integer(pos))]
  dt[,CpG := as.character(CpG)]
  dt[,alu_region := as.character(alu_region)]
  dt[, meth := suppressWarnings(as.numeric(meth))] # keep each CpG-Alu row's own methylation value
  dt[!is.na(pos)]
}

# 2. Binning + per-sample summary

# assigns every row to a 1 Mb bin based on its position
assign_bins <- function(dt, bin_size = 1e6) {
  bin_size <- as.integer(bin_size)
  # work out the genomic bin boundaries
  dt[, bin_start := (pos - 1L) %/% bin_size * bin_size + 1L] # get position within the bin size
  dt[, bin_end   := bin_start + bin_size - 1L]
  dt
}

# takes one sample's already-read methylation table and summarizes it into per-bin counts
summarize_bins <- function(dt, bin_size = 1e6) {
  if (nrow(dt) == 0L) {
    return(data.table(chr = character(0), bin_start = integer(0), bin_end = integer(0), n_CpGs = integer(0), n_Alus = integer(0),  mean_meth = numeric(0), sum_meth = numeric(0), pct_meth = numeric(0)))}
  dt <- assign_bins(dt, bin_size)
  
  has_cpg <- !is.na(dt$CpG) & nzchar(dt$CpG)
  has_alu <- !is.na(dt$alu_region) & nzchar(dt$alu_region)
  
  # step 1: average methylation of each individual Alu region (its own CpGs only), per bin
  per_alu <- dt[has_cpg & has_alu, .(alu_mean_meth = mean(meth, na.rm = TRUE) # e.g. Alu A: mean(80,50,100) = 76.7%; Alu B: mean(0,60) = 30%
  ), by = .(chr, bin_start, bin_end, alu_region)]
  
  # step 2: bin-level counts/sums, still computed per CpG-Alu observation
  bin_counts <- dt[has_cpg, .(n_CpGs = .N, # every row in the bin with a non-empty CpG counts once, i.e. once per CpG-Alu pair (a CpG in 2 Alus counts twice)
    n_Alus = uniqueN(alu_region[!is.na(alu_region) & nzchar(alu_region)]), # distinct alu regions touched in the bin (doesn't change n_CpGs/meth above)
    sum_meth = sum(meth, na.rm = TRUE),  # total methylation signal across every individual CpG-Alu observation
    pct_meth  = 100 * sum(!is.na(meth) & meth > 0) / sum(!is.na(meth)) # % of CpG-Alu observations with a positive methylation value
  ), by = .(chr, bin_start, bin_end)]
  
  # step 3: bin's mean_meth = average of each Alu's own average, so every Alu counts once no matter how many CpGs it has (Alu A 76.7%, Alu B 30% -> bin = (76.7+30)/2 = 53.3%)
  alu_bin_means <- per_alu[, .(mean_meth = mean(alu_mean_meth, na.rm = TRUE)), by = .(chr, bin_start, bin_end)]
  bin_counts <- merge(bin_counts, alu_bin_means, by = c("chr", "bin_start", "bin_end"), all.x = TRUE)
  setcolorder(bin_counts, c("chr", "bin_start", "bin_end", "n_CpGs", "n_Alus", "mean_meth", "sum_meth", "pct_meth"))
  
  # clean up any NaNs that pop out of empty means
  
  bin_counts[is.nan(mean_meth), mean_meth := NA_real_]
  bin_counts[is.nan(pct_meth), pct_meth := NA_real_]
  setorder(bin_counts, chr, bin_start)
  bin_counts[]
}

# full pipeline for one sample, start to finish
bin_methylation_file <- function(filepath, bin_size = 1e6) {
  dt <- read_methylation_table(filepath)
  summarize_bins(dt, bin_size = bin_size)
}

# 3. Sample discovery: loose files and/or samples bundled inside .tar.gz archives

# turns a raw file/entry name into a short, readable sample ID
derive_sample_id <- function(entry_path) {
  fname <- basename(entry_path)
  fname <- sub("\\.(tsv|txt)$", "", fname, ignore.case = TRUE) # drop .tsv or .txt extension
  fname <- sub("\\.RGok.*$", "", fname) # trim off anything after .RGok
  fname <- sub("_collapsed$", "", fname) # remove _collapsed from the name
  fname <- sub("_sorted$", "", fname) # remove _sorted from the name
  fname
}

# processes every sample file bundled inside one .tar.gz archive
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

# scans a directory for raw methylation input: loose files and/or .tar.gz archives
list_raw_methylation_files <- function(input_dir, tar_pattern = "\\.tar\\.gz$|_tar\\.gz$|\\.tgz$", flat_pattern = "\\.tsv$|\\.txt$") { 
  tar_files  <- list.files(input_dir, pattern = tar_pattern, full.names = TRUE, ignore.case = TRUE) # any compressed archives in the folder
  flat_files <- list.files(input_dir, pattern = flat_pattern, full.names = TRUE, ignore.case = TRUE) # any loose files sitting directly in the folder
  list(tar_files = tar_files, flat_files = flat_files) # bundle both together
}

# main entry point: processes every raw sample found in input_dir

process_all_raw_methylation <- function(input_dir, output_dir = NULL, bin_size = 1e6, cleanup = TRUE) {
  found <- list_raw_methylation_files(input_dir)
  all_results <- list()
  
  # process compressed archives first
  for (tar_path in found$tar_files) {
    message("Processing archive: ", basename(tar_path))
    res <- process_tar_archive(tar_path, bin_size = bin_size, cleanup = cleanup)
    all_results[names(res)] <- res
  }
  
  # then the standalone files
  for (flat_path in found$flat_files) {
    sample_id <- derive_sample_id(flat_path)
    message("Processing file: ", sample_id)
    all_results[[sample_id]] <- bin_methylation_file(flat_path, bin_size = bin_size)
  }
  
  # save everything out if an output folder was given
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
  sample_tables <- sample_tables[vapply(sample_tables, nrow, integer(1)) > 0] # skip any empty samples
  if (length(sample_tables) == 0) {
    return(data.table(chr = character(0), bin_start = integer(0), bin_end = integer(0), n_cpg = integer(0), n_alu = integer(0), mean_meth = numeric(0), pct_meth = numeric(0)))
  }
  combined <- rbindlist(sample_tables, idcol = "sample_id", fill = TRUE)
  agg_fun <- switch(fun, max = function(x) max(x, na.rm = TRUE),# largest value across samples for a bin
                    sum = function(x) sum(x, na.rm = TRUE),# adds values up across samples for a bin
                    mean = function(x) mean(x, na.rm = TRUE),# averages values across samples for a bin
                    median = function(x) stats::median(x, na.rm = TRUE) # middle value across samples for a bin
  )
  
  bin_summary <- combined[, .(
    n_cpg = as.integer(round(agg_fun(n_CpGs))), # combine each sample's CpG-Alu count for this bin into one number
    n_alu = as.integer(round(agg_fun(n_Alus))), # same idea but for distinct alus
    mean_meth = agg_fun(mean_meth), # combine each sample's per-observation average methylation for this bin
    pct_meth  = agg_fun(pct_meth)), # combine each sample's % methylated CpG-Alu observations for this bin
    by = .(chr, bin_start, bin_end)]
  setorder(bin_summary, chr, bin_start)  # sort the final table
  bin_summary[ ]
}