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
    # No header row: fall back to the documented column layouts by position, chosen automatically from the number of tab-separated fields
    ncol_detected <- length(first_fields)
    if (ncol_detected == 8L) {
      # chr start end CpG meth_pct meth unmeth alu_region
      idx_chr <- 1L; idx_end <- 3L; idx_cpg <- 4L; idx_alu <- 8L
    } else if (ncol_detected == 10L) {
      # chr start end CpG strand meth_pct meth unmeth compartment alu_region  (run134)
      idx_chr <- 1L; idx_end <- 3L; idx_cpg <- 4L; idx_alu <- 10L
    } else {
      # Unknown layout
      idx_chr <- 1L; idx_end <- 3L; idx_cpg <- 4L; idx_alu <- ncol_detected
      warning("Unexpected column count (", ncol_detected, ") with no header in ", filepath)
    }
    dt <- fread(filepath, sep = "\t", header = FALSE,select = c(idx_chr, idx_end, idx_cpg, idx_alu),
                na.strings = c("", "NA", "."), showProgress = FALSE)
    setnames(dt, c("chr", "pos", "CpG", "alu_region"))
  }

  dt[,chr:=as.character(chr)]
  dt[, pos := suppressWarnings(as.integer(pos))]
  dt[,CpG := as.character(CpG)]
  dt[,alu_region := as.character(alu_region)]
  dt[!is.na(pos)]
}
