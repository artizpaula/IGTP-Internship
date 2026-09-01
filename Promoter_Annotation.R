# Promoter annotation: match promoter reference coordinates against bins.

# Before running: download a BED file from UCSC Table Browser (hg19, track = NCBI RefSeq, output format = BED, "upstream by 1000 bases"), then convert
# it to a TSV with columns promoter_name / chr / start / end and save it as Promoter_reference_GRCh37.tsv in your Archivo dataset folder.

# 1. Read the promoter reference file

read_promoter_reference <- function(promoter_tsv_path) {
  promoters <- read.delim(promoter_tsv_path, sep = "\t", stringsAsFactors = FALSE,
                          check.names = FALSE, quote = "")
  required_cols <- c("promoter_name", "chr", "start", "end")
  missing_cols <- setdiff(required_cols, names(promoters))
  if (length(missing_cols) > 0) {
    stop("Promoter reference is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  promoters$promoter_name <- trimws(promoters$promoter_name)
  promoters$chr <- trimws(as.character(promoters$chr))
  promoters$chr <- sub("^chr", "", promoters$chr, ignore.case = TRUE) # normalize "chr1" -> "1" to match bin_table$chr
  promoters$start <- suppressWarnings(as.numeric(promoters$start))
  promoters$end   <- suppressWarnings(as.numeric(promoters$end))
  
  usable <- !is.na(promoters$start) & !is.na(promoters$end) &
    nzchar(promoters$promoter_name) & nzchar(promoters$chr)
  skipped_missing_coords <- sum(!usable)
  promoter_reference <- promoters[usable, c("promoter_name", "chr", "start", "end")]
  
  # guard against reversed start/end so the overlap rule (start <= end) always holds
  swap <- promoter_reference$start > promoter_reference$end
  if (any(swap)) {
    tmp <- promoter_reference$start[swap]
    promoter_reference$start[swap] <- promoter_reference$end[swap]
    promoter_reference$end[swap] <- tmp
  }
  
  attr(promoter_reference, "skipped_rows") <- skipped_missing_coords
  promoter_reference
}

# 2. Bin <-> promoter coordinate-overlap join

annotate_bins_with_promoters <- function(bin_table, promoter_reference) {
  required_bin_cols <- c("bin_id", "chr", "bin_start", "bin_end")
  missing_bin_cols <- setdiff(required_bin_cols, names(bin_table))
  if (length(missing_bin_cols) > 0) {
    stop("bin_table is missing required column(s): ", paste(missing_bin_cols, collapse = ", "))
  }
  required_prom_cols <- c("promoter_name", "chr", "start", "end")
  missing_prom_cols <- setdiff(required_prom_cols, names(promoter_reference))
  if (length(missing_prom_cols) > 0) {
    stop("promoter_reference is missing required column(s): ", paste(missing_prom_cols, collapse = ", "))
  }
  
  bins <- bin_table[, c("bin_id", "chr", "bin_start", "bin_end")]
  bins$chr <- trimws(as.character(bins$chr))
  proms <- promoter_reference
  proms$chr <- trimws(as.character(proms$chr))
  
  # Cartesian join within each chromosome
  candidates <- merge(bins, proms, by = "chr", all = FALSE)
  overlap <- candidates$start <= candidates$bin_end & candidates$end >= candidates$bin_start
  promoter_hits <- candidates[overlap, c("bin_id", "chr", "bin_start", "bin_end",
                                         "promoter_name", "start", "end")]
  rownames(promoter_hits) <- NULL
  promoter_hits <- promoter_hits[order(promoter_hits$bin_id, promoter_hits$promoter_name), ]
  
  if (nrow(promoter_hits) > 0) {
    promoter_summary <- aggregate(promoter_name ~ bin_id, data = promoter_hits,
                                  FUN = function(x) length(unique(x)))
    names(promoter_summary)[2] <- "promoter_count"
  } else {
    promoter_summary <- data.frame(bin_id = character(0), promoter_count = integer(0))
  }
  
  # drop any old/stale version of this column before merging in the new one
  bin_table$promoter_count <- NULL
  
  bin_table <- merge(bin_table, promoter_summary, by = "bin_id", all.x = TRUE, sort = FALSE)
  bin_table$promoter_count[is.na(bin_table$promoter_count)] <- 0L
  bin_table$promoter_count <- as.integer(bin_table$promoter_count)
  
  list(bin_table = bin_table, promoter_hits = promoter_hits)
}

# 3. Validation

# Cross-checks the annotation for internal consistency and prints a short summary of coverage
validate_promoter_bin_annotation <- function(bin_table, promoter_reference, promoter_hits) {
  bin_chr_lookup <- setNames(as.character(bin_table$chr), bin_table$bin_id)
  bad_chr <- promoter_hits$chr != bin_chr_lookup[promoter_hits$bin_id]
  if (any(bad_chr)) {
    stop(sum(bad_chr), " promoter-bin match(es) span different chromosomes - this should never happen.")
  }
  
  bad_overlap <- !(promoter_hits$start <= promoter_hits$bin_end & promoter_hits$end >= promoter_hits$bin_start)
  if (any(bad_overlap)) {
    stop(sum(bad_overlap), " promoter-bin pair(s) violate the overlap rule ",
         "(start <= bin_end AND end >= bin_start).")
  }
  
  n_bins <- nrow(bin_table)
  n_bins_with_promoter <- sum(bin_table$promoter_count > 0)
  n_ref_promoters <- nrow(promoter_reference)
  matched_promoters <- unique(promoter_hits$promoter_name)
  
  message("Promoter-bin annotation validated: ", nrow(promoter_hits),
          " overlapping pair(s) found, 0 rule violations.")
  message(n_bins_with_promoter, " / ", n_bins, " bins have >=1 promoter; ",
          length(matched_promoters), " / ", n_ref_promoters, " reference promoters matched >=1 bin.")
  skipped <- attr(promoter_reference, "skipped_rows")
  if (!is.null(skipped) && skipped > 0) {
    message(skipped, " row(s) in the promoter reference were skipped (missing name/coordinates).")
  }
  
  invisible(list(n_bins = n_bins, n_bins_with_promoter = n_bins_with_promoter,
                 n_reference_promoters = n_ref_promoters, n_matched_promoters = length(matched_promoters)))
}