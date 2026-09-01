# Week 4
# Gene annotation: match Cancer Gene Census (COSMIC) coordinates

read_gene_reference <- function(cosmic_tsv_path, crc_gene_list_path = NULL) {
  cosmic <- read.delim(cosmic_tsv_path, sep = "\t", stringsAsFactors = FALSE,
                       check.names = FALSE, quote = "")
  required_cols <- c("GENE_SYMBOL", "COSMIC_GENE_ID", "CHROMOSOME", "GENOME_START", "GENOME_STOP")
  missing_cols <- setdiff(required_cols, names(cosmic))
  if (length(missing_cols) > 0) {
    stop("COSMIC census is missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  cosmic$GENE_SYMBOL<- trimws(cosmic$GENE_SYMBOL)
  cosmic$CHROMOSOME <- trimws(as.character(cosmic$CHROMOSOME))
  cosmic$GENOME_START <- suppressWarnings(as.numeric(cosmic$GENOME_START))
  cosmic$GENOME_STOP  <- suppressWarnings(as.numeric(cosmic$GENOME_STOP))
  
  use_curated_subset <- !is.null(crc_gene_list_path) && !is.na(crc_gene_list_path) &&
    nzchar(crc_gene_list_path)
  
  if (use_curated_subset) {
    wanted_genes <- readLines(crc_gene_list_path, warn = FALSE)
    wanted_genes <- trimws(wanted_genes)
    wanted_genes <- wanted_genes[nzchar(wanted_genes)]
    if (length(wanted_genes) > 0 && toupper(wanted_genes[1]) == "GENE") wanted_genes <- wanted_genes[-1] # drop header row if there is one
    wanted_genes <- unique(toupper(wanted_genes))
  } else {
    # full census: every distinct gene symbol in the COSMIC table
    wanted_genes <- unique(toupper(cosmic$GENE_SYMBOL))
  }
  
  cosmic_upper<- toupper(cosmic$GENE_SYMBOL)
  matched_idx <- match(wanted_genes, cosmic_upper)
  
  found <- !is.na(matched_idx)
  has_coords <- found
  has_coords[found] <- !is.na(cosmic$GENOME_START[matched_idx[found]]) &
    !is.na(cosmic$GENOME_STOP[matched_idx[found]])
  
  usable <- found & has_coords
  ref_idx <- matched_idx[usable]
  
  gene_reference <- data.frame(gene_id = cosmic$COSMIC_GENE_ID[ref_idx],
    gene_name = cosmic$GENE_SYMBOL[ref_idx],
    chr= cosmic$CHROMOSOME[ref_idx],
    gene_start = cosmic$GENOME_START[ref_idx],
    gene_end  = cosmic$GENOME_STOP[ref_idx],
    stringsAsFactors = FALSE)
  
  # Guard against reversed start/end so the overlap rule (which assumes gene_start <= gene_end) always holds
  swap <- gene_reference$gene_start > gene_reference$gene_end
  if (any(swap)) {
    tmp <- gene_reference$gene_start[swap]
    gene_reference$gene_start[swap] <- gene_reference$gene_end[swap]
    gene_reference$gene_end[swap] <- tmp
  }
  
  attr(gene_reference, "skipped_genes") <- list(
    not_in_census = wanted_genes[!found],
    missing_coordinates = wanted_genes[found & !has_coords])
  
  gene_reference
}

# 2. Bin <-> gene coordinate-overlap join

annotate_bins_with_genes <- function(bin_table, gene_reference, id_sep = ";") {
  required_bin_cols <- c("bin_id", "chr", "bin_start", "bin_end")
  missing_bin_cols <- setdiff(required_bin_cols, names(bin_table))
  if (length(missing_bin_cols) > 0) {
    stop("bin_table is missing required column(s): ", paste(missing_bin_cols, collapse = ", "))
  }
  required_gene_cols <- c("gene_id", "gene_name", "chr", "gene_start", "gene_end")
  missing_gene_cols <- setdiff(required_gene_cols, names(gene_reference))
  if (length(missing_gene_cols) > 0) {
    stop("gene_reference is missing required column(s): ", paste(missing_gene_cols, collapse = ", "))
  }
  
  bins <- bin_table[, c("bin_id", "chr", "bin_start", "bin_end")]
  bins$chr <- trimws(as.character(bins$chr))
  genes <- gene_reference
  genes$chr <- trimws(as.character(genes$chr))
  
  # Cartesian join within each chromosome
  candidates <- merge(bins, genes, by = "chr", all = FALSE)
  overlap <- candidates$gene_start <= candidates$bin_end & candidates$gene_end >= candidates$bin_start
  gene_hits <- candidates[overlap, c("bin_id", "chr", "bin_start", "bin_end",
                                     "gene_id", "gene_name", "gene_start", "gene_end")]
  rownames(gene_hits) <- NULL
  gene_hits <- gene_hits[order(gene_hits$bin_id, gene_hits$gene_name), ]
  
  if (nrow(gene_hits) > 0) {
    per_bin <- aggregate(gene_name ~ bin_id, data = gene_hits,
                         FUN = function(x) length(unique(x)))
    names(per_bin)[2] <- "gene_count"
    ids_per_bin <- aggregate(gene_id ~ bin_id, data = gene_hits,
                             FUN = function(x) paste(sort(unique(x)), collapse = id_sep))
    names(ids_per_bin)[2] <- "gene_ids"
    names_per_bin <- aggregate(gene_name ~ bin_id, data = gene_hits,
                               FUN = function(x) paste(sort(unique(x)), collapse = id_sep))
    names(names_per_bin)[2] <- "gene_names"
    gene_summary <- Reduce(function(a, b) merge(a, b, by = "bin_id", all = TRUE),
                           list(per_bin, ids_per_bin, names_per_bin))
  } else {
    gene_summary <- data.frame(bin_id = character(0), gene_count = integer(0),
                               gene_ids = character(0), gene_names = character(0))}
  
  # drop any old/stale versions of these columns before merging in the new ones
  bin_table$gene_count <- NULL
  bin_table$gene_ids <- NULL
  bin_table$gene_names <- NULL
  
  bin_table <- merge(bin_table, gene_summary, by = "bin_id", all.x = TRUE, sort = FALSE)
  bin_table$gene_count[is.na(bin_table$gene_count)] <- 0L
  bin_table$gene_count <- as.integer(bin_table$gene_count)
  # gene_ids / gene_names stay NA on purpose for bins with no matching gene
  
  bin_table$genes <- bin_table$gene_names
  list(bin_table = bin_table, gene_hits = gene_hits)
}

# 3. Validation

# Cross-checks the annotation for internal consistency and prints a short summary of coverage on both sides
validate_gene_bin_annotation <- function(bin_table, gene_reference, gene_hits) {
  bin_chr_lookup <- setNames(as.character(bin_table$chr), bin_table$bin_id)
  bad_chr <- gene_hits$chr != bin_chr_lookup[gene_hits$bin_id]
  if (any(bad_chr)) {
    stop(sum(bad_chr), " gene-bin match(es) span different chromosomes - this should never happen.")
  }
  
  bad_overlap <- !(gene_hits$gene_start <= gene_hits$bin_end & gene_hits$gene_end >= gene_hits$bin_start)
  if (any(bad_overlap)) {
    stop(sum(bad_overlap), " gene-bin pair(s) violate the overlap rule ",
         "(gene_start <= bin_end AND gene_end >= bin_start).")
  }
  
  listed_counts <- ifelse(
    is.na(bin_table$gene_names) | !nzchar(bin_table$gene_names), 0L,
    lengths(strsplit(bin_table$gene_names, ";", fixed = TRUE))
  )
  mismatched <- bin_table$gene_count != listed_counts
  if (any(mismatched)) {
    stop(sum(mismatched), " bin(s) have gene_count that disagrees with the ",
         "number of genes actually listed in gene_names.")
  }
  
  n_bins <- nrow(bin_table)
  n_bins_with_gene <- sum(bin_table$gene_count > 0)
  n_ref_genes <- nrow(gene_reference)
  matched_genes <- unique(gene_hits$gene_name)
  unmatched_genes <- setdiff(gene_reference$gene_name, matched_genes)
  skipped <- attr(gene_reference, "skipped_genes")
  
  message("Gene-bin annotation validated: ", nrow(gene_hits),
          " overlapping pair(s) found, 0 rule violations.")
  message(n_bins_with_gene, " / ", n_bins, " bins have >=1 gene; ",
          length(matched_genes), " / ", n_ref_genes, " reference genes matched >=1 bin.")
  if (length(unmatched_genes) > 0) {
    message("Reference gene(s) with coordinates but no overlapping bin ",
            "(check bin coverage for that region): ", paste(unmatched_genes, collapse = ", "))
  }
  if (!is.null(skipped) && length(skipped$not_in_census) > 0) {
    message("Gene(s) requested but not found in the COSMIC census (skipped): ",
            paste(skipped$not_in_census, collapse = ", "))
  }
  if (!is.null(skipped) && length(skipped$missing_coordinates) > 0) {
    message("Gene(s) found in the census but missing coordinates (skipped): ",
            paste(skipped$missing_coordinates, collapse = ", "))
  }
  
  invisible(list(n_bins = n_bins, n_bins_with_gene = n_bins_with_gene,n_reference_genes = n_ref_genes, n_matched_genes = length(matched_genes), unmatched_genes = unmatched_genes, skipped_genes = skipped))
}