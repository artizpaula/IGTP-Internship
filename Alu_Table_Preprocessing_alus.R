# Alu-level preprocessing (associated with Data_Preprocessing.R)

library(data.table)

# Paths (change to personal pathways, same convention as Data_Preprocessing.R)
alu_matrix_path <- "alu_methylation_matrix_all_runs_comb_norm.tsv"
data_app_path   <- "data_app.rds"
output_path     <- "data_app_alus.rds"
cosmic_tsv_path    <- "Cosmic_CancerGeneCensus_v101_GRCh37.tsv"
promoter_tsv_path  <- "Promoter_reference_GRCh37.tsv"

source("Gene_Annotation.R")
source("Promoter_Annotation.R")

chrom_list <- c(as.character(1:22), "X", "Y")

# 1. Read the Alu matrix and reshape wide

message("Reading Alu methylation matrix...")
alu_wide <- fread(alu_matrix_path, sep = "\t", header = TRUE, na.strings = c("", "NA"))
setnames(alu_wide, 1, "alu_id")

sample_cols <- setdiff(names(alu_wide), "alu_id")

alu_long <- melt(alu_wide, id.vars = "alu_id", measure.vars = sample_cols,
                 variable.name = "sample_id", value.name = "methylation", na.rm = TRUE)
alu_long[, sample_id := tolower(as.character(sample_id))]
rm(alu_wide); gc()

# 2. Parse "chr:start-end" into chr/alu_start/alu_end

alu_coords <- unique(data.table(alu_id = alu_long$alu_id))
alu_coords[, c("chr", "range") := tstrsplit(alu_id, ":", fixed = TRUE)]
alu_coords[, c("alu_start", "alu_end") := tstrsplit(range, "-", fixed = TRUE)]
alu_coords[, range := NULL]
alu_coords[, alu_start := as.integer(alu_start)]
alu_coords[, alu_end   := as.integer(alu_end)]

# 3. Bring in Tumor/Normal + patient_id from metadata

data_app <- readRDS(data_app_path)
metadata <- data_app$metadata
bin_table <- data_app$bin_table

meta_small <- as.data.table(metadata[, c("sample_id", "Type", "patient_id")])
meta_small[, sample_id := tolower(sample_id)]

alu_long <- merge(alu_long, meta_small, by = "sample_id", all.x = TRUE)
if (any(is.na(alu_long$Type))) {
  warning(sum(is.na(alu_long$Type)), " methylation value(s) belong to a sample_id ",
          "not found in metadata and were dropped.")
  alu_long <- alu_long[!is.na(Type)]
}

n_total_Tumor  <- sum(metadata$Type == "Tumor")
n_total_Normal <- sum(metadata$Type == "Normal")

message("Summarizing per-Alu Tumor/Normal statistics...")
alu_stats_by_type <- alu_long[, .(
  n_present = .N,
  mean_meth = mean(methylation, na.rm = TRUE),
  sd_meth   = sd(methylation, na.rm = TRUE)
), by = .(alu_id, Type)]

wide_mean <- dcast(alu_stats_by_type, alu_id ~ Type, value.var = "mean_meth")
setnames(wide_mean, c("Tumor", "Normal"), c("mean_methylation_tumor", "mean_methylation_normal"))
wide_sd <- dcast(alu_stats_by_type, alu_id ~ Type, value.var = "sd_meth")
setnames(wide_sd, c("Tumor", "Normal"), c("sd_methylation_tumor", "sd_methylation_normal"))
wide_n <- dcast(alu_stats_by_type, alu_id ~ Type, value.var = "n_present")
setnames(wide_n, c("Tumor", "Normal"), c("n_present_Tumor", "n_present_Normal"))

alu_table <- Reduce(function(a, b) merge(a, b, by = "alu_id", all = TRUE),
                    list(alu_coords, wide_mean, wide_sd, wide_n))

for (col in c("n_present_Tumor", "n_present_Normal")) {
  alu_table[[col]][is.na(alu_table[[col]])] <- 0L
}
alu_table[, tumor_normal_diff := mean_methylation_tumor - mean_methylation_normal]
alu_table[, n_total_Tumor := n_total_Tumor]
alu_table[, n_total_Normal := n_total_Normal]
alu_table[, prevalence_tumor  := ifelse(n_total_Tumor  > 0, n_present_Tumor  / n_total_Tumor,  NA_real_)]
alu_table[, prevalence_normal := ifelse(n_total_Normal > 0, n_present_Normal / n_total_Normal, NA_real_)]
alu_table[, n_samples_total  := n_total_Tumor + n_total_Normal]
alu_table[, n_present_total  := n_present_Tumor + n_present_Normal]
alu_table[, prevalence_total := ifelse(n_samples_total > 0, n_present_total / n_samples_total, NA_real_)]

# 4. Assign each Alu to its parent 1 Mb bin

bin_size <- 1e6
alu_table[, bin_start := as.integer((alu_start - 1L) %/% bin_size * bin_size + 1L)]
alu_table[, bin_end   := as.integer(bin_start + bin_size - 1L)]
alu_table[, bin_id    := paste(chr, bin_end, sep = "_")]

# 5. Gene / promoter annotation

annotate_alus_directly <- file.exists(cosmic_tsv_path) || file.exists(promoter_tsv_path)

if (annotate_alus_directly) {
  alu_for_annot <- as.data.frame(alu_table[, .(bin_id = alu_id, chr, bin_start = alu_start, bin_end = alu_end)])

  if (file.exists(cosmic_tsv_path)) {
    gene_reference <- read_gene_reference(cosmic_tsv_path)
    gene_result <- annotate_bins_with_genes(alu_for_annot, gene_reference)
    alu_for_annot <- gene_result$bin_table
  } else {
    warning("Gene annotation source file not found (", cosmic_tsv_path, "); skipping direct Alu gene annotation.")
    alu_for_annot$gene_count <- 0L; alu_for_annot$gene_ids <- NA_character_; alu_for_annot$gene_names <- NA_character_
  }

  if (file.exists(promoter_tsv_path)) {
    promoter_reference <- read_promoter_reference(promoter_tsv_path)
    promoter_result <- annotate_bins_with_promoters(alu_for_annot, promoter_reference)
    alu_for_annot <- promoter_result$bin_table
  } else {
    warning("Promoter annotation source file not found (", promoter_tsv_path, "); skipping direct Alu promoter annotation.")
    alu_for_annot$promoter_count <- 0L
  }

  gene_promoter_by_alu <- setNames(
    alu_for_annot[, c("bin_id", "gene_count", "gene_ids", "gene_names", "promoter_count")],
    c("alu_id", "gene_count", "gene_ids", "gene_names", "promoter_count"))
  alu_table <- merge(alu_table, as.data.table(gene_promoter_by_alu), by = "alu_id", all.x = TRUE)

} else {
  message("COSMIC/promoter reference files not found next to this script; ",
          "inheriting gene/promoter annotation from each Alu's parent 1 Mb bin instead ",
          "(bin_table already has this precomputed in data_app.rds).")
  bin_annot <- as.data.table(bin_table[, c("bin_id", "gene_count", "gene_ids", "gene_names", "promoter_count")])
  alu_table <- merge(alu_table, bin_annot, by = "bin_id", all.x = TRUE)
  alu_table[, gene_count := ifelse(is.na(gene_count), 0L, gene_count)]
  alu_table[, promoter_count := ifelse(is.na(promoter_count), 0L, promoter_count)]
}

# 6. Final column order/tidy-up (mirrors bin_table's column order)

alu_table <- alu_table[, .(alu_id, chr, alu_start, alu_end, bin_id,
                          mean_methylation_tumor, mean_methylation_normal,
                          sd_methylation_tumor, sd_methylation_normal, tumor_normal_diff,
                          gene_count, gene_ids, gene_names, promoter_count,
                          n_present_Normal, n_present_Tumor, n_total_Normal, n_total_Tumor,
                          prevalence_tumor, prevalence_normal,
                          n_samples_total, n_present_total, prevalence_total)]

alu_table <- as.data.frame(alu_table)
alu_table$chr <- factor(alu_table$chr, levels = chrom_list)
alu_table <- alu_table[order(alu_table$chr, alu_table$alu_start), ]
alu_table$chr <- as.character(alu_table$chr)
rownames(alu_table) <- NULL

# human-readable "start\u2013end" label
add_thousands_sep <- function(x) gsub("\\B(?=(\\d{3})+(?!\\d))", ",", format(x, scientific = FALSE, trim = TRUE), perl = TRUE)
alu_table$alu_coordinates <- paste0(add_thousands_sep(alu_table$alu_start), "\u2013", add_thousands_sep(alu_table$alu_end))

message("Built alu_table with ", format(nrow(alu_table), big.mark = ","), " Alu elements.")

saveRDS(list(alu_table = alu_table), output_path)
message("Saved ", output_path)
