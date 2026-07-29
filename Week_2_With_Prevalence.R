# Week 2
# Structure/Calculate values that will be represented in the table
# (e.g. average methylation, tumor-normal difference, associated genes,
# functional annotations, prevalence data) prior to starting the application generation.

# Missing: alus and CpGs

library(jsonlite) # JSON reading

# Directory paths (personal paths)
dataset_metadata <- "/Users/paulaartizduenas/Desktop/Internship - IGTP/Dataset/Metadata"
dataset_bins <- "/Users/paulaartizduenas/Desktop/Internship - IGTP/Dataset/Bins"
metadata_dir <- file.path(dataset_metadata, "Metadata_all_runs_combined.csv")
bins_dir <- file.path(dataset_bins, "counts_bins_norm_mean")
output_dir  <- file.path(dirname(dataset_metadata), "Data Processed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Metadata
metadata <- read.delim(metadata_dir, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE) 
colnames(metadata) <- trimws(colnames(metadata)) 
if (ncol(metadata) <= 1) stop("Only 1 column detected - delimiter is wrong.")

mutation_cols <- c("KRAS", "BRAF", "TP53")
for (col in mutation_cols) {
  metadata[[col]] <- ifelse(is.na(metadata[[col]]), "unknown", as.character(metadata[[col]]))
} 

metadata$sample_id <- tolower(trimws(metadata$Sample))

id_extraction <- function(x) {
  out <- character(length(x))
  for (i in seq_along(x)) {
    val <- x[i]
    if (is.na(val)) {
      out[i] <- NA_character_
    } else {
      match_pos <- regexpr("[0-9]+", val)
      if (match_pos == -1) {
        out[i] <- NA_character_
      } else {
        out[i] <- regmatches(val, match_pos)
      }
    }
  }
  out
}
metadata$patient_id <- as.integer(id_extraction(metadata$Sample2))
write.csv(metadata, file.path(output_dir, "Metadata_clean.csv"), row.names = FALSE)

# 2. Loading and reshaping of bins
bin_files <- list.files(bins_dir, pattern = "^counts_.*\\.txt$", full.names = TRUE)

read_one_sample <- function(path) {
  file_name <- basename(path)
  sample_id <- sub("^counts_", "", file_name)
  sample_id <- sub("\\.txt$", "", sample_id)
  sample_id <- tolower(sample_id)
  
  bins <- fromJSON(path)
  chr_names <- names(bins)
  rows_per_chr <- vector("list", length(chr_names))
  
  for (i in seq_along(chr_names)) {
    chr <- chr_names[i]
    values <- bins[[chr]]
    bin_position <- as.integer(names(values))
    
    methylation <- numeric(length(values))
    for (j in seq_along(values)) {
      v <- values[[j]]
      methylation[j] <- if (is.null(v)) NA_real_ else as.numeric(v)
    }
    
    rows_per_chr[[i]] <- data.frame(
      sample_id = sample_id,
      chr = chr,
      bin_position = bin_position,
      methylation = methylation,
      stringsAsFactors = FALSE)
  }
  
  do.call(rbind, rows_per_chr)
}

samples_list <- lapply(bin_files, read_one_sample)
methylation_long <- do.call(rbind, samples_list)
write.csv(methylation_long, file.path(output_dir, "Methylation_long.csv"), row.names = FALSE)

# 3. Null values classifications
n_samples_per_bin <- aggregate(
  list(n_samples = rep(1L, nrow(methylation_long))),
  by = list(chr = methylation_long$chr, bin_position = methylation_long$bin_position),
  FUN = sum
)
n_na_per_bin <- aggregate(
  list(n_na = as.integer(is.na(methylation_long$methylation))),
  by = list(chr = methylation_long$chr, bin_position = methylation_long$bin_position),
  FUN = sum
)
bin_na_summary <- merge(n_samples_per_bin, n_na_per_bin, by = c("chr", "bin_position"))

bin_na_summary$bin_status <- ifelse(
  bin_na_summary$n_na == 0, "complete",
  ifelse(bin_na_summary$n_na == bin_na_summary$n_samples, "structural_gap", "sample_specific_missing")
)

# 4. Sample ID matching validation
samples_meta <- unique(metadata$sample_id)
samples_bins <- unique(methylation_long$sample_id)

only_in_meta <- setdiff(samples_meta, samples_bins)
only_in_bins <- setdiff(samples_bins, samples_meta)

if (length(only_in_meta) > 0) warning("Samples in metadata with no methylation file: ", paste(only_in_meta, collapse = ", "))
if (length(only_in_bins) > 0) warning("Methylation files with no metadata entry: ", paste(only_in_bins, collapse = ", "))

# 5. Bin annotation
bin_annotation <- unique(bin_na_summary[, c("chr", "bin_position")])
bin_annotation$bin_id <- paste(bin_annotation$chr, bin_annotation$bin_position, sep = "_")
bin_annotation$bin_start <- bin_annotation$bin_position - 999999L
bin_annotation$bin_end <- bin_annotation$bin_position
bin_annotation$n_cpg <- NA_integer_
bin_annotation$n_alu <- NA_integer_
bin_annotation$genes <- NA_character_ 
bin_annotation$functional_annotation <- NA_character_

bin_annotation <- bin_annotation[, c(
  "bin_id", "chr", "bin_position", "bin_start", "bin_end",
  "n_cpg", "n_alu", "genes", "functional_annotation"
)]

write.csv(bin_annotation, file.path(output_dir, "Bin_annotation_template.csv"), row.names = FALSE)

# 6. Tumor/Normal pairing and paired differences
merged <- merge(methylation_long, metadata[, c("sample_id", "Type", "patient_id")], by = "sample_id")

merged_small <- merged[, c("chr", "bin_position", "patient_id", "Type", "methylation")]

tumor_side <- merged_small[merged_small$Type == "Tumor", c("chr", "bin_position", "patient_id", "methylation")]
names(tumor_side)[names(tumor_side) == "methylation"] <- "Tumor"

normal_side <- merged_small[merged_small$Type == "Normal", c("chr", "bin_position", "patient_id", "methylation")]
names(normal_side)[names(normal_side) == "methylation"] <- "Normal"

paired_diff <- merge(tumor_side, normal_side, by = c("chr", "bin_position", "patient_id"), all = TRUE)
paired_diff$diff <- paired_diff$Tumor - paired_diff$Normal

# 7. Prevalence values
merged$chr <- as.character(merged$chr)
merged$bin_position <- as.integer(merged$bin_position)
merged$Type <- as.character(merged$Type)

n_present_long <- aggregate(
  list(n_present = as.integer(!is.na(merged$methylation))),
  by = list(chr = merged$chr, bin_position = merged$bin_position, Type = merged$Type),
  FUN = sum
)
n_total_long <- aggregate(
  list(n_total = rep(1L, nrow(merged))),
  by = list(chr = merged$chr, bin_position = merged$bin_position, Type = merged$Type),
  FUN = sum
)
presence_long <- merge(n_present_long, n_total_long, by = c("chr", "bin_position", "Type"))

tumor_presence <- presence_long[presence_long$Type == "Tumor", c("chr", "bin_position", "n_present", "n_total")]
names(tumor_presence)[3:4] <- c("n_present_Tumor", "n_total_Tumor")

normal_presence <- presence_long[presence_long$Type == "Normal", c("chr", "bin_position", "n_present", "n_total")]
names(normal_presence)[3:4] <- c("n_present_Normal", "n_total_Normal")

presence_summary <- merge(tumor_presence, normal_presence, by = c("chr", "bin_position"), all = TRUE)

expected_cols <- c("n_present_Tumor", "n_present_Normal", "n_total_Tumor", "n_total_Normal")
for (col in expected_cols) {
  if (!col %in% names(presence_summary)) presence_summary[[col]] <- 0L
  presence_summary[[col]][is.na(presence_summary[[col]])] <- 0L
}

prevalence_summary <- data.frame(
  chr = presence_summary$chr,
  bin_position = presence_summary$bin_position,
  n_present_Normal = presence_summary$n_present_Normal,
  n_present_Tumor = presence_summary$n_present_Tumor,
  n_total_Normal = presence_summary$n_total_Normal,
  n_total_Tumor = presence_summary$n_total_Tumor,
  stringsAsFactors = FALSE
)

prevalence_summary$prevalence_tumor <- ifelse(
  prevalence_summary$n_total_Tumor > 0,
  prevalence_summary$n_present_Tumor / prevalence_summary$n_total_Tumor,
  NA_real_
)
prevalence_summary$prevalence_normal <- ifelse(
  prevalence_summary$n_total_Normal > 0,
  prevalence_summary$n_present_Normal / prevalence_summary$n_total_Normal,
  NA_real_
)
prevalence_summary$n_samples_total <- prevalence_summary$n_total_Tumor + prevalence_summary$n_total_Normal
prevalence_summary$n_present_total <- prevalence_summary$n_present_Tumor + prevalence_summary$n_present_Normal
prevalence_summary$prevalence_total <- ifelse(
  prevalence_summary$n_samples_total > 0,
  prevalence_summary$n_present_total / prevalence_summary$n_samples_total,
  NA_real_
)

write.csv(prevalence_summary, file.path(output_dir, "Bin_prevalence_detection.csv"), row.names = FALSE)

# 8. Final Table

# 8.1 Mean and SD computation
mean_sd_summary <- aggregate(
  list(
    mean_methylation_tumor  = paired_diff$Tumor,
    mean_methylation_normal = paired_diff$Normal,
    sd_methylation_tumor    = paired_diff$Tumor,
    sd_methylation_normal   = paired_diff$Normal,
    tumor_normal_diff       = paired_diff$diff
  ),
  by = list(chr = paired_diff$chr, bin_position = paired_diff$bin_position),
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
)

bin_methylation_summary <- data.frame(
  chr                     = mean_sd_summary$chr,
  bin_position            = mean_sd_summary$bin_position,
  mean_methylation_tumor  = mean_sd_summary$mean_methylation_tumor[, "mean"],
  mean_methylation_normal = mean_sd_summary$mean_methylation_normal[, "mean"],
  sd_methylation_tumor    = mean_sd_summary$sd_methylation_tumor[, "sd"],
  sd_methylation_normal   = mean_sd_summary$sd_methylation_normal[, "sd"],
  tumor_normal_diff       = mean_sd_summary$tumor_normal_diff[, "mean"]
)

bin_methylation_summary$bin_id <- paste(bin_methylation_summary$chr, bin_methylation_summary$bin_position, sep = "_")

# 8.2 Merges with rest of info
bin_annotation$chr <- as.character(bin_annotation$chr)
bin_annotation$bin_position <- as.integer(bin_annotation$bin_position)
bin_na_summary$chr <- as.character(bin_na_summary$chr)
bin_na_summary$bin_position <- as.integer(bin_na_summary$bin_position)

bin_table <- merge(bin_annotation, bin_na_summary, by = c("chr", "bin_position"), all.x = TRUE, sort = FALSE)
bin_methylation_summary_no_keys <- bin_methylation_summary[, !(names(bin_methylation_summary) %in% c("chr", "bin_position"))]
bin_table <- merge(bin_table, bin_methylation_summary_no_keys, by = "bin_id", all.x = TRUE, sort = FALSE)
bin_table <- merge(bin_table, prevalence_summary, by = c("chr", "bin_position"), all.x = TRUE, sort = FALSE)

# 8.3 Sd included in final table
bin_table <- bin_table[, c(
  "bin_id", "chr", "bin_position", "bin_start", "bin_end", "bin_status",
  "mean_methylation_tumor", "mean_methylation_normal", 
  "sd_methylation_tumor", "sd_methylation_normal", "tumor_normal_diff",
  "n_cpg", "n_alu", "genes", "functional_annotation",
  "n_present_Normal", "n_present_Tumor", "n_total_Normal", "n_total_Tumor",
  "prevalence_tumor", "prevalence_normal",
  "n_samples_total", "n_present_total", "prevalence_total"
)]

write.csv(bin_table, file.path(output_dir, "Bin_table.csv"), row.names = FALSE)

# RDS with full table
saveRDS(
  list(metadata = metadata, methylation_long = methylation_long, bin_table = bin_table),
  file.path(output_dir, "data_app.rds")
)
