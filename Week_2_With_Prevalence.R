# Week 2
# Structure/Calculate values that will be represented in the table
# (e.g. average methylation, tumor-normal difference, associated genes,
# functional annotations, prevalence data) prior to starting the application generation.

# Missing: CpG Values, Alu Values, functional and gene annotation

library(jsonlite) # JSON reading

# Directory paths (personal paths)
dataset_metadata <- "/Users/paulaartizduenas/Desktop/Internship - IGTP/Dataset/Metadata"
dataset_bins <- "/Users/paulaartizduenas/Desktop/Internship - IGTP/Dataset/Bins"
metadata_dir <- file.path(dataset_metadata, "Metadata_all_runs_combined.csv")
bins_dir <- file.path(dataset_bins, "counts_bins_norm_mean")
output_dir  <- file.path(dirname(dataset_metadata), "Data Processed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Metadata

# 1.1 File delimiter: file is tab-separated, not comma-separated
metadata <- read.delim(metadata_dir, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE) # Read .csv, tab delimited 
colnames(metadata) <- trimws(colnames(metadata)) # Remove Leading/Trailing Whitespace
if (ncol(metadata) <= 1) stop("Only 1 column detected - delimiter is wrong.")

# 1.2 NA handling: NA = not tested, must stay "unknown", never defined as wild-type/0
mutation_cols <- c("KRAS", "BRAF", "TP53")
for (col in mutation_cols) {
  metadata[[col]] <- ifelse(is.na(metadata[[col]]), "unknown", as.character(metadata[[col]]))
} 

# normalize id for joining (disk files mix "run.." / "Run.." casing)
metadata$sample_id <- tolower(trimws(metadata$Sample))

# Loop -> extract the first numeric run of digits from Sample2 with a loop so that missing/odd values simply become NA instead of crashing (ex. Sample2 = N1 -> 1)
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

bin_files <- list.files(bins_dir, pattern = "^counts_.*\\.txt$", full.names = TRUE) # List all values with specified pattern: "^counts_.*\\.txt$" (regular expression)

read_one_sample <- function(path) {
  # recover sample id from file name: strip "counts_" prefix and ".txt" suffix
  file_name <- basename(path)
  sample_id <- sub("^counts_", "", file_name)
  sample_id <- sub("\\.txt$", "", sample_id)
  sample_id <- tolower(sample_id)
  
  bins <- fromJSON(path) # each file stores a nested json: chromosome -> bin -> value {chr: {bin_position: value}}
  
  # flatten the nested list into one row per bin, chromosome by chromosome
  chr_names <- names(bins)
  rows_per_chr <- vector("list", length(chr_names))
  
  for (i in seq_along(chr_names)) {
    chr <- chr_names[i]
    values <- bins[[chr]]
    bin_position <- as.integer(names(values))
    
    methylation <- numeric(length(values))
    for (j in seq_along(values)) {
      v <- values[[j]]
      # convert null json values into r missing values
      methylation[j] <- if (is.null(v)) NA_real_ else as.numeric(v)
    }
    
    rows_per_chr[[i]] <- data.frame(
      sample_id = sample_id,
      chr = chr,
      bin_position = bin_position,
      methylation = methylation,
      stringsAsFactors = FALSE)
  }
  
  do.call(rbind, rows_per_chr) # Saves data into a long format table
}

samples_list <- lapply(bin_files, read_one_sample)
methylation_long <- do.call(rbind, samples_list)
write.csv(methylation_long, file.path(output_dir, "Methylation_long.csv"), row.names = FALSE)

# 3. Null values classifications

# count, per bin, how many samples exist in total and how many of those are NA
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

# Distinguish biological absence from sample-specific missing values
# complete: no missing values (no NA), sample_specific_missing: missing only in some samples, structural_gap: missing in every sample
bin_na_summary$bin_status <- ifelse(
  bin_na_summary$n_na == 0, "complete",
  ifelse(bin_na_summary$n_na == bin_na_summary$n_samples, "structural_gap", "sample_specific_missing")
)

print(table(bin_na_summary$bin_status))

# 4. Sample ID matching validation

samples_meta <- unique(metadata$sample_id)
samples_bins <- unique(methylation_long$sample_id)

only_in_meta <- setdiff(samples_meta, samples_bins)
only_in_bins <- setdiff(samples_bins, samples_meta)

if (length(only_in_meta) > 0) warning("Samples in metadata with no methylation file: ", paste(only_in_meta, collapse = ", ")) # samples in metadata missing a methylation file
if (length(only_in_bins) > 0) warning("Methylation files with no metadata entry: ", paste(only_in_bins, collapse = ", ")) # methylation files with no metadata entry
if (length(only_in_meta) == 0 && length(only_in_bins) == 0) message("All sample IDs match between metadata and methylation files.") # all sample IDs match

# 5. Bin annotation

bin_annotation <- unique(bin_na_summary[, c("chr", "bin_position")])
bin_annotation$bin_id <- paste(bin_annotation$chr, bin_annotation$bin_position, sep = "_")
bin_annotation$bin_start <- bin_annotation$bin_position - 999999L # bins represent 1 mb genomic windows
bin_annotation$bin_end <- bin_annotation$bin_position
bin_annotation$n_cpg <- NA_integer_ # Information not given yet
bin_annotation$n_alu <- NA_integer_ # Information not given yet
bin_annotation$genes <- NA_character_ 
bin_annotation$functional_annotation <- NA_character_

bin_annotation <- bin_annotation[, c(
  "bin_id", "chr", "bin_position", "bin_start", "bin_end",
  "n_cpg", "n_alu", "genes", "functional_annotation"
)]

write.csv(bin_annotation, file.path(output_dir, "Bin_annotation_template.csv"), row.names = FALSE)

# 6. Tumor/Normal pairing and paired differences

merged <- merge(methylation_long, metadata[, c("sample_id", "Type", "patient_id")], by = "sample_id")

# check that every patient has exactly one Tumor and one Normal sample
type_per_patient <- unique(merged[, c("patient_id", "Type")])
pair_counts <- as.data.frame(table(type_per_patient$patient_id))
names(pair_counts) <- c("patient_id", "n")
pairing_check <- pair_counts[pair_counts$n != 2, ]
if (nrow(pairing_check) > 0) warning("Patients without a complete Tumor/Normal pair: ", paste(pairing_check$patient_id, collapse = ", "))

# spread Type (Tumor/Normal) into two columns "by hand": split then merge back together
merged_small <- merged[, c("chr", "bin_position", "patient_id", "Type", "methylation")]

tumor_side <- merged_small[merged_small$Type == "Tumor", c("chr", "bin_position", "patient_id", "methylation")]
names(tumor_side)[names(tumor_side) == "methylation"] <- "Tumor"

normal_side <- merged_small[merged_small$Type == "Normal", c("chr", "bin_position", "patient_id", "methylation")]
names(normal_side)[names(normal_side) == "methylation"] <- "Normal"

paired_diff <- merge(tumor_side, normal_side, by = c("chr", "bin_position", "patient_id"), all = TRUE)
paired_diff$diff <- paired_diff$Tumor - paired_diff$Normal # positive values indicate higher tumor methylation

# 7. Prevalence values
# prevalence = in how many samples (from Tumor and/or Normal) this bin has a
# non-NA methylation value, with respect to the total samples in this group.
# It is calculated separately for Tumor and Normal, and also combined.

merged$chr <- as.character(merged$chr)
merged$bin_position <- as.integer(merged$bin_position)
merged$Type <- as.character(merged$Type)

# count, per bin and per Type, how many samples have a value and how many exist in total
n_present_long <- aggregate(
  list(n_present = as.integer(!is.na(merged$methylation))), # count samples with available methylation values
  by = list(chr = merged$chr, bin_position = merged$bin_position, Type = merged$Type),
  FUN = sum
)
n_total_long <- aggregate(
  list(n_total = rep(1L, nrow(merged))),
  by = list(chr = merged$chr, bin_position = merged$bin_position, Type = merged$Type),
  FUN = sum
)
presence_long <- merge(n_present_long, n_total_long, by = c("chr", "bin_position", "Type"))

# spread Tumor/Normal into their own columns, same "split and merge" trick as before
tumor_presence <- presence_long[presence_long$Type == "Tumor", c("chr", "bin_position", "n_present", "n_total")]
names(tumor_presence)[3:4] <- c("n_present_Tumor", "n_total_Tumor")

normal_presence <- presence_long[presence_long$Type == "Normal", c("chr", "bin_position", "n_present", "n_total")]
names(normal_presence)[3:4] <- c("n_present_Normal", "n_total_Normal")

presence_summary <- merge(tumor_presence, normal_presence, by = c("chr", "bin_position"), all = TRUE)

# ensure expected columns exist even if one group is absent, and fill any gaps with 0
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
# fraction of Tumor samples with a non-NA methylation value at that bin.
prevalence_summary$prevalence_tumor <- ifelse(
  prevalence_summary$n_total_Tumor > 0,
  prevalence_summary$n_present_Tumor / prevalence_summary$n_total_Tumor,
  NA_real_
)
# fraction of Normal samples with a non-NA methylation value at that bin.
prevalence_summary$prevalence_normal <- ifelse(
  prevalence_summary$n_total_Normal > 0,
  prevalence_summary$n_present_Normal / prevalence_summary$n_total_Normal,
  NA_real_
)
prevalence_summary$n_samples_total <- prevalence_summary$n_total_Tumor + prevalence_summary$n_total_Normal
prevalence_summary$n_present_total <- prevalence_summary$n_present_Tumor + prevalence_summary$n_present_Normal
# fraction of all samples (Tumor + Normal combined) with a non-NA methylation value at that bin.
prevalence_summary$prevalence_total <- ifelse(
  prevalence_summary$n_samples_total > 0,
  prevalence_summary$n_present_total / prevalence_summary$n_samples_total,
  NA_real_
)

write.csv(prevalence_summary, file.path(output_dir, "Bin_prevalence_detection.csv"), row.names = FALSE)

# 8. Final Table

bin_methylation_summary <- aggregate(
  list(
    mean_methylation_tumor = paired_diff$Tumor,   # average methylation across paired samples
    mean_methylation_normal = paired_diff$Normal,
    tumor_normal_diff = paired_diff$diff
  ),
  by = list(chr = paired_diff$chr, bin_position = paired_diff$bin_position),
  FUN = function(x) mean(x, na.rm = TRUE)
)
bin_methylation_summary$bin_id <- paste(bin_methylation_summary$chr, bin_methylation_summary$bin_position, sep = "_")

bin_annotation$chr <- as.character(bin_annotation$chr)
bin_annotation$bin_position <- as.integer(bin_annotation$bin_position)
bin_na_summary$chr <- as.character(bin_na_summary$chr)
bin_na_summary$bin_position <- as.integer(bin_na_summary$bin_position)

bin_table <- merge(bin_annotation, bin_na_summary, by = c("chr", "bin_position"), all.x = TRUE, sort = FALSE)
bin_methylation_summary_no_keys <- bin_methylation_summary[, !(names(bin_methylation_summary) %in% c("chr", "bin_position"))]
bin_table <- merge(bin_table, bin_methylation_summary_no_keys, by = "bin_id", all.x = TRUE, sort = FALSE)
bin_table <- merge(bin_table, prevalence_summary, by = c("chr", "bin_position"), all.x = TRUE, sort = FALSE)

bin_table <- bin_table[, c(
  "bin_id", "chr", "bin_position", "bin_start", "bin_end", "bin_status",
  "mean_methylation_tumor", "mean_methylation_normal", "tumor_normal_diff",
  "n_cpg", "n_alu", "genes", "functional_annotation",
  "n_present_Normal", "n_present_Tumor", "n_total_Normal", "n_total_Tumor",
  "prevalence_tumor", "prevalence_normal",
  "n_samples_total", "n_present_total", "prevalence_total"
)]

write.csv(bin_table, file.path(output_dir, "Bin_table.csv"), row.names = FALSE)

# Save all precomputed objects for the shiny app
saveRDS(
  list(metadata = metadata, methylation_long = methylation_long, bin_table = bin_table),
  file.path(output_dir, "data_app.rds")
)
