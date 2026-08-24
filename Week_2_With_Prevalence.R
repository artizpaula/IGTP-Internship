# Week 2
# Structure/Calculate values that will be represented in the table (e.g. average methylation, tumor-normal difference, associated genes, functional annotations, prevalence data)

library(jsonlite)   # for reading the JSON files
library(data.table) # fast reading/aggregation for the big raw per-CpG methylation tables

setwd("/Users/paulaartizduenas/Desktop/Project/R Scripts")
# CpG/Alu-per-bin counting functions, used in section 5b below
source("Alus_and_CpGs.R")
# gene annotation functions (bin <-> gene coordinate overlap), used in section 5c below
source("Gene_Annotation.R")

# Directory paths (personal paths)
# Paths
dataset_metadata <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Metadata"
dataset_bins <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Bins"
dataset_alus_cpgs <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Alus_CpGs"
dataset_gene_annotation <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Archivo"

# Directories
metadata_dir <- file.path(dataset_metadata, "Metadata_all_runs_combined.csv")
bins_dir <- file.path(dataset_bins, "counts_bins_norm_mean")
output_dir  <- file.path(dirname(dataset_metadata), "Data Processed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Gene annotation source file (from archive.zip: COSMIC Cancer Gene Census)
cosmic_tsv_path    <- file.path(dataset_gene_annotation, "Cosmic_CancerGeneCensus_Tsv_v101_GRCh37", "Cosmic_CancerGeneCensus_v101_GRCh37.tsv")
crc_gene_list_path <- file.path(dataset_gene_annotation, "CRC_curated_genes.txt")

# 1. Metadata

# 1.1 File delimiter: it's tab-separated, not comma-separated
metadata <- read.delim(metadata_dir, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE) # read the csv, tab delimited 
colnames(metadata) <- trimws(colnames(metadata)) # remove leading/trailing whitespace from col names
if (ncol(metadata) <= 1) stop("Only 1 column detected - delimiter is wrong.")

# 1.2 NA handling: NA means not tested, so it has to stay "unknown", never treated as wild-type/0
mutation_cols <- c("KRAS", "BRAF", "TP53")
for (col in mutation_cols) {
  metadata[[col]] <- ifelse(is.na(metadata[[col]]), "unknown", as.character(metadata[[col]]))} 

# normalize id for joining later (files on disk mix "run.." / "Run.." casing)
metadata$sample_id <- tolower(trimws(metadata$Sample))

# loop to pull the first run of digits out of Sample2, so missing/weird values just become NA
# instead of crashing the whole thing (ex. Sample2 = N1 -> 1)
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

bin_files <- list.files(bins_dir, pattern = "^counts_.*\\.txt$", full.names = TRUE) # list all files matching the pattern "^counts_.*\\.txt$"

read_one_sample <- function(path) {
  # get sample id back from the file name: strip "counts_" prefix and ".txt" suffix
  file_name <- basename(path)
  sample_id <- sub("^counts_", "", file_name)
  sample_id <- sub("\\.txt$", "", sample_id)
  sample_id <- tolower(sample_id)
  
  bins <- fromJSON(path) # each file is a nested json: chromosome -> bin -> value {chr: {bin_position: value}}
  
  # flatten the nested list into one row per bin, going chromosome by chromosome
  chr_names <- names(bins)
  rows_per_chr <- vector("list", length(chr_names))
  
  for (i in seq_along(chr_names)) {
    chr <- chr_names[i]
    values <- bins[[chr]]
    bin_position <- as.integer(names(values))
    
    methylation <- numeric(length(values))
    for (j in seq_along(values)) {
      v <- values[[j]]
      # turn null json values into r's NA
      methylation[j] <- if (is.null(v)) NA_real_ else as.numeric(v)
    }
    
    rows_per_chr[[i]] <- data.frame(
      sample_id = sample_id,
      chr = chr,
      bin_position = bin_position,
      methylation = methylation,
      stringsAsFactors = FALSE)
  }
  
  do.call(rbind, rows_per_chr) # stack it all into one long-format table
}

samples_list <- lapply(bin_files, read_one_sample)
methylation_long <- do.call(rbind, samples_list)
write.csv(methylation_long, file.path(output_dir, "Methylation_long.csv"), row.names = FALSE)

# 3. Null values classifications

# count, per bin, how many samples exist in total and how many of those are NA
n_samples_per_bin <- aggregate(
  list(n_samples = rep(1L, nrow(methylation_long))),
  by = list(chr = methylation_long$chr, bin_position = methylation_long$bin_position),
  FUN= sum
)
n_na_per_bin <- aggregate(
  list(n_na = as.integer(is.na(methylation_long$methylation))),
  by = list(chr = methylation_long$chr, bin_position = methylation_long$bin_position),
  FUN= sum
)
bin_na_summary <- merge(n_samples_per_bin, n_na_per_bin, by = c("chr", "bin_position"))

# distinguish biological absence from sample-specific missing values:
# complete = no NAs at all, sample_specific_missing = missing only in some samples, structural_gap = missing in every sample
bin_na_summary$bin_status <- ifelse(bin_na_summary$n_na == 0, "complete",
  ifelse(bin_na_summary$n_na == bin_na_summary$n_samples, "structural_gap", "sample_specific_missing")
)

print(table(bin_na_summary$bin_status))

# 4. Sample ID matching validation

samples_meta <- unique(metadata$sample_id)
samples_bins <- unique(methylation_long$sample_id)

only_in_meta <- setdiff(samples_meta, samples_bins)
only_in_bins <- setdiff(samples_bins, samples_meta)

if (length(only_in_meta) > 0) warning("Samples in metadata with no methylation file: ", paste(only_in_meta, collapse = ", ")) # metadata samples missing a methylation file
if (length(only_in_bins) > 0) warning("Methylation files with no metadata entry: ", paste(only_in_bins, collapse = ", ")) # methylation files with no metadata entry
if (length(only_in_meta) == 0 && length(only_in_bins) == 0) message("All sample IDs match between metadata and methylation files.") # everything matches up

# 5. Bin annotation

bin_annotation <- unique(bin_na_summary[, c("chr", "bin_position")])
bin_annotation$bin_id <- paste(bin_annotation$chr, bin_annotation$bin_position, sep = "_")
bin_annotation$bin_start <- bin_annotation$bin_position - 999999L # bins are 1 mb genomic windows
bin_annotation$bin_end <- bin_annotation$bin_position
bin_annotation$genes <- NA_character_ 
bin_annotation$functional_annotation <- NA_character_

# 5b. CpG & Alu counts per 1 Mb bin

cpg_alu_dir <- file.path(output_dir, "CpG_Alu_per_sample")
raw_meth_files <- list_raw_methylation_files(dataset_alus_cpgs)

if (length(raw_meth_files$tar_files) == 0 && length(raw_meth_files$flat_files) == 0) {
  warning("No raw methylation files found in ", dataset_alus_cpgs)
} else {
  cpg_alu_samples <- process_all_raw_methylation(input_dir  = dataset_alus_cpgs,output_dir = cpg_alu_dir, bin_size   = 1e6)
  cpg_alu_annotation <- aggregate_bin_annotation(cpg_alu_samples, fun = "max")
  write.csv(cpg_alu_annotation, file.path(output_dir, "CpG_Alu_bin_annotation.csv"), row.names = FALSE)
  
  # fill in n_cpg / n_alu / %methylation by matching on chr + bin_start + bin_end
  key_bins <- paste(bin_annotation$chr, bin_annotation$bin_start, bin_annotation$bin_end)
  key_cpg_alu <- paste(cpg_alu_annotation$chr, cpg_alu_annotation$bin_start, cpg_alu_annotation$bin_end)
  idx <- match(key_bins, key_cpg_alu)
  bin_annotation$n_cpg <- cpg_alu_annotation$n_cpg[idx]
  bin_annotation$n_alu <- cpg_alu_annotation$n_alu[idx]
  # mean_meth = per-Alu average %methylation (each Alu counted once); pct_meth = % of
  # CpG-Alu observations that are methylated. see Alus_and_CpGs.R::summarize_bins()
  bin_annotation$alu_methylation_pct <- cpg_alu_annotation$mean_meth[idx]
  bin_annotation$cpg_methylation_pct <- cpg_alu_annotation$pct_meth[idx]
}

# guard columns so the script still runs (with NA %methylation) if no raw Alu/CpG files were found above
if (!"alu_methylation_pct" %in% names(bin_annotation)) bin_annotation$alu_methylation_pct <- NA_real_
if (!"cpg_methylation_pct" %in% names(bin_annotation)) bin_annotation$cpg_methylation_pct <- NA_real_

# 5c. Gene annotation (full COSMIC Cancer Gene Census, matched to bins by coordinate

if (!file.exists(cosmic_tsv_path)) {
  warning("Gene annotation source file not found (looked for ", cosmic_tsv_path,
          "); bin_annotation$gene_count/gene_ids/gene_names will be empty.")
  bin_annotation$gene_count <- 0L
  bin_annotation$gene_ids <- NA_character_
  bin_annotation$gene_names <- NA_character_
} else {
  gene_reference <- read_gene_reference(cosmic_tsv_path) # full census; all associated genes per bin, not just the curated CRC subset
  gene_annotation_result <- annotate_bins_with_genes(bin_annotation, gene_reference)
  bin_annotation <- gene_annotation_result$bin_table
  validate_gene_bin_annotation(bin_annotation, gene_reference, gene_annotation_result$gene_hits)
  write.csv(gene_annotation_result$gene_hits, file.path(output_dir, "Gene_bin_overlaps.csv"), row.names = FALSE)
}

bin_annotation <- bin_annotation[, c("bin_id", "chr", "bin_position", "bin_start", "bin_end", "n_cpg", "n_alu", "alu_methylation_pct", "cpg_methylation_pct", "gene_count", "gene_ids", "gene_names", "genes", "functional_annotation")]
write.csv(bin_annotation, file.path(output_dir, "Bin_annotation_template.csv"), row.names = FALSE)

# 6. Tumor/Normal pairing and paired differences

merged <- merge(methylation_long, metadata[, c("sample_id", "Type", "patient_id")], by = "sample_id")

# check that every patient has exactly one Tumor and one Normal sample
type_per_patient <- unique(merged[, c("patient_id", "Type")])
pair_counts <- as.data.frame(table(type_per_patient$patient_id))
names(pair_counts) <- c("patient_id", "n")
pairing_check <- pair_counts[pair_counts$n != 2, ]
if (nrow(pairing_check) > 0) warning("Patients without a complete Tumor/Normal pair: ", paste(pairing_check$patient_id, collapse = ", "))

# spread Type (Tumor/Normal) into two columns "by hand" - split then merge back together
merged_small <- merged[, c("chr", "bin_position", "patient_id", "Type", "methylation")]

tumor_side <- merged_small[merged_small$Type == "Tumor", c("chr", "bin_position", "patient_id", "methylation")]
names(tumor_side)[names(tumor_side) == "methylation"] <- "Tumor"

normal_side <- merged_small[merged_small$Type == "Normal", c("chr", "bin_position", "patient_id", "methylation")]
names(normal_side)[names(normal_side) == "methylation"] <- "Normal"

paired_diff <- merge(tumor_side, normal_side, by = c("chr", "bin_position", "patient_id"), all = TRUE)
paired_diff$diff <- paired_diff$Tumor - paired_diff$Normal # positive = higher tumor methylation

# 7. Prevalence values
# prevalence = in how many samples (Tumor and/or Normal) this bin has a non-NA methylation value, relative to the total samples in that group
merged$chr <- as.character(merged$chr)
merged$bin_position <- as.integer(merged$bin_position)
merged$Type <- as.character(merged$Type)

# count, per bin and per Type, how many samples have a value vs how many exist in total
n_present_long <- aggregate(
  list(n_present = as.integer(!is.na(merged$methylation))), # samples that actually have a methylation value
  by = list(chr = merged$chr, bin_position = merged$bin_position, Type = merged$Type),
  FUN= sum
)
n_total_long <- aggregate(
  list(n_total = rep(1L, nrow(merged))),
  by = list(chr = merged$chr, bin_position = merged$bin_position, Type = merged$Type),
  FUN= sum
)
presence_long <- merge(n_present_long, n_total_long, by = c("chr", "bin_position", "Type"))

# spread Tumor/Normal into their own columns, same split-and-merge trick as before
tumor_presence <- presence_long[presence_long$Type == "Tumor", c("chr", "bin_position", "n_present", "n_total")]
names(tumor_presence)[3:4] <- c("n_present_Tumor", "n_total_Tumor")

normal_presence <- presence_long[presence_long$Type == "Normal", c("chr", "bin_position", "n_present", "n_total")]
names(normal_presence)[3:4] <- c("n_present_Normal", "n_total_Normal")

presence_summary <- merge(tumor_presence, normal_presence, by = c("chr", "bin_position"), all = TRUE)

# make sure expected columns exist even if a whole group is missing, and fill any gaps with 0
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
  stringsAsFactors = FALSE)

# fraction of Tumor samples with a non-NA methylation value at that bin
prevalence_summary$prevalence_tumor <- ifelse(
  prevalence_summary$n_total_Tumor > 0,
  prevalence_summary$n_present_Tumor / prevalence_summary$n_total_Tumor,
  NA_real_
)
# fraction of Normal samples with a non-NA methylation value at that bin
prevalence_summary$prevalence_normal <- ifelse(
  prevalence_summary$n_total_Normal > 0,
  prevalence_summary$n_present_Normal / prevalence_summary$n_total_Normal,
  NA_real_
)
prevalence_summary$n_samples_total <- prevalence_summary$n_total_Tumor + prevalence_summary$n_total_Normal
prevalence_summary$n_present_total <- prevalence_summary$n_present_Tumor + prevalence_summary$n_present_Normal
# fraction of all samples (Tumor + Normal combined) with a non-NA methylation value at that bin
prevalence_summary$prevalence_total <- ifelse(
  prevalence_summary$n_samples_total > 0,
  prevalence_summary$n_present_total / prevalence_summary$n_samples_total,
  NA_real_
)

write.csv(prevalence_summary, file.path(output_dir, "Bin_prevalence_detection.csv"), row.names = FALSE)

# 8. Final Table

# 8.1 Mean and SD computation
mean_sd_summary <- aggregate(
  list(mean_methylation_tumor  = paired_diff$Tumor,
    mean_methylation_normal = paired_diff$Normal,
    sd_methylation_tumor = paired_diff$Tumor,
    sd_methylation_normal = paired_diff$Normal,
    tumor_normal_diff = paired_diff$diff
  ),
  by = list(chr = paired_diff$chr, bin_position = paired_diff$bin_position),
  FUN= function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
)

bin_methylation_summary <- data.frame(
  chr = mean_sd_summary$chr,
  bin_position = mean_sd_summary$bin_position,
  mean_methylation_tumor  = mean_sd_summary$mean_methylation_tumor[, "mean"],
  mean_methylation_normal = mean_sd_summary$mean_methylation_normal[, "mean"],
  sd_methylation_tumor = mean_sd_summary$sd_methylation_tumor[, "sd"],
  sd_methylation_normal = mean_sd_summary$sd_methylation_normal[, "sd"],
  tumor_normal_diff = mean_sd_summary$tumor_normal_diff[, "mean"]
)

bin_methylation_summary$bin_id <- paste(bin_methylation_summary$chr, bin_methylation_summary$bin_position, sep = "_")

# 8.2 Merge everything else back together
bin_annotation$chr <- as.character(bin_annotation$chr)
bin_annotation$bin_position <- as.integer(bin_annotation$bin_position)
bin_na_summary$chr <- as.character(bin_na_summary$chr)
bin_na_summary$bin_position <- as.integer(bin_na_summary$bin_position)

bin_table <- merge(bin_annotation, bin_na_summary, by = c("chr", "bin_position"), all.x = TRUE, sort = FALSE)
bin_methylation_summary_no_keys <- bin_methylation_summary[, !(names(bin_methylation_summary) %in% c("chr", "bin_position"))]
bin_table <- merge(bin_table, bin_methylation_summary_no_keys, by = "bin_id", all.x = TRUE, sort = FALSE)
bin_table <- merge(bin_table, prevalence_summary, by = c("chr", "bin_position"), all.x = TRUE, sort = FALSE)

# 8.3 SD included in the final table
bin_table <- bin_table[, c("bin_id", "chr", "bin_position", "bin_start", "bin_end", "bin_status",
  "mean_methylation_tumor", "mean_methylation_normal", 
  "sd_methylation_tumor", "sd_methylation_normal", "tumor_normal_diff",
  "n_cpg", "n_alu", "alu_methylation_pct", "cpg_methylation_pct",
  "gene_count", "gene_ids", "gene_names", "genes", "functional_annotation",
  "n_present_Normal", "n_present_Tumor", "n_total_Normal", "n_total_Tumor",
  "prevalence_tumor", "prevalence_normal",
  "n_samples_total", "n_present_total", "prevalence_total")]

write.csv(bin_table, file.path(output_dir, "Bin_table.csv"), row.names = FALSE)

# RDS with the full table, this is what app.R actually loads
saveRDS(list(metadata = metadata, methylation_long = methylation_long, bin_table = bin_table), file.path(output_dir, "data_app.rds"))