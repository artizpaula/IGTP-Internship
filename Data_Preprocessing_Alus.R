# Data_Preprocessing_Alus.R

suppressMessages({ library(data.table) })

raw_tsv_path<- "alu_methylation_matrix_all_runs_comb_norm.tsv"
data_app_rds_path <- "data_app.rds"
output_path <- "data_app_alus.rds"

chrom_list <- c(as.character(1:22), "X", "Y")

# 1: read the raw wide matrix and melt to long format
t0 <- Sys.time()
message("Reading raw Alu matrix (wide)...")
dt <- fread(raw_tsv_path, sep = "\t", header = TRUE, na.strings = "NA")
setnames(dt, "Alu", "alu_region")

# parse "chr:start-end" -> chr / alu_start / alu_end
m <- regmatches(dt$alu_region, regexec("^([^:]+):([0-9]+)-([0-9]+)$", dt$alu_region))
bad <- vapply(m, length, integer(1)) != 4
if (any(bad)) message("Warning: ", sum(bad), " alu_region values didn't parse as chr:start-end (dropped).")
dt <- dt[!bad]; m <- m[!bad]
dt[, chr := vapply(m, `[`, character(1), 2)]
dt[, alu_start := as.integer(vapply(m, `[`, character(1), 3))]
dt[, alu_end   := as.integer(vapply(m, `[`, character(1), 4))]
dt[, alu_id := paste(chr, alu_start, alu_end, sep = "_")]

sample_cols <- setdiff(names(dt), c("alu_region", "chr", "alu_start", "alu_end", "alu_id"))
message("Melting to long format (", nrow(dt), " Alus x ", length(sample_cols), " samples)...")
long <- melt(dt, id.vars = c("alu_id", "chr", "alu_start", "alu_end"), measure.vars = sample_cols,
             variable.name = "sample_id", value.name = "methylation", na.rm = TRUE)
long[, sample_id := tolower(as.character(sample_id))]
alu_coords <- dt[, .(alu_id, chr, alu_start, alu_end)]
rm(dt); gc()
message("Long rows (non-NA): ", nrow(long), "  t=", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s")

# 2: join patient metadata, build alu_table summary stats
metadata <- as.data.table(readRDS(data_app_rds_path)$metadata)
meta_small <- metadata[, .(sample_id, patient_id, Type)]

setkey(long, sample_id); setkey(meta_small, sample_id)
long[meta_small, on = "sample_id", `:=`(Type = i.Type, patient_id = i.patient_id)]
long <- long[!is.na(Type)]  # keep only samples that exist in metadata
long[, patient_id := as.character(patient_id)]
message("Rows after metadata join: ", nrow(long), "  t=", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s")

n_total_by_type <- meta_small[, .N, by = Type]
n_total_tumor  <- n_total_by_type[Type == "Tumor",  N]
n_total_normal <- n_total_by_type[Type == "Normal", N]

alu_table <- long[, .(
  mean_methylation_tumor  = mean(methylation[Type == "Tumor"],  na.rm = TRUE),
  mean_methylation_normal = mean(methylation[Type == "Normal"], na.rm = TRUE),
  sd_methylation_tumor    = sd(methylation[Type == "Tumor"],   na.rm = TRUE),
  sd_methylation_normal   = sd(methylation[Type == "Normal"],  na.rm = TRUE),
  n_present_Tumor  = sum(Type == "Tumor"),
  n_present_Normal = sum(Type == "Normal")
), by = .(alu_id, chr, alu_start, alu_end)]
alu_table[, tumor_normal_diff := mean_methylation_tumor - mean_methylation_normal]
alu_table[, n_total_Tumor  := n_total_tumor]
alu_table[, n_total_Normal := n_total_normal]
alu_table[, prevalence_tumor  := n_present_Tumor  / n_total_Tumor]
alu_table[, prevalence_normal := n_present_Normal / n_total_Normal]
alu_table[, n_samples_total := n_total_Tumor + n_total_Normal]
alu_table[, n_present_total := n_present_Tumor + n_present_Normal]
alu_table[, prevalence_total := n_present_total / n_samples_total]
alu_table[, alu_status := fifelse(n_present_total == 0L, "structural_gap",
                            fifelse(n_present_total == n_samples_total, "complete", "sample_specific_missing"))]

# gene/promoter annotation fallback (see header note) - overwrite these if you wire in the reference files
alu_table[, gene_count := 0L]; alu_table[, gene_ids := NA_character_]; alu_table[, gene_names := NA_character_]
alu_table[, promoter_count := 0L]; alu_table[, functional_annotation := NA_character_]

message("alu_table rows: ", nrow(alu_table), "  t=", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s")

# 3: per (patient, Alu) Tumor/Normal shift + per-Alu significance test
pat_alu_type <- long[, .(methylation = mean(methylation, na.rm = TRUE)), by = .(patient_id, alu_id, chr, Type)]
tumor_vals  <- pat_alu_type[Type == "Tumor",  .(patient_id, alu_id, chr, meth_tumor = methylation)]
normal_vals <- pat_alu_type[Type == "Normal", .(patient_id, alu_id, meth_normal = methylation)]
rm(pat_alu_type)
patient_alu_shift <- merge(tumor_vals, normal_vals, by = c("patient_id", "alu_id"))
patient_alu_shift[, shift := meth_tumor - meth_normal]
rm(tumor_vals, normal_vals)

alu_stats <- patient_alu_shift[, .(mean_shift = mean(shift, na.rm = TRUE), sd_shift = sd(shift, na.rm = TRUE), n = sum(is.finite(shift))), by = alu_id]
alu_stats[, t_stat := ifelse(n > 1 & sd_shift > 0, mean_shift / (sd_shift / sqrt(n)), NA_real_)]
alu_stats[, p_value := ifelse(!is.na(t_stat) & n > 2, 2 * pt(-abs(t_stat), df = pmax(n - 1, 1)), NA_real_)]
alu_stats[, q_value := p.adjust(p_value, method = "BH")]
alu_table <- merge(alu_table, alu_stats, by = "alu_id", all.x = TRUE)
rm(alu_stats); gc()
message("patient_alu_shift rows: ", nrow(patient_alu_shift), "  t=", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s")

# 4: shrink for memory (factors + drop redundant columns), order, save
alu_table[, chr := factor(chr, levels = chrom_list)]
setorder(alu_table, chr, alu_start)
alu_table <- as.data.frame(alu_table)

long[, alu_id := factor(alu_id, levels = alu_table$alu_id)]
long[, sample_id := factor(sample_id)]
methylation_long_alus <- as.data.frame(long[, .(sample_id, alu_id, methylation)]) # chr/alu_start/alu_end dropped, join from alu_table by alu_id if needed
rm(long); gc()

patient_alu_shift[, alu_id := factor(alu_id, levels = alu_table$alu_id)]
patient_alu_shift[, patient_id := factor(patient_id)]
patient_alu_shift[, chr := factor(chr, levels = chrom_list)]
patient_alu_shift <- as.data.frame(patient_alu_shift)

saveRDS(list(alu_table = alu_table, methylation_long_alus = methylation_long_alus, patient_alu_shift = patient_alu_shift),
        output_path)
message("Saved ", output_path, " (", round(file.size(output_path) / 1e6, 1), " MB). Total time: ",
        round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s")
