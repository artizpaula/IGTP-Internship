# Exploration of Epigenomic Data in Colorectal Cancer
# Paula Artiz Dueñas - UPC

library(shiny)
library(bslib)
library(bsicons) # icons for value boxes
library(DT)
library(plotly)
library(ggplot2)
library(patchwork) # aligns the annotation strip + heatmap with independent legends

data <- readRDS("/Users/paulaartizduenas/Desktop/Project/Dataset/Data Processed/data_app.rds")

metadata <- data$metadata
bin_table <- data$bin_table
methylation_long <- data$methylation_long
chrom_list <- c(as.character(1:22), "X", "Y")

# bin choices for the sidebar dropdown, grouped by chromosome and sorted by position
bin_table <- bin_table[order(match(bin_table$chr, chrom_list), bin_table$bin_position), ]
bin_choices_by_chr <- split(bin_table$bin_id, factor(bin_table$chr, levels = chrom_list))

# Bin Table tab: extra columns we need + a couple of lookup tables

# turns a bin's start/end into something readable, e.g. "1 MB-2 MB"
format_mb_label <- function(bp) {
  mb <- bp/1e6
  if (abs(mb - round(mb)) < 1e-9) paste0(round(mb)," MB")
  else paste0(format(round(mb, 1), nsmall = 1)," MB")
}
format_bin_coordinates <- function(start, end) {
  paste0(format_mb_label(start), "\u2013", format_mb_label(end))
}
bin_table$bin_coordinates <- mapply(format_bin_coordinates, bin_table$bin_start, bin_table$bin_end)

# SD wasn't saved in data_app.rds so just fill it with NA if it's missing
if (!"sd_methylation_tumor" %in% names(bin_table)) bin_table$sd_methylation_tumor <- NA_real_
if (!"sd_methylation_normal" %in% names(bin_table)) bin_table$sd_methylation_normal <- NA_real_

# tumor - normal mean difference
bin_table$mean_diff_tumor_normal <- bin_table$mean_methylation_tumor - bin_table$mean_methylation_normal

# Alu/CpG counts per bin, already computed upstream in Week_2_With_Prevalence.R
if (!"n_alu" %in% names(bin_table)) bin_table$n_alu <- NA_real_
if (!"n_cpg" %in% names(bin_table)) bin_table$n_cpg <- NA_real_

# gene annotation columns (gene_count / gene_ids / gene_names)
gene_annotation_present <- all(c("gene_count", "gene_ids", "gene_names") %in% names(bin_table)) &&
  any(!is.na(bin_table$gene_names))

if (!gene_annotation_present) {
  
  find_first_existing <- function(paths) {
    hit <- paths[file.exists(paths)]
    if (length(hit) > 0) hit[[1]] else NA_character_
  }
  # candidate locations for the two source files the annotation needs, tried in order
  # (matches the "Archivo" folder layout used by Week_2_With_Prevalence.R)
  cosmic_tsv_path <- find_first_existing(c(
    "Archivo/Cosmic_CancerGeneCensus_Tsv_v101_GRCh37/Cosmic_CancerGeneCensus_v101_GRCh37.tsv",
    "../Archivo/Cosmic_CancerGeneCensus_Tsv_v101_GRCh37/Cosmic_CancerGeneCensus_v101_GRCh37.tsv",
    "Data/Archivo/Cosmic_CancerGeneCensus_Tsv_v101_GRCh37/Cosmic_CancerGeneCensus_v101_GRCh37.tsv"
  ))
  crc_gene_list_path <- find_first_existing(c(
    "Archivo/CRC_curated_genes.txt",
    "../Archivo/CRC_curated_genes.txt",
    "Data/Archivo/CRC_curated_genes.txt"
  ))
  gene_annotation_script <- find_first_existing(c(
    "Gene_and_Functional_Annotation.R",
    "../Gene_and_Functional_Annotation.R"
  ))
  
  if (!is.na(cosmic_tsv_path) && !is.na(crc_gene_list_path) && !is.na(gene_annotation_script)) {
    source(gene_annotation_script, local = TRUE)
    gene_reference <- read_gene_reference(cosmic_tsv_path, crc_gene_list_path)
    gene_annotation_result <- annotate_bins_with_genes(bin_table, gene_reference)
    bin_table <- gene_annotation_result$bin_table
    # re-apply the same chromosome/position ordering used above, in case merge() reshuffled rows
    bin_table <- bin_table[order(match(bin_table$chr, chrom_list), bin_table$bin_position), ]
    message("Bin Table: gene annotation columns were missing from data_app.rds and have been ",
            "recomputed on load from the COSMIC Cancer Gene Census + curated CRC gene list.")
  } else {
    warning("Bin Table: gene annotation source files not found (looked for the COSMIC census TSV, ",
            "the curated CRC gene list, and Gene_and_Functional_Annotation.R next to app.R or one ",
            "level up). Gene Count / Genes / Gene IDs will show as empty until data_app.rds is ",
            "regenerated with gene annotation, or these files are made available alongside app.R.")
    if (!"gene_count" %in% names(bin_table)) bin_table$gene_count <- 0L
    if (!"gene_ids" %in% names(bin_table)) bin_table$gene_ids <- NA_character_
    if (!"gene_names" %in% names(bin_table)) bin_table$gene_names <- NA_character_
    if (!"genes" %in% names(bin_table)) bin_table$genes <- bin_table$gene_names
  }
}

# whether any bin actually ended up with gene annotation
gene_annotation_available <- any(!is.na(bin_table$gene_names))

# sorted list of every distinct gene symbol annotated on any bin
all_gene_names <- sort(unique(trimws(unlist(strsplit(na.omit(bin_table$gene_names), ";", fixed = TRUE)))))

# Alu/CpG methylation percentage columns, computed upstream in Week_2_With_Prevalence.Rre)
if (!"alu_methylation_pct" %in% names(bin_table)) {
  bin_table$alu_methylation_pct <- if ("mean_meth" %in% names(bin_table)) bin_table$mean_meth else NA_real_
}
if (!"cpg_methylation_pct" %in% names(bin_table)) {
  bin_table$cpg_methylation_pct <- if ("pct_meth" %in% names(bin_table)) bin_table$pct_meth else NA_real_
}

# columns to show in the Bin Table tab, in the order we want them
bin_table_columns <- list(
  list(id = "bin_id",                 label = "Bin ID",              group = "Bin"),
  list(id = "bin_coordinates",        label = "Coordinates",         group = "Bin"),
  list(id = "mean_methylation_tumor", label = "Tumor Mean",          group = "Methylation (\u03b2-value)", digits = 3),
  list(id = "mean_methylation_normal",label = "Normal Mean",         group = "Methylation (\u03b2-value)", digits = 3),
  list(id = "sd_methylation_tumor",   label = "Tumor SD",            group = "Methylation (\u03b2-value)", digits = 3),
  list(id = "sd_methylation_normal",  label = "Normal SD",           group = "Methylation (\u03b2-value)", digits = 3),
  list(id = "mean_diff_tumor_normal", label = "\u0394 (Tumor \u2212 Normal)", group = "Methylation (\u03b2-value)", digits = 3),
  list(id = "n_alu",                  label = "Alu Count",           group = "Genomic Content"),
  list(id = "n_cpg",                  label = "CpG Count",           group = "Genomic Content"),
  list(id = "alu_methylation_pct",    label = "Alu Meth %",          group = "Genomic Content", digits = 2),
  list(id = "cpg_methylation_pct",    label = "CpG Meth %",          group = "Genomic Content", digits = 2),
  list(id = "gene_count",             label = "Gene Count",          group = "Gene Annotation (COSMIC CRC census)"),
  list(id = "gene_names",             label = "Genes",               group = "Gene Annotation (COSMIC CRC census)"),
  list(id = "gene_ids",               label = "Gene IDs (COSMIC)",   group = "Gene Annotation (COSMIC CRC census)"))

# two-row header for the Bin Table DT
bin_table_header_sketch <- local({
  col_groups <- vapply(bin_table_columns, function(c) c$group, character(1))
  group_rle  <- rle(col_groups)
  htmltools::withTags(table(
    class = "display",
    thead(
      tr(
        lapply(seq_along(group_rle$lengths), function(i) {
          th(
            colspan = group_rle$lengths[i],
            style = "text-align:center; background:#eef2f5; color:#16324f; border-bottom:2px solid #d8e0e6; font-size:11.5px; text-transform:uppercase; letter-spacing:0.4px;",
            group_rle$values[i]
          )
        })
      ),
      tr(
        lapply(bin_table_columns, function(c) th(style = "white-space:nowrap;", c$label))
      )
    )
  ))
})

# numeric metrics we expose as filters
bin_table_filters <- list(
  list(id = "mean_methylation_tumor", label = "Mean Meth (Tumor)", op = ">=", step = 0.01),
  list(id = "mean_methylation_normal", label = "Mean Meth (Normal)", op = ">=", step = 0.01),
  list(id = "sd_methylation_tumor", label = "SD Meth (Tumor)", op = ">=", step = 0.01),
  list(id = "sd_methylation_normal",label = "SD Meth (Normal)", op = ">=", step = 0.01),
  list(id = "mean_diff_tumor_normal", label = "Mean Difference", op = ">=", step = 0.01),
  list(id = "n_cpg", label = "CpG Count", op = ">=", step = 1),
  list(id = "n_alu", label = "Alu Count", op = ">=", step = 1),
  list(id = "alu_methylation_pct", label = "Alu Methylation %", op = ">=", step = 0.1),
  list(id = "cpg_methylation_pct", label = "CpG Methylation %", op = ">=", step = 0.1),
  list(id = "gene_count", label = "Gene Count", op = ">=", step = 1)
)

# builds the grid of slider inputs for the Bin Table filters
build_bin_table_filter_inputs <- function(filters, df) {
  lapply(filters, function(f) {
    vals <- df[[f$id]]
    vals <- vals[is.finite(vals)]
    step <- if (!is.null(f$step)) f$step else 0.01
    rng <- if (length(vals) > 0) range(vals) else c(0, 1)
    
    if (step >= 1) {
      lo <- floor(rng[1])
      hi <- ceiling(rng[2])
    } else {
      lo <- floor(rng[1] * 100) / 100
      hi <- ceiling(rng[2] * 100) / 100
    }
    if (hi <= lo) hi <- lo + step
    
    div(
      style = "flex: 1 1 200px; min-width: 180px; max-width: 250px;",
      sliderInput(
        inputId = paste0("binfilter_", f$id),
        label = paste0(f$label, " ", f$op),
        min = lo, max = hi, value = lo, step = step,
        width = "100%"
      )
    )
  })
}

# applies all the configured filters to a bin_table subset
apply_bin_table_filters <- function(df, filters, input) {
  # chromosome filter
  chr_selected <- input[["binfilter_chr"]]
  if (!is.null(chr_selected)) {
    df <- df[as.character(df$chr) %in% chr_selected, , drop = FALSE]
  }
  
  for (f in filters) {
    threshold <- input[[paste0("binfilter_", f$id)]]
    if (is.null(threshold)) next
    col <- df[[f$id]]
    if (is.null(col) || all(is.na(col))) next
    keep <- switch(f$op,
                   ">=" = !is.na(col) & col >= threshold,
                   ">"  = !is.na(col) & col >  threshold,
                   rep(TRUE, nrow(df))
    )
    df <- df[keep, , drop = FALSE]
  }
  
  # gene search: exact match
  gene_query <- input[["binfilter_gene"]]
  if (!is.null(gene_query) && nzchar(trimws(gene_query)) && "gene_names" %in% names(df)) {
    pattern <- paste0("(^|;)\\s*", toupper(trimws(gene_query)), "\\s*($|;)")
    keep <- !is.na(df$gene_names) & grepl(pattern, toupper(df$gene_names))
    df <- df[keep, , drop = FALSE]
  }
  
  df
}


# builds the data frame that actually gets displayed, based on bin_table_columns
build_display_bin_table <- function(df, columns = bin_table_columns, plain = FALSE) {
  col_ids <- vapply(columns, function(c) c$id, character(1))
  out <- df[, col_ids, drop = FALSE]
  for (i in seq_along(columns)) {
    if (!is.null(columns[[i]]$digits)) {
      out[[i]] <- round(out[[i]], columns[[i]]$digits)
    } else if (is.character(out[[i]])) {
      out[[i]][is.na(out[[i]])] <- "\u2014"
    }
  }
  names(out) <- vapply(columns, function(c) c$label, character(1))
  out
}

# length of each chromosome
chrom_lengths <- tapply(bin_table$bin_end, bin_table$chr, max)
chrom_lengths <- setNames(as.numeric(chrom_lengths[chrom_list]), chrom_list)

# running offset so we can lay all chromosomes out on one axis
chr_offset <- setNames(numeric(length(chrom_list)), chrom_list)
running <- 0
for (chr in chrom_list) {
  chr_offset[chr] <- running
  running <- running + chrom_lengths[[chr]]
}
chr_mid <- chr_offset + chrom_lengths / 2

# genome-wide bin table, used for the manhattan-style plots

genome_wide_bins <- bin_table
genome_wide_bins$chr <- factor(genome_wide_bins$chr, levels = chrom_list)
genome_wide_bins <- genome_wide_bins[order(genome_wide_bins$chr, genome_wide_bins$bin_position), ]
genome_wide_bins$genome_x <- chr_offset[as.character(genome_wide_bins$chr)] + genome_wide_bins$bin_position
genome_wide_bins$overall_meth <- rowMeans(
  genome_wide_bins[, c("mean_methylation_tumor", "mean_methylation_normal")], na.rm = TRUE)
genome_wide_bins$chr_parity <- factor(as.integer(genome_wide_bins$chr) %% 2)

# dashed lines marking where each chromosome starts, reused on both plots
chr_boundary_shapes <- lapply(as.numeric(chr_offset)[-1], function(x) {
  list(type = "line", x0 = x, x1 = x, y0 = 0, y1 = 1, yref = "paper",
       line = list(color = "rgba(150,150,150,0.4)", width = 1, dash = "dot"))
})

# Feature x Bin Heatmap: getting the per-patient data ready
bin_lookup <- bin_table[, c("chr", "bin_position", "bin_id")]
bin_lookup$chr <- as.character(bin_lookup$chr)
methylation_long$chr <- as.character(methylation_long$chr)
ml_annot <- merge(methylation_long, bin_lookup, by = c("chr", "bin_position"))
ml_annot <- merge(ml_annot, metadata[, c("sample_id", "patient_id", "Type")], by = "sample_id")
patient_bin_type <- aggregate(methylation ~ patient_id + bin_id + chr + Type, data = ml_annot, FUN = mean, na.rm = TRUE)

tumor_vals <- patient_bin_type[patient_bin_type$Type == "Tumor",  c("patient_id", "bin_id", "chr", "methylation")]
normal_vals <- patient_bin_type[patient_bin_type$Type == "Normal", c("patient_id", "bin_id", "methylation")]
names(tumor_vals)[4] <- "meth_tumor"
names(normal_vals)[3] <- "meth_normal"
patient_bin_shift <- merge(tumor_vals, normal_vals, by = c("patient_id", "bin_id"))
patient_bin_shift$shift <- patient_bin_shift$meth_tumor - patient_bin_shift$meth_normal
patient_bin_shift$bin_id <- factor(patient_bin_shift$bin_id, levels = bin_table$bin_id)

# genome-wide Tumor vs Normal significance testing 
bin_stats <- aggregate(
  shift ~ bin_id, data = patient_bin_shift,
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = sum(is.finite(x)))
)
bin_stats <- do.call(data.frame, bin_stats)
names(bin_stats) <- c("bin_id", "mean_shift", "sd_shift", "n")
bin_stats$bin_id <- as.character(bin_stats$bin_id)
bin_stats$t_stat <- ifelse(bin_stats$n > 1 & bin_stats$sd_shift > 0,
                           bin_stats$mean_shift / (bin_stats$sd_shift / sqrt(bin_stats$n)), NA_real_)
bin_stats$p_value <- ifelse(!is.na(bin_stats$t_stat) & bin_stats$n > 2,
                            2 * stats::pt(-abs(bin_stats$t_stat), df = pmax(bin_stats$n - 1, 1)), NA_real_)
bin_stats$q_value <- stats::p.adjust(bin_stats$p_value, method = "BH")

# one row per patient with the clinical fields we can use to annotate the heatmap
# (Recaiguda/BRAF/KRAS/TP53/MSS/sexe/estadi2 are the same for a patient's Tumor and Normal sample)
patient_annotation <- unique(metadata[metadata$Type == "Tumor", c(
  "patient_id", "Recaiguda", "BRAF", "KRAS", "TP53", "MSS", "sexe", "estadi2"
)])
patient_annotation$patient_id <- as.character(patient_annotation$patient_id)
patient_bin_shift$patient_id <- as.character(patient_bin_shift$patient_id)

# merge the per-bin significance stats computed above into the genome-wide bin table
genome_wide_bins <- merge(
  genome_wide_bins, bin_stats[, c("bin_id", "mean_shift", "sd_shift", "n", "p_value", "q_value")],
  by = "bin_id", all.x = TRUE
)
genome_wide_bins$chr <- factor(genome_wide_bins$chr, levels = chrom_list)
genome_wide_bins <- genome_wide_bins[order(genome_wide_bins$chr, genome_wide_bins$bin_position), ]

manhattan_sig_threshold <- 0.05
manhattan_outlier_threshold <- stats::quantile(abs(genome_wide_bins$mean_diff_tumor_normal), 0.95, na.rm = TRUE)

genome_wide_bins$manhattan_category <- with(genome_wide_bins, {
  sig <- !is.na(q_value) & q_value < manhattan_sig_threshold
  outlier <- !is.na(mean_diff_tumor_normal) & abs(mean_diff_tumor_normal) >= manhattan_outlier_threshold
  ifelse(sig & mean_diff_tumor_normal > 0, "Significant hypermethylation",
         ifelse(sig & mean_diff_tumor_normal < 0, "Significant hypomethylation",
                ifelse(outlier, "Extreme outlier", "Not significant")))
})
genome_wide_bins$manhattan_category <- factor(
  genome_wide_bins$manhattan_category,
  levels = c("Not significant", "Extreme outlier", "Significant hypomethylation", "Significant hypermethylation")
)

# options for the "which feature" dropdown on the heatmap tab
heatmap_feature_choices <- c(
  "Methylation Values (No Association)" = "none",
  "Relapse"  = "Recaiguda",
  "BRAF status" = "BRAF",
  "KRAS status" = "KRAS",
  "TP53 status" = "TP53",
  "MSI/MSS status"= "MSS",
  "Sex" = "sexe",
  "Stage" = "estadi2"
)

# turns a 0/1 (or other two-level) coded column into readable labels
recode_binary <- function(x, labels = c("0" = "Wild-type", "1" = "Mutant")) {
  x_chr <- as.character(x)
  out <- ifelse(!is.na(x_chr) & x_chr %in% names(labels), labels[x_chr], x_chr)
  out
}

# short p-value label used on plot annotations across the app
format_pvalue <- function(p) {
  if (is.null(p) || !is.finite(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", format(round(p, 3), nsmall = 3))
}

# runs a Wilcoxon rank-sum test or Welch's t-test between two numeric vectors
run_group_test <- function(x, y, method = c("wilcox", "ttest")) {
  method <- match.arg(method)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) {
    return(list(p_value = NA_real_, label = "Not enough samples for a statistical test",
                n_x = length(x), n_y = length(y), method = method))
  }
  res <- tryCatch(
    if (identical(method, "ttest")) stats::t.test(x, y) else stats::wilcox.test(x, y),
    error = function(e) NULL
  )
  if (is.null(res)) {
    return(list(p_value = NA_real_, label = "Statistical test unavailable",
                n_x = length(x), n_y = length(y), method = method))
  }
  method_label <- if (identical(method, "ttest")) "Welch's t-test" else "Wilcoxon rank-sum"
  list(p_value = res$p.value, label = paste0(method_label, ": ", format_pvalue(res$p.value)),
       n_x = length(x), n_y = length(y), method = method)
}

# Clinical Explorer: one-row-per-patient clinical profile table
patient_level_cols <- local({
  candidate_cols <- setdiff(names(metadata), c("sample_id", "Type"))
  keep <- vapply(candidate_cols, function(cn) {
    all(tapply(metadata[[cn]], metadata$patient_id, function(v) length(unique(v[!is.na(v)])) <= 1))
  }, logical(1))
  candidate_cols[keep]
})
if (!"patient_id" %in% patient_level_cols) patient_level_cols <- c("patient_id", patient_level_cols)

clinical_profile_table <- unique(metadata[, patient_level_cols, drop = FALSE])
clinical_profile_table$patient_id <- as.character(clinical_profile_table$patient_id)
clinical_profile_table <- clinical_profile_table[!duplicated(clinical_profile_table$patient_id), ]

# recode the fields that use a 0/1 or Catalan-language convention into readable labels
if ("sexe" %in% names(clinical_profile_table))
  clinical_profile_table$sexe <- recode_binary(clinical_profile_table$sexe, c("Dona" = "Female", "Home" = "Male"))
for (gene_col in intersect(c("BRAF", "KRAS", "TP53"), names(clinical_profile_table)))
  clinical_profile_table[[gene_col]] <- recode_binary(clinical_profile_table[[gene_col]])
if ("Recaiguda" %in% names(clinical_profile_table))
  clinical_profile_table$Recaiguda <- recode_binary(clinical_profile_table$Recaiguda, c("0" = "No", "1" = "Yes"))

# per-patient mean methylation (Tumor/Normal, across all samples/bins) and sample counts
sample_mean_meth <- aggregate(methylation ~ sample_id, data = methylation_long, FUN = mean, na.rm = TRUE)
sample_mean_meth <- merge(sample_mean_meth, metadata[, c("sample_id", "patient_id", "Type")], by = "sample_id")
sample_mean_meth$patient_id <- as.character(sample_mean_meth$patient_id)

patient_type_mean <- aggregate(methylation ~ patient_id + Type, data = sample_mean_meth, FUN = mean, na.rm = TRUE)
patient_type_mean <- reshape(patient_type_mean, idvar = "patient_id", timevar = "Type", direction = "wide")
names(patient_type_mean) <- sub("^methylation\\.", "", names(patient_type_mean))
if (!"Tumor"  %in% names(patient_type_mean)) patient_type_mean$Tumor  <- NA_real_
if (!"Normal" %in% names(patient_type_mean)) patient_type_mean$Normal <- NA_real_
names(patient_type_mean)[names(patient_type_mean) == "Tumor"]  <- "MeanMethTumor"
names(patient_type_mean)[names(patient_type_mean) == "Normal"] <- "MeanMethNormal"
patient_type_mean$MethShift <- patient_type_mean$MeanMethTumor - patient_type_mean$MeanMethNormal

sample_counts <- as.data.frame.matrix(table(metadata$patient_id, metadata$Type))
sample_counts$patient_id <- as.character(rownames(sample_counts))
names(sample_counts)[names(sample_counts) == "Tumor"]  <- "nTumorSamples"
names(sample_counts)[names(sample_counts) == "Normal"] <- "nNormalSamples"

clinical_profile_table <- merge(clinical_profile_table, patient_type_mean[, c("patient_id", "MeanMethTumor", "MeanMethNormal", "MethShift")], by = "patient_id", all.x = TRUE)
clinical_profile_table <- merge(clinical_profile_table, sample_counts[, intersect(c("patient_id", "nTumorSamples", "nNormalSamples"), names(sample_counts))], by = "patient_id", all.x = TRUE)
clinical_profile_table <- clinical_profile_table[order(as.numeric(suppressWarnings(clinical_profile_table$patient_id)), clinical_profile_table$patient_id), ]

# nicer display names for the clinical table's columns; anything not listed here
clinical_column_labels <- c(
  patient_id = "Patient ID", sexe = "Sex", estadi2 = "Stage", MSS = "MSI/MSS Status",
  BRAF = "BRAF", KRAS = "KRAS", TP53 = "TP53", Recaiguda = "Relapse",
  MeanMethTumor = "Mean Methylation (Tumor)", MeanMethNormal = "Mean Methylation (Normal)",
  MethShift = "Methylation Shift (Tumor \u2212 Normal)",
  nTumorSamples = "Tumor Samples", nNormalSamples = "Normal Samples"
)
clinical_label_for <- function(col_id) {
  if (col_id %in% names(clinical_column_labels)) return(unname(clinical_column_labels[col_id]))
  tools::toTitleCase(gsub("_", " ", col_id))
}

# builds a patient x bin matrix of methylation shift, only used to order rows via clustering when no clinical feature is picked
build_shift_matrix <- function(df, bin_ids) {
  m <- tapply(df$shift, list(df$patient_id, as.character(df$bin_id)), FUN = identity)
  keep_cols <- bin_ids[bin_ids %in% colnames(m)]
  m[, keep_cols, drop = FALSE]
}

# reshapes a long (id, key, value) data frame into a numeric id x key matrix
build_wide_matrix <- function(df, id_col, key_col, value_col) {
  df <- df[is.finite(df[[value_col]]), c(id_col, key_col, value_col)]
  if (nrow(df) == 0) return(NULL)
  wide <- reshape(df, idvar = id_col, timevar = key_col, direction = "wide")
  row_ids <- as.character(wide[[id_col]])
  wide[[id_col]] <- NULL
  names(wide) <- sub(paste0("^", value_col, "\\."), "", names(wide))
  mat <- as.matrix(wide)
  rownames(mat) <- row_ids
  mat <- mat[, colSums(!is.finite(mat)) == 0, drop = FALSE]
  mat <- mat[rowSums(!is.finite(mat)) == 0, , drop = FALSE]
  if (nrow(mat) == 0 || ncol(mat) == 0) return(NULL)
  mat
}

# categorical colour palette that adapts to however many levels there are, reused wherever a clinical field maps to colour
categorical_palette <- function(values) {
  base_palette <- c("#0e7c86", "#d1495b", "#3aa9c9", "#e0a339", "#2fae66", "#16324f", "#7d92a3")
  levels_now <- sort(unique(values))
  setNames(base_palette[((seq_along(levels_now) - 1) %% length(base_palette)) + 1], levels_now)
}

# shared ggplot2 theme so every static plot in the app looks consistent
theme_app <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.subtitle = element_text(size = rel(0.82), color = "#7d92a3", margin = margin(b = 10)),
      plot.caption= element_text(size = rel(0.68), color = "#7d92a3", hjust = 0),
      axis.title = element_text(size = rel(0.88), color = "#16324f"),
      axis.text = element_text(size = rel(0.82), color = "#3a5470"),
      legend.title= element_text(size = rel(0.85), color = "#16324f", face = "bold"),
      legend.text= element_text(size = rel(0.82), color = "#3a5470"),
      strip.text  = element_text(size = rel(0.88), face = "bold", color = "#16324f"),
      strip.background = element_rect(fill = "#eef2f5", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#e3e8ec", linewidth = 0.4),
      panel.spacing = unit(14, "pt"),
      plot.margin = margin(10, 14, 6, 6)
    )
}

# orders category labels sensibly instead of just sorting alphabetically
natural_level_order <- function(x) {
  vals <- unique(as.character(x))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  num <- suppressWarnings(as.numeric(vals))
  if (length(vals) > 0 && !anyNA(num)) return(vals[order(num)])
  roman_rank <- c("I" = 1, "II" = 2, "III" = 3, "IV" = 4, "V" = 5)
  root <- toupper(sub("^([IVX]+).*$", "\\1", vals))
  if (length(vals) > 0 && all(root %in% names(roman_rank))) {
    suffix <- toupper(sub("^[IVX]+", "", vals))
    return(vals[order(roman_rank[root], suffix)])
  }
  sort(vals)
}

# static lookup table of hg19/GRCh37 coordinates for colorectal-cancer driver genes.
# gene set expanded (82 genes) using https://www.intogen.org/search and UCSC Genome Browser 

gene_lookup_table <- data.frame(
  gene = c(
    "APC", "TP53", "KRAS", "PIK3CA", "SMAD4", "FBXW7", "SOX9", "AMER1",
    "BRAF", "TCF7L2", "LRP1B", "ATM", "ARID1A", "NRAS", "KMT2C", "GRIN2A",
    "PCBP1", "ERBB2", "MTOR", "CTNNB1", "BCL9L", "SMAD2", "CYLD", "ERBB3",
    "PTEN", "BCL9", "TAF1L", "SMAD3", "RNF43", "ERBB4", "PIK3R1", "ACVR2A",
    "MAP2K4", "ROBO2", "ARID1B", "TBX3", "ANK1", "PRKD1", "EPHA7", "SMARCA4",
    "TGIF1", "DUSP16", "NF1", "SIRPA", "RBM10", "ACVR1B", "MYH9", "EGFR",
    "BCOR", "NIN", "USP6", "EP300", "PTPN11", "AKT1", "PBRM1", "ELF3",
    "PARP4", "MAP3K1", "TSC1", "ING1", "GNAS", "AFDN", "RRN3", "CSMD3",
    "NBEA", "FAT4", "FAT3", "CDKN2A", "HRAS", "RNF6", "NCOR2", "SLC34A2",
    "BIRC6", "RSPO2", "TGFBR2", "CTNNA1", "CEP89", "PPP2R1A", "DICER1", "CARD11",
    "DCC", "BCORL1"),
  chr   = c("5", "17", "12", "3", "18", "4", "17", "X",
            "7", "10", "2", "11", "1", "1", "7", "16",
            "2", "17", "1", "3", "11", "18", "16", "12",
            "10", "1", "9", "15", "17", "2", "5", "2",
            "17", "3", "6", "12", "8", "14", "6", "19",
            "18", "12", "17", "20", "X", "12", "22", "7",
            "X", "14", "17", "22", "12", "14", "3", "1",
            "13", "5", "9", "13", "20", "6", "16", "8",
            "13", "4", "11", "9", "11", "13", "12", "4",
            "2", "8", "3", "5", "19", "19", "14", "7",
            "18", "X"),
  start = c(112073582, 7571739, 25358180, 178866145, 48556583, 153241696, 70117161, 63404997,
            140419137, 114710006, 140988992, 108093794,27022506, 115247090, 151832010, 9847265,
            70314609, 37856317, 11166592, 41240996, 118766845, 45335328, 50775997, 56473949,
            89623382, 147013271, 32629452, 67357940, 56431037, 212240442, 67511584, 148602598,
            11924194, 77089250, 157098512, 115108060, 41510744, 30045685, 93949738, 11071706,
            3450170, 12626216, 29421995, 1876053, 47004620, 52345483, 36677326, 55086710, 
            39910504, 51186481, 5019327, 41488596, 112856751, 105235686, 52579383, 201979715,
            24995069, 56111376, 135766736, 111366047, 57414803, 168227591, 15153890, 113235157,
            35516407,  126236110, 91957984, 21967751, 532242, 26786912, 124808961, 25657473,
            32582091, 108911544,  30648093, 138089114, 33366831, 52693305, 95552565, 2945776,
            49866567, 129116611),
  end   = c(112181936, 7590808, 25403863, 178957881, 48611412, 153457244, 70122557, 63425588,
            140624729, 114927437, 142888585, 108239829, 27108595, 115259392, 152133088, 10276285,
            70316335, 37884911, 11322608, 41281934, 118796635, 45457030, 50835846, 56497289,
            89731687, 147098017, 32635667, 67487507, 56494895, 213403526, 67597649, 148688391,
            12047145, 77699115, 157531913, 115121980, 41655140, 30397053, 94129277, 11172949,
            3459976, 12715797, 29704693, 1921238, 47046212, 52390862, 36784012, 55279321,
            39957211, 51297880, 5078286, 41576081, 112947722, 105262085, 52719929, 201986311,
            25086916, 56191979, 135820003, 111375686, 57486247, 168372703, 15188129, 114449168,
            36246873, 126414087, 92629639, 21994391, 535576, 26796451, 125052158, 25680370,
            32843945, 109095848, 30735634, 138270723, 33462864, 52732771, 95623866, 3083501,
            51062269, 129192046),
  stringsAsFactors = FALSE
)

# parses whatever the user typed into the genome browser search box into a (chr, start, end) target
parse_genomic_search <- function(query, bin_table, chrom_list, gene_lookup_table, pad = 2e6) {
  q <- trimws(query)
  if (!nzchar(q)) return(list(status = "error", message = "Please enter a coordinate, bin ID, or gene name."))
  q_nospace <- gsub(",", "", q)
  
  # case 1: chr:start[-end] coordinate format
  coord_pattern <- "^(chr|Chr|CHR)?([0-9]{1,2}|[XxYy])\\s*:\\s*([0-9]+)\\s*(-\\s*([0-9]+))?$"
  if (grepl(coord_pattern, q_nospace)) {
    m <- regmatches(q_nospace, regexec(coord_pattern, q_nospace))[[1]]
    chr <- toupper(m[3])
    if (!(chr %in% chrom_list)) return(list(status = "error", message = paste0("Unknown chromosome '", chr, "'.")))
    start <- as.numeric(m[4])
    end <- if (nzchar(m[6])) as.numeric(m[6]) else start
    if (end < start) { tmp <- start; start <- end; end <- tmp }
    return(list(
      status = "ok", chr = chr, start = start, end = end,
      message = paste0("Jumped to chr", chr, ":", format(start, big.mark = ",", scientific = FALSE),
                       "-", format(end, big.mark = ",", scientific = FALSE))
    ))
  }
  
  # case 2: bin ID format "chr_position"
  bin_pattern <- "^([0-9]{1,2}|[XxYy])_([0-9]+)$"
  if (grepl(bin_pattern, q_nospace)) {
    m <- regmatches(q_nospace, regexec(bin_pattern, q_nospace))[[1]]
    chr <- toupper(m[2])
    pos <- as.numeric(m[3])
    if (chr %in% chrom_list) {
      hit <- bin_table[bin_table$chr == chr & bin_table$bin_position == pos, ]
      if (nrow(hit) == 1) {
        return(list(status = "ok", chr = chr, start = hit$bin_start, end = hit$bin_end,
                    message = paste0("Jumped to bin ", hit$bin_id)))
      }
      return(list(status = "ok", chr = chr, start = max(1, pos - 999999), end = pos,
                  message = paste0("No bin exactly at ", q, "; showing the nearest 1 Mb window.")))
    }
  }
  
  # case 3: gene symbol already annotated on one or more bins
  gene_query <- toupper(q)
  has_annotation <- !is.na(bin_table$genes) &
    grepl(paste0("(^|[,;])\\s*", gene_query, "\\s*($|[,;])"), toupper(bin_table$genes))
  in_annotation <- bin_table[has_annotation,]
  if (nrow(in_annotation) >= 1) {
    return(list(status = "ok", chr = in_annotation$chr[1],
                start = min(in_annotation$bin_start), end = max(in_annotation$bin_end),
                message = paste0("Jumped to gene ", gene_query)))
  }
  
  # case 4: gene symbol not in bin_table, so fall back to the static gene_lookup_table
  in_lookup <- gene_lookup_table[toupper(gene_lookup_table$gene) == gene_query, ]
  if (nrow(in_lookup) == 1) {
    return(list(status = "ok", chr = in_lookup$chr,
                start = max(1, in_lookup$start-pad), end = in_lookup$end + pad,
                message = paste0("Jumped to ", gene_query, " (", in_lookup$chr, ":",
                                 format(in_lookup$start, big.mark = ",", scientific = FALSE), "-",
                                 format(in_lookup$end, big.mark = ",", scientific = FALSE), ")")
    ))
  }
  list(status = "error", message = paste0("'", q, "' was not recognized as a coordinate, bin ID, or gene name."))
}

# parses an uploaded bin file (.csv, .tsv, .txt, or anything else) into a vector of bin_id values that exist in bin_table

parse_bin_file <- function(filepath, filename, bin_table, chrom_list) {
  ext <- tolower(tools::file_ext(filename))
  raw_ids <- character(0)
  tbl <- NULL
  
  # only try reading it as a proper table if it's csv/tsv
  if (ext %in% c("csv", "tsv")) {
    sep <- if (ext == "tsv") "\t" else ","
    tbl <- tryCatch(
      utils::read.csv(filepath, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL)
  }
  
  if (!is.null(tbl) && ncol(tbl) >= 1 && nrow(tbl) >= 1) {
    nm <- tolower(trimws(names(tbl)))
    id_col  <- which(nm %in% c("bin_id", "bin", "id"))
    chr_col <- which(nm %in% c("chr", "chromosome", "chrom"))
    pos_col <- which(nm %in% c("bin_position", "position", "pos", "start"))
    if (length(id_col) >= 1) {
      raw_ids <- as.character(tbl[[id_col[1]]])
    } else if (length(chr_col) >= 1 && length(pos_col) >= 1) {
      raw_ids <- paste0(tbl[[chr_col[1]]], "_", tbl[[pos_col[1]]])
    } else if (ncol(tbl) == 1) {
      # single unlabeled column - could be a header row that's actually the first ID
      raw_ids <- c(names(tbl)[1], as.character(tbl[[1]]))
    } else {
      raw_ids <- as.character(unlist(tbl))
    }
  }
  
  # otherwise just treat the file as plain unstructured text
  if (length(raw_ids) == 0) {
    txt <- tryCatch(readLines(filepath, warn = FALSE), error = function(e) character(0))
    txt <- paste(txt, collapse = "\n")
    raw_ids <- strsplit(txt, "[,;\\s]+", perl = TRUE)[[1]]
  }
  
  raw_ids <- trimws(raw_ids)
  raw_ids <- raw_ids[nzchar(raw_ids)]
  n_entries <- length(unique(raw_ids))
  
  # normalize formatting differences: drop the "chr" prefix, turn ":"/"-" separators into "_"
  norm <- toupper(raw_ids)
  norm <- sub("^CHR", "", norm)
  norm <- gsub("[:\\-]", "_", norm)
  norm <- gsub("_+", "_", norm)
  norm <- unique(norm)
  
  bin_lookup <- toupper(bin_table$bin_id)
  hit <- match(norm, bin_lookup)
  matched <- unique(bin_table$bin_id[hit[!is.na(hit)]])
  unmatched <- norm[is.na(hit)]
  
  list(matched = matched, unmatched = unmatched, n_entries = n_entries)
}

# zooms in and out on the genome browser
zoom_view <- function(start, end, chr_len, factor, min_width = 2e6) {
  width <- end - start
  center <- (start + end)/2
  new_width <- max(min_width, min(chr_len, width*factor)) # factor < 1 zooms in, > 1 zooms out
  new_start <- center - new_width/2
  new_end <- center + new_width/2
  if (new_start < 1) {new_end <- new_end + (1 - new_start); new_start <- 1}
  if (new_end > chr_len) {new_start <- new_start - (new_end - chr_len); new_end <- chr_len}
  new_start <- max(1, new_start)
  list(start = round(new_start), end = round(new_end))}

# UI Definition

ui <- page_sidebar(
  title = div(
    style = "
      display:flex;
      justify-content:space-between;
      align-items:center;
      width:100%;
      padding:18px 10px;
      border-bottom:1px solid rgba(255,255,255,0.10);
    ",
    
    # Left side
    div(
      h2(
        "Exploration of Epigenomic Data in Colorectal Cancer",
        style = "
          margin:0;
          font-weight:600;
          font-size:28px;
          line-height:1.2;
        "
      ),
      
      tags$div(
        "Institut Germans Trias i Pujol (IGTP) · Universitat Politècnica de Catalunya (UPC)",
        style = "
          font-size:13px;
          color:#9fb3c8;
          margin-top:4px;
        "
      )
    ),
    
    # Right side
    div(
      style = "display:flex; gap:14px; align-items:center; margin-left:auto;",
      
      div(
        style = "display:flex; align-items:center; gap:8px; border:1px solid #3a5470; border-radius:8px; padding:8px 14px;",
        bs_icon("person", size = "2em"),
        textOutput("n_patients", inline = TRUE),
        " patients"
      ),
      
      div(
        style = "display:flex; align-items:center; gap:8px; border:1px solid #3a5470; border-radius:8px; padding:8px 14px;",
        bs_icon("clipboard-data", size = "2em"),
        textOutput("n_samples", inline = TRUE),
        " samples"
      )
    )
  ),
  theme = bs_theme(
    version = 5,
    base_font = font_google("IBM Plex Sans"),
    heading_font = font_google("Libre Franklin"),
    bg = "#f4f7f9", fg = "#0b2436",
    primary = "#0e7c86", secondary = "#16324f",
    success  = "#2fae66", info = "#3aa9c9", warning = "#e0a339", danger = "#d1495b",
    "navbar-bg"  = "#0b2436",
    base_font_size_scale = 0.98
  ),
  
  sidebar = sidebar(
    width = 300,
    tags$div("DATA SELECTION", style = "font-size:14px; font-weight:700; letter-spacing:0.8px; color:#7d92a3; margin-bottom:0px;"),
    selectInput("chr", "Chromosome:", choices = chrom_list),
    selectizeInput("bins", "Selected bins:", choices = NULL, multiple = TRUE, options = list(placeholder = "Find and select bins (ex: 1_1000000)...", plugins = list("remove_button"))),
    actionButton("add_chr_bins", "Add all bins on the chromosome", icon = bs_icon("plus-circle"), class = "btn-sm class=btn-outline-light w-100"),
    actionButton("clear_bins", "Clear bin selection", icon = bs_icon("x-circle"), class = "btn-sm class=btn-outline-light w-100"),
    fileInput("bins_file", "Or upload bins to select:", accept = c(".csv", ".tsv", ".txt", "text/csv", "text/tab-separated-values", "text/plain"), placeholder = "No file selected", buttonLabel = "Browse..."),
    tags$div("Accepts .csv/.tsv/.txt or a plain list of bin IDs like 1_1000000, one per line or comma-separated).", style = "font-size:11px; color:#7d92a3; margin-top:-10px; margin-bottom:8px;"),
    hr(style = "margin:1px 0; border-color:#3a5470;")),
  
  navset_card_underline(
    title = "Exploration and Visualization",
    
    # 1. Overview
    nav_panel("Overview", uiOutput("overview_content")),
    
    # 2. Genome Browser
    nav_panel("Genome Browser", div(style = "font-size:13px; color:#6c757d; margin-bottom:10px;", "Browse methylation across the whole genome: pick a chromosome, zoom in/out, or search by coordinates, bin ID, or gene name."), card(card_header("All chromosomes, click a point to jump to that chromosome"), plotlyOutput("genome_overview_plot", height = "420px")), card(card_header(textOutput("browser_position_header", inline = TRUE)), div(style = "display:flex; flex-wrap:wrap; gap:0px; align-items:center", div(style = "display:flex; gap:6px; align-items:center; flex:1 1 320px; min-width:280px;", textInput("browser_search", label = NULL, placeholder = "e.g. 12:25000000-26000000, 12_25000000, or KRAS", width = "100%"), actionButton("browser_search_go", NULL, icon = bs_icon("search"), class = "btn-sm btn-outline-secondary")), div(style = "display:flex; gap:4px; align-items:center;", actionButton("prev_chr", NULL, icon = bs_icon("chevron-left"), class = "btn-sm btn-outline-secondary"), selectInput("browser_chr", NULL, choices = chrom_list, selected = "1", width = "90px"), actionButton("next_chr", NULL, icon = bs_icon("chevron-right"), class = "btn-sm btn-outline-secondary")), div(style = "display:flex; gap:4px;", actionButton("zoom_in", NULL, icon = bs_icon("zoom-in"), class = "btn-sm btn-outline-secondary"), actionButton("zoom_out", NULL, icon = bs_icon("zoom-out"), class = "btn-sm btn-outline-secondary"), actionButton("zoom_reset", "Whole chromosome", icon = bs_icon("arrow-counterclockwise"), class = "btn-sm btn-outline-secondary"))), plotlyOutput("browser_chr_plot", height = "560px"), tags$div("Tip: bins currently in your sidebar selection are outlined on the plot above. You can also drag directly on the plot to zoom, and double-click it to reset.", style = "font-size:11px; color:#7d92a3; margin-top:6px;"))),
    
    # 3. Tumor vs Normal
    nav_panel("Tumor vs. Normal", layout_columns(col_widths = c(6, 6, 12), card(card_header("Methylation density (Tumor vs Normal)"), plotOutput("tn_density")), card(card_header("PCA / UMAP (interactive)"), radioButtons("proj_method", NULL, choices = c("PCA", "UMAP"), inline = TRUE), plotlyOutput("tn_projection", height = "400px")), card(card_header("Patient similarity network"), plotOutput("tn_network", height = "500px")))),
    
    # 4. Genome-wide Profile
    nav_panel(
      "Genome-wide Profile",
      div(
        style = "font-size:13px; color:#6c757d; margin-bottom:10px;",
        "Genome-wide comparison of Tumor and Normal methylation across all 1 Mb bins. Each point is a bin; a paired test (per patient, Tumor vs. Normal) flags bins that differ significantly (Benjamini-Hochberg adjusted q < 0.05), and bins in the top 5% by absolute difference are flagged as extreme outliers even when not significant. Hover a point for the bin, chromosome, difference, and q-value; drag to zoom, double-click to reset, and use the camera icon to export."
      ),
      card(
        card_header("Manhattan-style plot: Tumor \u2212 Normal methylation difference"),
        plotlyOutput("manhattan_plot", height = "550px"),
        div(
          style = "display:flex; gap:20px; flex-wrap:wrap; margin-top:12px; padding-top:10px; border-top:1px solid #e3e8ec; font-size:12px; color:#3a5470;",
          tags$span(tags$span(style = "display:inline-block; width:10px; height:10px; border-radius:50%; background:#d1495b; margin-right:6px;"), "Significant hypermethylation (q < 0.05)"),
          tags$span(tags$span(style = "display:inline-block; width:10px; height:10px; border-radius:50%; background:#3aa9c9; margin-right:6px;"), "Significant hypomethylation (q < 0.05)"),
          tags$span(tags$span(style = "display:inline-block; width:10px; height:10px; border-radius:50%; background:#e0a339; margin-right:6px;"), "Extreme outlier (top 5% |\u0394|, not significant)"),
          tags$span(tags$span(style = "display:inline-block; width:10px; height:10px; border-radius:50%; background:#b7c1c9; margin-right:6px;"), "Not significant")
        )
      )
    ),
    
    # 5. Feature x Chromosome Heatmap
    nav_panel("Feature × Chromosome Heatmap", card(card_header("Diverging heatmap of methylation shift (feature × chromosome)"), div(style = "font-size:13px; color:#6c757d; margin-bottom:10px;", "Each row is a patient, each column a chromosome, and the colour is the Tumor \u2212 Normal methylation shift for that patient/chromosome. Select a feature below to add a colour strip showing that clinical variable alongside the heatmap."),
                                                   div(style = "max-width:500px; margin-bottom:10px;", selectInput("heatmap_feature", "Select Feature:", choices = heatmap_feature_choices, selected = "none")),
                                                   plotOutput("feature_heatmap", height = "1100px"))),
    
    # 6. Clinical Explorer
    nav_panel(
      "Clinical Explorer",
      div(
        style = "font-size:13px; color:#6c757d; margin-bottom:10px;",
        "Compare tumor methylation by mutation status, and cross-reference patients against their full clinical profile. Click a row in the table to highlight that patient in the plot; use the search box or column filters to narrow the table down."
      ),
      uiOutput("selected_patient_banner"),
      layout_columns(
        col_widths = c(5, 7),
        card(
          card_header("Methylation by mutation status"),
          div(
            style = "display:flex; flex-wrap:wrap; gap:18px; align-items:flex-end; margin-bottom:6px;",
            div(style = "min-width:160px;", selectInput("mutation_gene", "Gene:", choices = c("KRAS", "BRAF", "TP53"))),
            div(style = "min-width:220px;", radioButtons("mutation_stat_test", "Statistical test:", choices = c("Wilcoxon rank-sum" = "wilcox", "Welch's t-test" = "ttest"), selected = "wilcox", inline = TRUE))
          ),
          plotOutput("clinical_boxplot", height = "420px")
        ),
        card(
          card_header("Clinical metadata explorer"),
          div(style = "font-size:12px; color:#7d92a3; margin-bottom:8px;", "One row per patient. Use the column filters below the headers to narrow results, click any column header to sort, and click a row to select that patient."),
          DTOutput("clinical_table")
        )
      )
    ),
    
    # 7. Bin Table
    nav_panel(
      "Bin Table",
      card(
        card_header(
          div(
            style = "display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:8px;",
            div(style = "display:flex; align-items:center; gap:8px;", bs_icon("table", size = "1.05em"), "Filterable Bin-level Data"),
              uiOutput("bintable_count_badge")
            )
          ),
        
        tags$details(
          open = "open",
          style = "background:#ffffff; border:1px solid #d8e0e6; border-radius:10px; padding:0; margin-bottom:20px; box-shadow:0 1px 4px rgba(22,50,79,0.08); overflow:hidden;",
          tags$summary(
            style = "font-weight:700; font-size:14.5px; color:#ffffff; background:#16324f; cursor:pointer; display:flex; align-items:center; gap:8px; padding:12px 18px; list-style:none;",
            bs_icon("sliders", size = "1.1em"),
            "Filter Bin Metrics",
            tags$span(style = "margin-left:auto; font-weight:400; font-size:11px; color:#c7d2da;", "Click to expand / collapse")
          ),
          div(
            style = "padding:18px 20px 20px 20px;",
            
            # Chromosome selection
            div(
              style = "margin-bottom:16px;",
              div(
                style = "display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:8px; margin-bottom:8px;",
                tags$span(
                  style = "font-weight:700; font-size:11.5px; text-transform:uppercase; letter-spacing:0.6px; color:#7d92a3; display:flex; align-items:center; gap:6px;",
                  bs_icon("bar-chart-steps"), "Chromosome"
                ),
                div(
                  style = "display:flex; gap:6px;",
                  actionButton("binfilter_chr_all", "Select all", class = "btn-sm btn-outline-secondary", style = "font-size:11px; padding:2px 10px; border-radius:5px;"),
                  actionButton("binfilter_chr_clear", "Clear", class = "btn-sm btn-outline-secondary", style = "font-size:11px; padding:2px 10px; border-radius:5px;")
                )
              ),
              div(
                style = "background:#f7f9fa; border:1px solid #e7ecef; border-radius:8px; padding:10px 12px;",
                checkboxGroupInput("binfilter_chr", label = NULL, choices = chrom_list, selected = chrom_list, inline = TRUE)
              )
            ),
            
            tags$hr(style = "border-top:1px solid #eef2f5; margin:14px 0;"),
            
            # Gene search
            div(
              style = "margin-bottom:16px;",
              div(
                style = "display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:8px; margin-bottom:8px;",
                tags$span(
                  style = "font-weight:700; font-size:11.5px; text-transform:uppercase; letter-spacing:0.6px; color:#7d92a3; display:flex; align-items:center; gap:6px;",
                  bs_icon("search"), "Gene"
                ),
                actionButton("binfilter_gene_clear", "Clear", class = "btn-sm btn-outline-secondary", style = "font-size:11px; padding:2px 10px; border-radius:5px;")
              ),
              div(
                style = "background:#f7f9fa; border:1px solid #e7ecef; border-radius:8px; padding:10px 12px;",
                selectizeInput(
                  "binfilter_gene", label = NULL, choices = NULL, selected = character(0),
                  options = list(placeholder = "Search a gene, e.g. TP53...", maxOptions = 200),
                  width = "100%"
                ),
                tags$div(
                  style = "font-size:11px; color:#7d92a3; margin-top:4px;",
                  "Shows every bin whose coordinates overlap this gene (gene_start \u2264 bin_end AND gene_end \u2265 bin_start)."
                )
              )
            ),
            
            tags$hr(style = "border-top:1px solid #eef2f5; margin:14px 0;"),
            
            # Numeric metric filters
            tags$div(
              style = "font-weight:700; font-size:11.5px; text-transform:uppercase; letter-spacing:0.6px; color:#7d92a3; margin-bottom:10px; display:flex; align-items:center; gap:6px;",
              bs_icon("sliders2"), "Methylation & Content Metrics"
            ),
            div(
              style = "display:flex; flex-wrap:wrap; gap:16px; align-items:flex-start;",
              build_bin_table_filter_inputs(bin_table_filters, bin_table)
            ),
            
            tags$hr(style = "border-top:1px solid #eef2f5; margin:18px 0 14px 0;"),
            
            div(
              style = "display:flex; justify-content:flex-end;",
              actionButton(
                "binfilter_reset", "Reset all filters",
                icon = bs_icon("arrow-counterclockwise"),
                class = "btn-sm btn-outline-danger",
                style = "font-size:12px; border-radius:6px; padding:5px 14px;"
              )
            )
          )
        ),
        div(
          style = "display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:12px; margin-bottom:15px;",
          downloadButton(
            "bintable_download", 
            "Download CSV Data", 
            icon = icon("download"),
            class = "btn-lg btn-primary",
            style = "font-weight: 600; padding: 10px 22px; font-size: 15px; border-radius: 6px;"
          ),
          uiOutput("bintable_gene_note")
        ),
        div(style = "margin-top:4px;", DTOutput("bintable"))
      )
    )
  ),
  
  div(style = "display:flex; align-items:center; justify-content:center; gap:18px; margin-top:24px; padding:16px 0; border-top:1px solid #dde3e8; color:#6c757d; font-size:13px;", tags$span("Paula Artiz Dueñas, Bioinformatics Student", style = "margin-left:10px;"))
)

# Server Definition

server <- function(input, output, session) {
  
  updateSelectizeInput(session, "bins", choices = bin_table$bin_id,selected = character(0), server = TRUE)  
  updateSelectizeInput(session, "binfilter_gene", choices = all_gene_names, selected = character(0), server = TRUE)
  observeEvent(input$add_chr_bins, {
    chr_bins <- bin_choices_by_chr[[input$chr]]
    updated <- union(input$bins, chr_bins)
    updateSelectizeInput(session, "bins", choices = bin_table$bin_id, selected = updated, server = TRUE)
  })
  
  observeEvent(input$clear_bins, {
    updateSelectizeInput(session, "bins", choices = bin_table$bin_id, selected = character(0), server = TRUE)
  })
  
  observeEvent(input$bins_file, {
    req(input$bins_file)
    res <- parse_bin_file(input$bins_file$datapath, input$bins_file$name, bin_table, chrom_list)
    if (length(res$matched) == 0) {
      showNotification("Couldn't match any bins in file.", type = "error", duration = 7)
      return()
    }
    updated <- union(input$bins, res$matched)
    updateSelectizeInput(session, "bins", choices = bin_table$bin_id, selected = updated, server = TRUE)
  })
  
  selected_bins <- reactive({
    if (length(input$bins) > 0) input$bins else bin_choices_by_chr[[input$chr]]
  })
  
  selected_bin_table <- reactive({
    bin_table[bin_table$bin_id %in% selected_bins(), ]
  })
  
  observeEvent(input$binfilter_chr_all, {
    updateCheckboxGroupInput(session, "binfilter_chr", selected = chrom_list)
  })
  
  observeEvent(input$binfilter_chr_clear, {
    updateCheckboxGroupInput(session, "binfilter_chr", selected = character(0))
  })
  
  observeEvent(input$binfilter_gene_clear, {
    updateSelectizeInput(session, "binfilter_gene", choices = all_gene_names, selected = character(0), server = TRUE)
  })
  
  # resets every Bin Table filter back to its default
  observeEvent(input$binfilter_reset, {
    updateCheckboxGroupInput(session, "binfilter_chr", selected = chrom_list)
    updateSelectizeInput(session, "binfilter_gene", choices = all_gene_names, selected = character(0), server = TRUE)
    for (f in bin_table_filters) {
      vals <- bin_table[[f$id]]
      vals <- vals[is.finite(vals)]
      step <- if (!is.null(f$step)) f$step else 0.01
      lo <- if (length(vals) > 0) {
        if (step >= 1) floor(min(vals)) else floor(min(vals) * 100) / 100
      } else 0
      updateSliderInput(session, paste0("binfilter_", f$id), value = lo)
    }
  })
  
  # Bin Table tab has its own filter panel
  filtered_bin_table <- reactive({
    apply_bin_table_filters(bin_table, bin_table_filters, input)
  })
  
  # live "N of TOTAL bins" badge shown in the Bin Table card header
  output$bintable_count_badge <- renderUI({
    n <- nrow(filtered_bin_table())
    total <- nrow(bin_table)
    span(
      style = "font-size:12px; font-weight:600; color:#16324f; background:#eef2f5; border:1px solid #d8e0e6; border-radius:20px; padding:4px 12px; white-space:nowrap;",
      paste0(format(n, big.mark = ","), " of ", format(total, big.mark = ","), " bins")
    )
  })
  
  output$n_samples <- renderText({ nrow(metadata) })
  output$n_patients <- renderText({ length(unique(metadata$patient_id)) })
  
  # 1. Overview
  
  overview_ml_selected <- reactive({
    df <- methylation_long
    if ("bin_id" %in% names(df)) {
      df <- df[df$bin_id %in% selected_bins(), ]
    } else {
      df <- df[df$chr %in% unique(selected_bin_table()$chr), ]
    }
    merge(df, metadata[, c("sample_id", "patient_id", "Type", "sexe")], by = "sample_id")
  })
  
  # one row per patient: mean methylation (Tumor/Normal) and the Tumor-Normal shift
  overview_patient_shift <- reactive({
    df <- overview_ml_selected()
    req(nrow(df) > 0)
    agg <- aggregate(methylation ~ patient_id + Type, data = df, FUN = mean, na.rm = TRUE)
    wide <- reshape(agg, idvar = "patient_id", timevar = "Type", direction = "wide")
    names(wide) <- sub("^methylation\\.", "", names(wide))
    wide$shift <- wide$Tumor - wide$Normal
    merge(wide, patient_annotation, by = "patient_id")
  })
  
  # small label reused on every Overview plot so it's clear which region is shown
  region_label <- reactive({
    df <- selected_bin_table()
    paste0(nrow(df), " bin(s) - chr", paste(unique(df$chr), collapse = ", "))
  })
  
  # Overview: click-to-enlarge for the 4 plots
  overview_plot_meta <- list(
    list(id = "overview_boxplot",title = "Tumor vs. Normal Methylation"),
    list(id = "overview_mutations", title = "Mutation Status vs. Methylation Shift"),
    list(id = "overview_composition", title = "Stage & MSI/MSS vs. Methylation Shift"),
    list(id = "overview_sex_comparison", title = "Sex Distribution")
  )
  
  # NULL means show the 2x2 grid, otherwise it's the output id of the zoomed-in plot
  overview_zoom <- reactiveVal(NULL)
  
  # clicking the expand icon, or clicking the plot itself, enlarges it
  lapply(overview_plot_meta, function(pm) {
    observeEvent(input[[paste0("expand_", pm$id)]], { overview_zoom(pm$id) }, ignoreInit = TRUE)
    observeEvent(input[[paste0(pm$id, "_click")]], { overview_zoom(pm$id) }, ignoreInit = TRUE)
  })
  
  observeEvent(input$overview_back, { overview_zoom(NULL) }, ignoreInit = TRUE)
  
  output$overview_content <- renderUI({
    zoom <- overview_zoom()
    
    if (is.null(zoom)) {
      # overview grid: same 2x2 layout as before, each card now has a small
      card_list <- lapply(overview_plot_meta, function(pm) {
        card(
          card_header(
            div(style = "display:flex; justify-content:space-between; align-items:center; width:100%;",
                tags$span(pm$title),
                actionButton(paste0("expand_", pm$id), NULL, icon = bs_icon("arrows-fullscreen"),
                             class = "btn-sm btn-outline-secondary", title = "Ampliar",
                             style = "padding:2px 7px; margin-left:auto; margin-right:-5px;")
            )
          ),
          plotOutput(pm$id, height = "300px", click = paste0(pm$id, "_click")),
          height = "350px"
        )
      })
      do.call(layout_columns, c(list(col_widths = c(6, 6), row_heights = c("425px")), card_list))
      
    } else {
      # enlarged single-plot view, with a button to go back to the 2x2 grid
      meta <- Filter(function(pm) identical(pm$id, zoom), overview_plot_meta)[[1]]
      tagList(
        div(style = "margin-bottom:10px;",
            actionButton("overview_back", "Go back to Overview", icon = bs_icon("arrow-left-circle"),
                         class = "btn-sm btn-outline-primary")),
        card(card_header(meta$title), plotOutput(meta$id, height = "600px"), height = "650px")
      )
    }
  })
  
  output$overview_sex_comparison <- renderPlot({
    # counted at the patient level (via overview_patient_shift), not per bin-row
    df <- overview_patient_shift()
    validate(need(!is.null(df) && nrow(df) > 0, "No data available for the current selection."))
    
    sex_labels <- c("Dona" = "Female", "Home" = "Male")
    sex_vals <- sex_labels[as.character(df$sexe)]
    sex_vals[is.na(sex_vals)] <- "Unknown"
    sex_counts <- as.data.frame(table(Sex = sex_vals))
    sex_counts <- sex_counts[sex_counts$Freq > 0, ]
    names(sex_counts)[2] <- "Count"
    validate(need(nrow(sex_counts) > 0, "No sex information available for the current selection."))
    sex_counts$Pct <- sex_counts$Count / sum(sex_counts$Count)
    sex_counts$Label <- paste0(sex_counts$Sex, "\n", round(100 * sex_counts$Pct), "%")
    
    ggplot(sex_counts, aes(x = 2, y = Count, fill = Sex)) +
      geom_col(width = 1, color = "white", linewidth = 1.2) +
      geom_text(aes(label = Label), position = position_stack(vjust = 0.5),
                size = 3.4, color = "white", fontface = "bold", lineheight = 0.9) +
      coord_polar(theta = "y") +
      xlim(0.2, 2.5) +
      scale_fill_manual(values = c("Female" = "#d1495b", "Male" = "#3aa9c9", "Unknown" = "#b7c1c9")) +
      labs(fill = NULL, subtitle = paste0(nrow(df), " patient(s) \u2013 ", region_label())) +
      theme_void(base_size = 12) +
      theme(
        legend.position = "bottom",
        plot.subtitle = element_text(size = rel(0.78), color = "#7d92a3", hjust = 0.5)
      )
  })
  
  output$overview_boxplot <- renderPlot({
    df <- selected_bin_table()
    validate(need(nrow(df) > 0, "No bins selected. Choose a chromosome or bins in the sidebar."))
    plot_df <- data.frame(
      Type = rep(c("Normal", "Tumor"), each = nrow(df)),
      Methylation = c(df$mean_methylation_normal, df$mean_methylation_tumor)
    )
    plot_df <- plot_df[is.finite(plot_df$Methylation), ]
    validate(need(nrow(plot_df) > 0, "No valid methylation values for the current selection."))
    plot_df$Type <- factor(plot_df$Type, levels = c("Normal", "Tumor"))
    
    ggplot(plot_df, aes(x = Type, y = Methylation, fill = Type)) +
      geom_jitter(width = 0.07, size = 0.9, alpha = 0.25, color = "#16324f") +
      geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.92, color = "#16324f", linewidth = 0.4) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4,
                   fill = "white", color = "#16324f", stroke = 0.8) +
      scale_fill_manual(values = c("Tumor" = "#d1495b", "Normal" = "#3aa9c9")) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.06))) +
      labs(x = NULL, y = "Mean methylation (per bin)", subtitle = region_label(),
           caption = "\u25c7 mean") +
      theme_app() +
      theme(legend.position = "none")
  })
  
  output$overview_mutations <- renderPlot({
    df <- overview_patient_shift()
    genes <- c("KRAS", "BRAF", "TP53")
    long_df <- do.call(rbind, lapply(genes, function(g) {
      data.frame(Gene = g, Status = as.character(df[[g]]), Shift = df$shift, stringsAsFactors = FALSE)
    }))
    long_df <- long_df[!is.na(long_df$Status) & is.finite(long_df$Shift), ]
    validate(need(nrow(long_df) > 0, "No mutation status data available for the current selection."))
    
    # convention for the KRAS/BRAF/TP53 fields: 0 = wild-type, 1 = mutant
    status_labels <- c("0" = "Wild-type", "1" = "Mutant")
    long_df$Status <- factor(
      ifelse(long_df$Status %in% names(status_labels), status_labels[long_df$Status], "Unknown"),
      levels = c("Wild-type", "Mutant", "Unknown")
    )
    n_by_gene <- table(long_df$Gene)
    long_df$GeneLabel <- factor(
      paste0(long_df$Gene, " (n=", n_by_gene[long_df$Gene], ")"),
      levels = paste0(genes, " (n=", n_by_gene[genes], ")")
    )
    
    ggplot(long_df, aes(x = Status, y = Shift, fill = Status)) +
      geom_hline(yintercept = 0, color = "#7d92a3", linewidth = 0.4, linetype = "dashed") +
      geom_jitter(width = 0.1, size = 0.9, alpha = 0.3, color = "#16324f") +
      geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.92, color = "#16324f", linewidth = 0.4) +
      scale_fill_manual(values = c("Wild-type" = "#3aa9c9", "Mutant" = "#d1495b", "Unknown" = "darkgray")) +
      facet_wrap(~ GeneLabel, nrow = 1, scales = "free_x") +
      labs(x = NULL, y = "Methylation shift (Tumor \u2212 Normal)", subtitle = region_label()) +
      theme_app() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
  })
  
  output$overview_composition <- renderPlot({
    df <- overview_patient_shift()
    feats <- c("estadi2", "MSS")
    long_df <- do.call(rbind, lapply(feats, function(f) {
      data.frame(Feature = f, Level = as.character(df[[f]]), Shift = df$shift, stringsAsFactors = FALSE)
    }))
    long_df <- long_df[!is.na(long_df$Level) & is.finite(long_df$Shift), ]
    validate(need(nrow(long_df) > 0, "No stage/MSS data available for the current selection."))
    
    agg <- aggregate(Shift ~ Feature + Level, data = long_df, FUN = mean)
    n_df <- aggregate(Shift ~ Feature + Level, data = long_df, FUN = length)
    names(n_df)[3] <- "n"
    agg <- merge(agg, n_df, by = c("Feature", "Level"))
    
    # reuse the nicer names already defined for the heatmap feature selector
    feature_display <- setNames(names(heatmap_feature_choices), heatmap_feature_choices)
    agg$FeatureLabel <- factor(feature_display[agg$Feature], levels = feature_display[feats])
    agg$Level <- factor(agg$Level, levels = natural_level_order(agg$Level))
    
    ggplot(agg, aes(x = Level, y = Shift, fill = Feature)) +
      geom_hline(yintercept = 0, color = "#7d92a3", linewidth = 0.4) +
      geom_col(width = 0.6, show.legend = FALSE) +
      geom_text(aes(y = Shift, label = paste0("n=", n), vjust = ifelse(Shift >= 0, -0.5, 1.3)),
                size = 2.7, color = "#16324f") +
      scale_fill_manual(values = c("estadi2" = "#d1495b", "MSS" = "#3aa9c9")) +
      scale_y_continuous(expand = expansion(mult = c(0.15, 0.18))) +
      facet_wrap(~ FeatureLabel, scales = "free_x", nrow = 1) +
      labs(x = NULL, y = "Mean methylation shift (Tumor \u2212 Normal)", subtitle = region_label()) +
      theme_app()
  })
  
  # 2. Genome Browser
  nav <- reactiveValues(chr = NULL, start = NULL, end = NULL)
  
  observeEvent(input$chr, {
    updateSelectInput(session, "browser_chr", selected = input$chr)
  })
  
  observeEvent(input$browser_chr, {
    if (!identical(nav$chr, input$browser_chr)) {
      nav$chr   <- input$browser_chr
      nav$start <- 1
      nav$end   <- chrom_lengths[[input$browser_chr]]
    }
  }, ignoreInit = FALSE)
  
  observeEvent(input$prev_chr, {
    idx <- match(nav$chr, chrom_list)
    if (!is.na(idx) && idx > 1) {
      updateSelectInput(session, "browser_chr", selected = chrom_list[idx - 1])
    } else {
      showNotification("Already at the first chromosome.", type = "message", duration = 2)
    }
  })
  
  observeEvent(input$next_chr, {
    idx <- match(nav$chr, chrom_list)
    if (!is.na(idx) && idx < length(chrom_list)) {
      updateSelectInput(session, "browser_chr", selected = chrom_list[idx + 1])
    } else {
      showNotification("Already at the last chromosome.", type = "message", duration = 2)
    }
  })
  
  observeEvent(input$zoom_in, {
    req(nav$chr)
    v <- zoom_view(nav$start, nav$end, chrom_lengths[[nav$chr]], factor = 0.5)
    nav$start <- v$start; nav$end <- v$end
  })
  
  observeEvent(input$zoom_out, {
    req(nav$chr)
    v <- zoom_view(nav$start, nav$end, chrom_lengths[[nav$chr]], factor = 2)
    nav$start <- v$start; nav$end <- v$end
  })
  
  observeEvent(input$zoom_reset, {
    req(nav$chr)
    nav$start <- 1
    nav$end <- chrom_lengths[[nav$chr]]
  })
  
  observeEvent(input$browser_search_go, {
    req(input$browser_search)
    res <- parse_genomic_search(input$browser_search, bin_table, chrom_list, gene_lookup_table)
    if (identical(res$status, "ok")) {
      nav$chr   <- res$chr # set this first so the browser_chr observer above sees it and skips its reset
      nav$start <- max(1, floor(res$start))
      nav$end   <- min(chrom_lengths[[res$chr]], ceiling(res$end))
      if (nav$end <= nav$start) nav$end <- min(chrom_lengths[[res$chr]], nav$start + 2e6)
      if (!identical(input$browser_chr, res$chr)) updateSelectInput(session, "browser_chr", selected = res$chr)
      showNotification(res$message, type = "message", duration = 4)
    } else {
      showNotification(res$message, type = "error", duration = 5)
    }
  })
  
  # clicking a point on the all-chromosomes overview jumps to that chromosome
  observeEvent(event_data("plotly_click", source = "genome_overview"), {
    ed <- event_data("plotly_click", source = "genome_overview")
    req(ed, ed$customdata[1] %in% chrom_list)
    updateSelectInput(session, "browser_chr", selected = ed$customdata[1])
  })
  
  output$browser_position_header <- renderText({
    req(nav$chr, nav$start, nav$end)
    paste0(
      "Genome Browser - chr", nav$chr, ": ",
      format(round(nav$start), big.mark = ",", scientific = FALSE), "-",
      format(round(nav$end), big.mark = ",", scientific = FALSE),
      " (", format_mb_label(nav$end - nav$start), ")"
    )
  })
  
  browser_view_bins <- reactive({
    req(nav$chr, nav$start, nav$end)
    df <- bin_table[bin_table$chr == nav$chr & bin_table$bin_end >= nav$start & bin_table$bin_start <= nav$end, ]
    df[order(df$bin_position), ]
  })
  
  output$genome_overview_plot <- renderPlotly({
    df <- genome_wide_bins
    p <- plot_ly(df, x = ~genome_x, y = ~overall_meth, color = ~chr_parity,
                 colors = c("#16324f", "#0e7c86"), customdata = ~as.character(chr),
                 type = "scatter", mode = "markers", marker = list(size = 4),
                 text = ~bin_id, hoverinfo = "text", source = "genome_overview") |>
      layout(
        showlegend = FALSE,
        xaxis = list(title = "", tickvals = as.numeric(chr_mid), ticktext = chrom_list, tickangle = 45),
        yaxis = list(title = "Mean methylation", range = c(0, 1)),
        shapes = chr_boundary_shapes) |>
      config(displayModeBar = FALSE) |>
      event_register("plotly_click")
  })
  
  output$browser_chr_plot <- renderPlotly({
    df <- browser_view_bins()
    if (nrow(df) == 0) {
      return(
        plotly_empty(type = "scatter", mode = "markers") |>
          layout(title = list(text = "No bins in this region", font = list(size = 14)))
      )
    }
    # only highlight bins the user actually picked in the sidebar
    hl <- df$bin_id %in% input$bins
    plot_ly(df, x = ~bin_position / 1e6) |>
      add_trace(y = ~mean_methylation_tumor, type = "scatter", mode = "lines+markers", name = "Tumor",
                line = list(color = "#d1495b"),
                marker = list(color = "#d1495b", size = ifelse(hl, 11, 6),
                              line = list(color = "#0b2436", width = ifelse(hl, 1.5, 0))),
                text = ~paste0(bin_id, "<br>Tumor mean: ", round(mean_methylation_tumor, 3)), hoverinfo = "text") |>
      add_trace(
        y = ~mean_methylation_normal, type = "scatter", mode = "lines+markers", name = "Normal",
        line = list(color = "#3aa9c9"),
        marker = list(color = "#3aa9c9", size = ifelse(hl, 11, 6),
                      line = list(color = "#0b2436", width = ifelse(hl, 1.5, 0))),
        text = ~paste0(bin_id, "<br>Normal mean: ", round(mean_methylation_normal, 3)), hoverinfo = "text") |>
      layout(
        xaxis = list(title = paste0("Position on chr", nav$chr, " (Mb)")),
        yaxis = list(title = "Mean methylation", range = c(0, 1)),
        legend = list(orientation = "h", x = 0, y = 1, yanchor = "bottom")
      )
  })
  
  # 3. Tumor vs. Normal
  
  # sample x bin matrix of raw methylation for the currently selected bins, used by PCA/UMAP
  projection_matrix <- reactive({
    bins_now <- selected_bins()
    req(length(bins_now) > 0)
    df <- ml_annot[ml_annot$bin_id %in% bins_now, c("sample_id", "bin_id", "methylation")]
    mat <- build_wide_matrix(df, "sample_id", "bin_id", "methylation")
    req_ok <- !is.null(mat) && ncol(mat) >= 2
    if (!req_ok) return(NULL)
    keep_cols <- apply(mat, 2, function(x) stats::sd(x) > 0)
    mat <- mat[, keep_cols, drop = FALSE]
    if (nrow(mat) < 3 || ncol(mat) < 2) return(NULL)
    mat
  })
  
  # patient x bin matrix of methylation shift (Tumor - Normal) for the selected bins, used by the network plot
  patient_shift_matrix <- reactive({
    bins_now <- selected_bins()
    req(length(bins_now) > 0)
    df <- patient_bin_shift[patient_bin_shift$bin_id %in% bins_now, c("patient_id", "bin_id", "shift")]
    mat <- build_wide_matrix(df, "patient_id", "bin_id", "shift")
    if (is.null(mat) || nrow(mat) < 3 || ncol(mat) < 2) return(NULL)
    mat
  })
  
  output$tn_density <- renderPlot({
    df <- selected_bin_table()
    validate(need(nrow(df) > 0, "No bins selected. Choose a chromosome or bins in the sidebar."))
    plot_df <- data.frame(
      Type = rep(c("Tumor", "Normal"), each = nrow(df)),
      Methylation = c(df$mean_methylation_tumor, df$mean_methylation_normal)
    )
    plot_df <- plot_df[is.finite(plot_df$Methylation), ]
    validate(need(nrow(plot_df) > 0, "No valid methylation values for the current selection."))
    
    ggplot(plot_df, aes(x = Methylation, fill = Type, color = Type)) +
      geom_density(alpha = 0.35, linewidth = 0.9, na.rm = TRUE) +
      scale_fill_manual(values = c("Tumor" = "#d1495b", "Normal" = "#3aa9c9")) +
      scale_color_manual(values = c("Tumor" = "#d1495b", "Normal" = "#3aa9c9")) +
      labs(x = "Mean methylation (per bin)", y = "Density", title = region_label(), fill = NULL, color = NULL) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "top")
  })
  
  output$tn_projection <- renderPlotly({
    method <- input$proj_method
    if (is.null(method)) method <- "PCA"
    mat <- projection_matrix()
    validate(need(!is.null(mat),
                  "Not enough samples with complete data for the selected bins to compute a projection. Try selecting more bins."))
    
    if (identical(method, "PCA")) {
      pca <- prcomp(mat, center = TRUE, scale. = TRUE)
      coords <- as.data.frame(pca$x[, 1:2])
      names(coords) <- c("Dim1", "Dim2")
      var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
      xlab <- paste0("PC1 (", var_exp[1], "%)")
      ylab <- paste0("PC2 (", var_exp[2], "%)")
    } else {
      validate(need(requireNamespace("uwot", quietly = TRUE),
                    "UMAP requires the 'uwot' package, which isn't installed. Install it with install.packages('uwot'), or switch to PCA."))
      set.seed(42)
      n_neighbors <- max(2, min(15, nrow(mat) - 1))
      um <- uwot::umap(mat, n_neighbors = n_neighbors, n_components = 2)
      coords <- as.data.frame(um)
      names(coords) <- c("Dim1", "Dim2")
      xlab <- "UMAP 1"
      ylab <- "UMAP 2"
    }
    
    coords$sample_id <- rownames(mat)
    coords <- merge(coords, metadata[, c("sample_id", "patient_id", "Type", "sexe")], by = "sample_id")
    
    plot_ly(
      coords, x = ~Dim1, y = ~Dim2, color = ~Type,
      colors = c("Tumor" = "#d1495b", "Normal" = "#3aa9c9"),
      type = "scatter", mode = "markers",
      marker = list(size = 10, line = list(color = "#0b2436", width = 1)),
      text = ~paste0("Sample: ", sample_id, "<br>Patient: ", patient_id, "<br>Type: ", Type),
      hoverinfo = "text"
    ) |>
      layout(
        xaxis = list(title = xlab), yaxis = list(title = ylab),
        legend = list(orientation = "h", x = 0, y = 1.08)
      )
  })
  
  output$tn_network <- renderPlot({
    mat <- patient_shift_matrix()
    validate(need(!is.null(mat),
                  "Not enough overlapping data across patients for the selected bins to build a similarity network. Try selecting more bins."))
    
    cor_mat <- stats::cor(t(mat))
    layout_xy <- tryCatch(stats::cmdscale(stats::as.dist(1 - cor_mat), k = 2), error = function(e) NULL)
    validate(need(!is.null(layout_xy), "Couldn't compute a layout for the current selection."))
    
    nodes <- data.frame(patient_id = rownames(mat), x = layout_xy[, 1], y = layout_xy[, 2], stringsAsFactors = FALSE)
    nodes <- merge(nodes, patient_annotation, by = "patient_id", all.x = TRUE)
    nodes$mss_label <- as.character(nodes$MSS)
    nodes$mss_label[is.na(nodes$mss_label)] <- "unknown"
    mss_palette <- categorical_palette(nodes$mss_label)
    
    # only draw an edge between the most similar pairs of patients (top 15% by correlation), keeps the plot readable
    ids <- rownames(mat)
    edge_df <- NULL
    if (length(ids) >= 2) {
      pairs <- utils::combn(ids, 2)
      edge_df <- data.frame(
        from = pairs[1, ], to = pairs[2, ],
        corr = cor_mat[cbind(pairs[1, ], pairs[2, ])],
        stringsAsFactors = FALSE
      )
      edge_df <- edge_df[is.finite(edge_df$corr), ]
      thresh <- stats::quantile(edge_df$corr, 0.85, na.rm = TRUE)
      edge_df <- edge_df[edge_df$corr >= thresh, ]
      edge_df <- merge(edge_df, setNames(nodes[, c("patient_id", "x", "y")], c("from", "x_from", "y_from")), by = "from")
      edge_df <- merge(edge_df, setNames(nodes[, c("patient_id", "x", "y")], c("to", "x_to", "y_to")), by = "to")
    }
    
    p <- ggplot()
    if (!is.null(edge_df) && nrow(edge_df) > 0) {
      p <- p + geom_segment(
        data = edge_df, aes(x = x_from, y = y_from, xend = x_to, yend = y_to, alpha = corr),
        color = "#7d92a3", linewidth = 0.4, show.legend = FALSE)
    }
    p +
      geom_point(data = nodes, aes(x = x, y = y, fill = mss_label), shape = 21, size = 6, color = "#0b2436", stroke = 0.6) +
      geom_text(data = nodes, aes(x = x, y = y, label = patient_id), vjust = -1.4, size = 3, color = "#16324f") +
      scale_fill_manual(values = mss_palette, name = "MSS status") +
      labs(
        x = "Dimension 1", y = "Dimension 2",
        title = "Patient similarity network (methylation-shift correlation)",
        subtitle = region_label()
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank())
  })
  
  # 4. Genome-wide Profile
  output$manhattan_plot <- renderPlotly({
    df <- genome_wide_bins
    validate(need(nrow(df) > 0, "No genome-wide methylation data available."))
    
    cat_colors <- c(
      "Not significant" = "#b7c1c9",
      "Extreme outlier" = "#e0a339",
      "Significant hypomethylation" = "#3aa9c9",
      "Significant hypermethylation" = "#d1495b"
    )
    
    plot_ly(
      df, x = ~genome_x, y = ~mean_diff_tumor_normal, color = ~manhattan_category,
      colors = cat_colors, type = "scattergl", mode = "markers",
      marker = list(size = 5, opacity = 0.85, line = list(width = 0)),
      text = ~paste0(
        "Bin: ", bin_id,
        "<br>Chr: ", as.character(chr),
        "<br>Diff (Tumor \u2212 Normal): ", round(mean_diff_tumor_normal, 3),
        "<br>q-value: ", ifelse(is.na(q_value), "NA", format(signif(q_value, 3), scientific = TRUE)),
        "<br>Patients tested: ", ifelse(is.na(n), 0, n)
      ),
      hoverinfo = "text"
    ) |>
      layout(
        xaxis = list(title = "", tickvals = as.numeric(chr_mid), ticktext = chrom_list, tickangle = 45),
        yaxis = list(title = "Methylation difference (Tumor \u2212 Normal)", zeroline = TRUE, zerolinewidth = 1.4, zerolinecolor = "#16324f"),
        shapes = chr_boundary_shapes,
        legend = list(orientation = "h", x = 0, y = 1.12),
        hovermode = "closest"
      ) |>
      config(
        displayModeBar = TRUE,
        toImageButtonOptions = list(format = "png", filename = "genome_wide_methylation_profile")
      )
  })
  
  # 5. Feature × Chromosome Heatmap
  output$feature_heatmap <- renderPlot({
    feature <- input$heatmap_feature
    if (is.null(feature)) feature <- "none"
    
    bins_now <- selected_bins()
    validate(need(length(bins_now) > 0, "Select at least one bin in the sidebar."))
    
    # restrict to the bins currently selected in the sidebar, kept in genome order
    ord_bins <- bin_table$bin_id[bin_table$bin_id %in% bins_now]
    shift_df <- patient_bin_shift[patient_bin_shift$bin_id %in% bins_now, ]
    shift_df$bin_id <- factor(shift_df$bin_id, levels = ord_bins)
    validate(need(nrow(shift_df) > 0, "No methylation shift data for the selected bins."))
    max_abs  <- max(abs(shift_df$shift), na.rm = TRUE)
    
    # row order: cluster by shift pattern if no feature is chosen, otherwise group by the feature value
    if (identical(feature, "none")) {mat <- build_shift_matrix(shift_df, ord_bins)
    mat[is.na(mat)] <- 0
    ord <- if (nrow(mat) > 1) rownames(mat)[hclust(dist(mat))$order] else rownames(mat)}
    else {
      ann_ord <- patient_annotation[order(patient_annotation[[feature]], as.numeric(patient_annotation$patient_id)), ]
      ord <- ann_ord$patient_id}
    
    shift_df$patient_id <- factor(shift_df$patient_id, levels = ord)
    
    p_main <- ggplot(shift_df, aes(x = bin_id, y = patient_id, fill = shift)) +
      geom_tile(color = "white", linewidth = 0.15) +
      scale_fill_gradient2(
        low = "#3aa9c9", mid = "white", high = "#d1495b", midpoint = 0,
        limits = c(-max_abs, max_abs), name = "Methylation\nshift\n(Tumor \u2212 Normal)") +
      scale_y_discrete(limits = rev(ord)) +
      labs(x = "Bin", y = "Patient") +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
    
    if (identical(feature, "none")) {
      return(p_main)
    }
    
    # colour strip for the chosen clinical feature, row-aligned with the main heatmap
    ann_df <- patient_annotation
    ann_df$patient_id <- factor(ann_df$patient_id, levels = ord)
    ann_df$x <- "Feature"
    ann_df$value <- as.character(ann_df[[feature]])
    ann_df$value[is.na(ann_df$value)] <- "unknown"
    
    levels_cat <- sort(unique(ann_df$value))
    base_palette <- c("#0e7c86", "#d1495b", "#3aa9c9", "#e0a339", "#2fae66", "#16324f", "#7d92a3")
    pal <- setNames(base_palette[seq_along(levels_cat)], levels_cat)
    
    feature_label <- names(heatmap_feature_choices)[heatmap_feature_choices == feature]
    
    p_ann <- ggplot(ann_df, aes(x = x, y = patient_id, fill = value)) +
      geom_tile(color = "white", linewidth = 0.15) +
      scale_fill_manual(values = pal, name = feature_label) +
      scale_y_discrete(limits = rev(ord)) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 12) +
      theme(axis.text.y = element_blank(), panel.grid = element_blank())
    
    p_ann + p_main + plot_layout(widths = c(1, 20)) # plot_layout() from patchwork package
  })
  
  # 6. Clinical Explorer
  
  # table data as its own reactive so DT row indices (input$clinical_table_rows_selected)
  clinical_table_data <- reactive({
    clinical_profile_table
  })
  
  clinical_table_proxy <- dataTableProxy("clinical_table")
  
  selected_patient_id <- reactive({
    sel <- input$clinical_table_rows_selected
    if (is.null(sel) || length(sel) == 0) return(NULL)
    as.character(clinical_table_data()$patient_id[sel])
  })
  
  observeEvent(input$clear_patient_selection, {
    selectRows(clinical_table_proxy, NULL)
  })
  
  output$clinical_boxplot <- renderPlot({
    gene <- input$mutation_gene
    req(gene)
    test_method <- input$mutation_stat_test
    if (is.null(test_method)) test_method <- "wilcox"
    
    df <- overview_patient_shift()
    validate(need(!is.null(df) && nrow(df) > 0, "No data available for the current bin selection."))
    validate(need(gene %in% names(df), paste0(gene, " status is not available in this dataset.")))
    
    status_labels <- c("0" = "Wild-type", "1" = "Mutant")
    status_raw <- as.character(df[[gene]])
    df$Status <- ifelse(status_raw %in% names(status_labels), status_labels[status_raw], NA_character_)
    df <- df[!is.na(df$Status) & is.finite(df$Tumor), ]
    validate(need(length(unique(df$Status)) == 2,
                  paste0("Need both Wild-type and Mutant tumor samples for ", gene, " in the current bin selection.")))
    df$Status <- factor(df$Status, levels = c("Wild-type", "Mutant"))
    
    n_by_status <- table(df$Status)
    test_res <- run_group_test(df$Tumor[df$Status == "Mutant"], df$Tumor[df$Status == "Wild-type"], method = test_method)
    y_max <- max(df$Tumor, na.rm = TRUE)
    
    p <- ggplot(df, aes(x = Status, y = Tumor, fill = Status)) +
      geom_jitter(width = 0.08, size = 1.1, alpha = 0.35, color = "#16324f") +
      geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.92, color = "#16324f", linewidth = 0.4) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = "white", color = "#16324f", stroke = 0.8) +
      scale_fill_manual(values = c("Wild-type" = "#3aa9c9", "Mutant" = "#d1495b")) +
      scale_x_discrete(labels = paste0(levels(df$Status), "\n(n=", as.numeric(n_by_status[levels(df$Status)]), ")")) +
      annotate("text", x = 1.5, y = y_max * 1.1, label = test_res$label, size = 3.5, color = "#16324f", fontface = "bold") +
      scale_y_continuous(limits = c(0, max(1, y_max * 1.18)), expand = expansion(mult = c(0.02, 0.02))) +
      labs(
        x = NULL, y = "Mean tumor methylation (per bin)",
        title = paste0(gene, " mutation status vs. tumor methylation"),
        subtitle = region_label()
      ) +
      theme_app() +
      theme(legend.position = "none")
    
    sel_id <- selected_patient_id()
    if (!is.null(sel_id) && sel_id %in% df$patient_id) {
      hl <- df[df$patient_id == sel_id, ]
      p <- p +
        geom_point(data = hl, aes(x = Status, y = Tumor), shape = 21, size = 4.8,
                   fill = "#e0a339", color = "#16324f", stroke = 1.1) +
        geom_text(data = hl, aes(x = Status, y = Tumor, label = paste0("Patient ", patient_id)),
                  vjust = -1.3, size = 3, color = "#16324f", fontface = "bold")
    }
    p
  })
  
  output$clinical_table <- renderDT({
    df <- clinical_table_data()
    display_df <- df
    names(display_df) <- vapply(names(df), clinical_label_for, character(1))
    num_cols <- vapply(display_df, is.numeric, logical(1))
    display_df[num_cols] <- lapply(display_df[num_cols], round, 3)
    
    datatable(
      display_df,
      filter = "top",
      selection = "single",
      rownames = FALSE,
      options = list(scrollX = TRUE, pageLength = 10, autoWidth = TRUE),
      class = "stripe hover compact"
    )
  })
  
  output$selected_patient_banner <- renderUI({
    sel_id <- selected_patient_id()
    if (is.null(sel_id)) {
      return(div(
        style = "font-size:12px; color:#7d92a3; background:#eef2f5; border:1px dashed #c3ccd3; border-radius:8px; padding:10px 16px; margin-bottom:14px;",
        bs_icon("info-circle"), " No patient selected \u2014 click a row in the clinical metadata table to highlight that patient in the plot."
      ))
    }
    row <- clinical_profile_table[clinical_profile_table$patient_id == sel_id, ]
    if (nrow(row) == 0) return(NULL)
    
    chip <- function(col_id, fallback = "\u2014") {
      raw <- if (col_id %in% names(row) && length(row[[col_id]]) == 1) row[[col_id]] else NA
      val <- if (is.na(raw)) fallback else if (is.numeric(raw)) format(round(raw, 3), nsmall = 3) else as.character(raw)
      div(
        style = "display:flex; flex-direction:column; gap:2px;",
        tags$span(clinical_label_for(col_id), style = "font-size:10px; color:#7d92a3; text-transform:uppercase; letter-spacing:0.4px;"),
        tags$span(val, style = "font-size:14px; font-weight:600; color:#16324f;")
      )
    }
    
    chip_cols <- intersect(c("sexe", "estadi2", "MSS", "BRAF", "KRAS", "TP53", "Recaiguda", "MethShift"), names(row))
    
    div(
      style = "display:flex; align-items:center; gap:26px; flex-wrap:wrap; background:#eef2f5; border:1px solid #d8e0e6; border-radius:8px; padding:12px 18px; margin-bottom:14px;",
      div(
        style = "display:flex; align-items:center; gap:8px; font-weight:700; color:#16324f;",
        bs_icon("person-check-fill", size = "1.3em"), paste0("Patient ", sel_id)
      ),
      lapply(chip_cols, chip),
      actionButton("clear_patient_selection", "Clear selection", icon = bs_icon("x-circle"),
                   class = "btn-sm btn-outline-secondary", style = "margin-left:auto;")
    )
  })
  
  # 7. Bin Table
  
  # note shown if gene annotation couldn't be computed at all
  output$bintable_gene_note <- renderUI({
    if (!gene_annotation_available) {
      return(div(
        style = "font-size:12px; color:#a34; background:#fdf1ef; border:1px solid #f0d6d1; border-radius:8px; padding:8px 14px; max-width:520px;",
        bs_icon("exclamation-triangle"),
        " Gene annotation source files weren't found when the app started, so Gene Count / Genes / Gene IDs are unavailable. See the app startup warning for details."
      ))
    }
  })
  
  output$bintable <- renderDT({
    display_df <- build_display_bin_table(filtered_bin_table())
    
    meth_cols  <- intersect(c("Tumor Mean", "Normal Mean", "Tumor SD", "Normal SD", "\u0394 (Tumor \u2212 Normal)"), names(display_df))
    pct_cols   <- intersect(c("Alu Meth %", "CpG Meth %"), names(display_df))
    count_cols <- intersect(c("Alu Count", "CpG Count", "Gene Count"), names(display_df))
    delta_col  <- "\u0394 (Tumor \u2212 Normal)"
    
    dt <- datatable(
      display_df,
      container = bin_table_header_sketch,
      rownames = FALSE,
      class = "stripe hover compact",
      options = list(
        scrollX = TRUE,
        pageLength = 15,
        lengthMenu = list(c(10, 15, 25, 50, -1), c("10", "15", "25", "50", "All")),
        dom = "ltip",
        language = list(search = "Quick search:", lengthMenu = "Show _MENU_ bins per page"),
        columnDefs = list(list(className = "dt-center", targets = "_all"))
      )
    )
    
    if (length(meth_cols) > 0) dt <- formatRound(dt, meth_cols, 3)
    if (length(pct_cols) > 0) dt <- formatRound(dt, pct_cols, 2)
    if (length(count_cols) > 0) dt <- formatRound(dt, count_cols, 0)
    
    # colour the tumor/normal shift so a hyper- vs hypo-methylated bin is visible at a glance
    if (delta_col %in% names(display_df)) {
      dt <- formatStyle(
        dt, delta_col,
        fontWeight = "600",
        color = styleInterval(0, c("#0e7c86", "#d1495b"))
      )
    }
    # make bins that actually carry a curated CRC gene stand out from the (expected) majority with none
    if ("Gene Count" %in% names(display_df)) {
      dt <- formatStyle(
        dt, "Gene Count",
        fontWeight = styleInterval(0, c("normal", "700")),
        color = styleInterval(0, c("#adb8c0", "#16324f")),
        backgroundColor = styleInterval(0, c("transparent", "#e8f5f3"))
      )
    }
    dt
  })
  
  output$bintable_download <- downloadHandler(
    filename = function() {"bintable.csv"},
    content = function(file) {
      write.csv(build_display_bin_table(filtered_bin_table(), plain = TRUE), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)