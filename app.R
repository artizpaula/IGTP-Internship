# Exploration of Epigenomic Data in Colorectal Cancer
# Paula Artiz Dueñas - UPC

library(shiny)
library(bslib)
library(bsicons) # icons for value boxes
library(DT)
library(plotly)
library(ggplot2)
library(patchwork) # aligns the annotation strip + heatmap with independent legends

data <- readRDS("data_app.rds") # this file has to be in the same folder as app.R

metadata <- data$metadata 
bin_table <- data$bin_table
methylation_long <- data$methylation_long
chrom_list <- c(as.character(1:22), "X", "Y")

# order the bins so the dropdown in the sidebar makes sense (by chromosome then position)
bin_table <- bin_table[order(match(bin_table$chr, chrom_list), bin_table$bin_position),]
bin_choices_by_chr <- split(bin_table$bin_id, factor(bin_table$chr, levels = chrom_list))

# Bin Table tab (extra things needed)

# turn a start/end pair into something like "1 MB-2 MB" so it's readable
format_mb_label <- function(bp) {
  mb <- bp/1e6
  if (abs(mb - round(mb)) < 1e-9) paste0(round(mb)," MB")
  else paste0(format(round(mb, 1), nsmall = 1)," MB")
}
format_bin_coordinates <- function(start, end) {
  paste0(format_mb_label(start), "\u2013", format_mb_label(end))
}
bin_table$bin_coordinates <- mapply(format_bin_coordinates, bin_table$bin_start, bin_table$bin_end)

# difference between tumor and normal mean
bin_table$mean_diff_tumor_normal <- bin_table$mean_methylation_tumor - bin_table$mean_methylation_normal

# whether any bin actually ended up with gene annotation
gene_annotation_available <- any(!is.na(bin_table$gene_names))

# list of every gene symbol that shows up somewhere in the table
all_gene_names <- sort(unique(trimws(unlist(strsplit(na.omit(bin_table$gene_names), ";", fixed = TRUE)))))

# whether any bin actually ended up with promoter annotation
promoter_annotation_available <- "promoter_count" %in% names(bin_table) && any(bin_table$promoter_count > 0, na.rm = TRUE)

# columns to show in the Bin Table tab, in the order I want them
bin_table_columns <- list(list(id = "bin_id", label = "Bin ID", group = "Bin"),
                          list(id = "bin_coordinates", label = "Coordinates",group = "Bin"),
                          list(id = "mean_methylation_tumor", label = "Tumor Mean", group = "Methylation Values", digits = 3),
                          list(id = "mean_methylation_normal",label = "Normal Mean", group = "Methylation Values", digits = 3),
                          list(id = "sd_methylation_tumor",label = "Tumor SD", group = "Methylation Values", digits = 3),
                          list(id = "sd_methylation_normal",label = "Normal SD", group = "Methylation Values", digits = 3),
                          list(id = "mean_diff_tumor_normal", label = "\u0394 (Tumor \u2212 Normal)", group = "Methylation Values", digits = 3),
                          list(id = "n_alu", label = "Alu Count", group = "Genomic Content"),
                          list(id = "n_cpg", label = "CpG Count", group = "Genomic Content"),
                          list(id = "alu_methylation_pct", label = "Alu Meth %",group = "Genomic Content", digits = 2),
                          list(id = "cpg_methylation_pct",label = "CpG Meth %", group = "Genomic Content", digits = 2),
                          list(id = "gene_count", label = "Gene Count",group = "Gene Annotation (COSMIC Cancer Gene Census)"),
                          list(id = "gene_names", label = "Genes", group = "Gene Annotation (COSMIC Cancer Gene Census)"),
                          list(id = "gene_ids", label = "Gene IDs (COSMIC)", group = "Gene Annotation (COSMIC Cancer Gene Census)"),
                          list(id = "promoter_count", label = "Promoter Count", group = "Promoter Annotation"),
                          list(id = "quick_links_html", label = "Quick Links", group = "External Resources", html = TRUE, plain_id = "quick_links_text"))

# two-row header for the Bin Table DT (htmltools to make it fancier)
bin_table_header_sketch <- local({
  col_groups <- vapply(bin_table_columns, function(c) c$group, character(1))
  group_rle  <- rle(col_groups)
  htmltools::withTags(table(
    class = "display",
    thead(tr(lapply(seq_along(group_rle$lengths), function(i) {
      th(colspan = group_rle$lengths[i],
         style = "text-align:center; background:#eef2f5; color:#16324f; border-bottom:2px solid #d8e0e6; font-size:11.5px; text-transform:uppercase; letter-spacing:0.4px;",
         group_rle$values[i])
    })),
    tr(lapply(bin_table_columns, function(c) th(style = "white-space:nowrap;", c$label))))))
})

# numeric metrics we expose as filters
bin_table_filters <- list(list(id = "mean_methylation_tumor", label = "Mean Meth (Tumor)", op = ">=", step = 0.01),
                          list(id = "mean_methylation_normal", label = "Mean Meth (Normal)", op = ">=", step = 0.01),
                          list(id = "sd_methylation_tumor", label = "SD Meth (Tumor)", op = ">=", step = 0.01),
                          list(id = "sd_methylation_normal",label = "SD Meth (Normal)", op = ">=", step = 0.01),
                          list(id = "mean_diff_tumor_normal", label = "Mean Difference", op = ">=", step = 0.01),
                          list(id = "n_cpg", label = "CpG Count", op = ">=", step = 1),
                          list(id = "n_alu", label = "Alu Count", op = ">=", step = 1),
                          list(id = "alu_methylation_pct", label = "Alu Methylation %", op = ">=", step = 0.1),
                          list(id = "cpg_methylation_pct", label = "CpG Methylation %", op = ">=", step = 0.1),
                          list(id = "gene_count", label = "Gene Count", op = ">=", step = 1),
                          list(id = "promoter_count", label = "Promoter Count", op = ">=", step = 1))

# makes one slider per filter, using min/max from the actual data so the range makes sense
build_bin_table_filter_inputs <- function(filters, df) {
  lapply(filters, function(f) {
    vals <- df[[f$id]]
    vals <- vals[is.finite(vals)]
    step <- if (!is.null(f$step)) f$step else 0.01
    rng <- if (length(vals) > 0) range(vals) else c(0, 1)
    if (step >= 1) {
      lo <- floor(rng[1])
      hi <- ceiling(rng[2])}
    else {
      lo <- floor(rng[1] * 100) / 100
      hi <- ceiling(rng[2] * 100) / 100}
    if (hi <= lo) hi <- lo + step
    
    div(style = "flex: 1 1 200px; min-width: 180px; max-width: 250px;",
        sliderInput(inputId = paste0("binfilter_", f$id),
                    label = paste0(f$label, " ", f$op),
                    min = lo, max = hi, value = lo, step = step,
                    width = "100%"))})
}

# applies all the configured filters to a bin_table subset
apply_bin_table_filters <- function(df, filters, input) {
  # chromosome filter
  chr_selected <- input[["binfilter_chr"]]
  if (!is.null(chr_selected)) {
    df <- df[as.character(df$chr) %in% chr_selected, , drop = FALSE]}
  
  for (f in filters) {
    threshold <- input[[paste0("binfilter_", f$id)]]
    if (is.null(threshold)) next
    col <- df[[f$id]]
    if (is.null(col) || all(is.na(col))) next
    keep <- switch(f$op,
                   ">=" = !is.na(col) & col >= threshold,
                   ">"  = !is.na(col) & col >  threshold,
                   rep(TRUE, nrow(df)))
    df <- df[keep, , drop = FALSE]}
  
  # gene search: exact match
  gene_query <- input[["binfilter_gene"]]
  if (!is.null(gene_query) && nzchar(trimws(gene_query)) && "gene_names" %in% names(df)) {
    pattern <- paste0("(^|;)\\s*", toupper(trimws(gene_query)), "\\s*($|;)")
    keep <- !is.na(df$gene_names) & grepl(pattern, toupper(df$gene_names))
    df <- df[keep, , drop = FALSE]
  }
  df}


# builds the data frame that actually gets shown in the Bin Table, based on bin_table_columns
build_display_bin_table <- function(df, columns = bin_table_columns, plain = FALSE) {
  col_ids <- vapply(columns, function(c) c$id, character(1))
  out <- df[, col_ids, drop = FALSE]
  for (i in seq_along(columns)) {
    if (plain && !is.null(columns[[i]]$plain_id)) {
      out[[i]] <- df[[columns[[i]]$plain_id]]
    }
    if (!is.null(columns[[i]]$digits)) {
      out[[i]] <- round(out[[i]], columns[[i]]$digits)
    } else if (is.character(out[[i]])) {
      out[[i]][is.na(out[[i]])] <- "\u2014"
    }}
  names(out) <- vapply(columns, function(c) c$label, character(1))
  out}

# length of each chromosome
chrom_lengths <- tapply(bin_table$bin_end, bin_table$chr, max)
chrom_lengths <- setNames(as.numeric(chrom_lengths[chrom_list]), chrom_list)

# running offset so we can lay all chromosomes out on one axis
chr_offset <- setNames(numeric(length(chrom_list)), chrom_list)
running <- 0
for (chr in chrom_list) {
  chr_offset[chr] <- running
  running <- running + chrom_lengths[[chr]]}
chr_mid <- chr_offset + chrom_lengths / 2

# genome-wide bin table, used for the manhattan-style plots

genome_wide_bins <- bin_table
genome_wide_bins$chr <- factor(genome_wide_bins$chr, levels = chrom_list)
genome_wide_bins <- genome_wide_bins[order(genome_wide_bins$chr, genome_wide_bins$bin_position), ]
genome_wide_bins$genome_x <- chr_offset[as.character(genome_wide_bins$chr)] + genome_wide_bins$bin_position
genome_wide_bins$overall_meth <- rowMeans(genome_wide_bins[, c("mean_methylation_tumor", "mean_methylation_normal")], na.rm = TRUE)
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

# quick per-bin t-test (Tumor vs Normal) done "by hand" with mean/sd instead of t.test() (t.test() is too slow)
bin_stats <- aggregate(shift ~ bin_id, data = patient_bin_shift,
                       FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE), n = sum(is.finite(x)))
)
bin_stats <- do.call(data.frame, bin_stats)
names(bin_stats) <- c("bin_id", "mean_shift", "sd_shift", "n")
bin_stats$bin_id <- as.character(bin_stats$bin_id)
bin_stats$t_stat <- ifelse(bin_stats$n > 1 & bin_stats$sd_shift > 0, bin_stats$mean_shift / (bin_stats$sd_shift / sqrt(bin_stats$n)), NA_real_)
bin_stats$p_value <- ifelse(!is.na(bin_stats$t_stat) & bin_stats$n > 2, 2 * stats::pt(-abs(bin_stats$t_stat), df = pmax(bin_stats$n - 1, 1)), NA_real_)
bin_stats$q_value <- stats::p.adjust(bin_stats$p_value, method = "BH")

# one row per patient with the clinical fields we can use to annotate the heatmap (Recaiguda/BRAF/KRAS/TP53/MSS/sexe/estadi2 are the same for a patient's Tumor and Normal sample).
patient_annotation <- metadata[, c("patient_id", "Recaiguda", "BRAF", "KRAS", "TP53", "MSS", "sexe", "estadi2")]
patient_annotation <- patient_annotation[!duplicated(patient_annotation$patient_id), ]
patient_annotation$patient_id <- as.character(patient_annotation$patient_id)
patient_bin_shift$patient_id <- as.character(patient_bin_shift$patient_id)

# one row per SAMPLE (Tumor and Normal kept separate) for the patient similarity network
sample_annotation <- unique(metadata[, c("sample_id", "patient_id", "Type", "Recaiguda", "BRAF", "KRAS", "TP53", "MSS", "sexe", "estadi2")])
sample_annotation$patient_id <- as.character(sample_annotation$patient_id)

# merge the per-bin significance stats computed above into the genome-wide bin table
genome_wide_bins <- merge(genome_wide_bins, bin_stats[, c("bin_id", "mean_shift", "sd_shift", "n", "p_value", "q_value")], by = "bin_id", all.x = TRUE)
genome_wide_bins$chr <- factor(genome_wide_bins$chr, levels = chrom_list)
genome_wide_bins <- genome_wide_bins[order(genome_wide_bins$chr, genome_wide_bins$bin_position), ]

manhattan_sig_threshold <- 0.05
manhattan_outlier_threshold <- stats::quantile(abs(genome_wide_bins$mean_diff_tumor_normal), 0.95, na.rm = TRUE)

genome_wide_bins$manhattan_category <- with(genome_wide_bins, {
  sig <- !is.na(q_value) & q_value < manhattan_sig_threshold
  outlier <- !is.na(mean_diff_tumor_normal) & abs(mean_diff_tumor_normal) >= manhattan_outlier_threshold
  ifelse(sig & mean_diff_tumor_normal > 0, "Significant hypermethylation",
         ifelse(sig & mean_diff_tumor_normal < 0, "Significant hypomethylation",
                ifelse(outlier, "Extreme outlier", "Not significant")))})
genome_wide_bins$manhattan_category <- factor(genome_wide_bins$manhattan_category,
                                              levels = c("Not significant", "Extreme outlier", "Significant hypomethylation", "Significant hypermethylation"))

# options for the "which feature" dropdown on the heatmap tab
heatmap_feature_choices <- c("Methylation Values (No Association)" = "none", "Relapse"  = "Recaiguda", "BRAF status" = "BRAF", "KRAS status" = "KRAS", "TP53 status" = "TP53",
                             "MSI/MSS status"= "MSS","Sex" = "sexe","Stage" = "estadi2")

# same clinical variables, but for the similarity network's node colour (plus Sample Type, which isn't offered on the heatmap tab)
network_color_choices <- c(heatmap_feature_choices[heatmap_feature_choices != "none"], "Sample Type" = "Type")

# turns the 0/1 codes into words like "Wild-type"/"Mutant" so the tables read better
recode_binary <- function(x, labels = c("0" = "Wild-type", "1" = "Mutant")) {
  x_chr <- as.character(x)
  out <- ifelse(!is.na(x_chr) & x_chr %in% names(labels), labels[x_chr], x_chr)
  out
}

# Short p-value label used on under the plots
format_pvalue <- function(p) {
  if (is.null(p) || !is.finite(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", format(round(p, 3), nsmall = 3))
}

# runs either a Wilcoxon test or a t-test between two groups, whichever the user picks
run_group_test <- function(x, y, method = c("wilcox", "ttest")) {
  method <- match.arg(method)
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) {
    return(list(p_value = NA_real_, label = "Not enough samples for a statistical test",
                n_x = length(x), n_y = length(y), method = method))
  }
  res <- tryCatch(
    if (identical(method, "ttest")) stats::t.test(x, y) else stats::wilcox.test(x, y),
    error = function(e) NULL)
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
    all(tapply(metadata[[cn]], metadata$patient_id, function(v) length(unique(v[!is.na(v)])) <= 1))}, logical(1))
  candidate_cols[keep]
})
if (!"patient_id" %in% patient_level_cols) patient_level_cols <- c("patient_id", patient_level_cols)

clinical_profile_table <- unique(metadata[, patient_level_cols, drop = FALSE])
clinical_profile_table$patient_id <- as.character(clinical_profile_table$patient_id)
clinical_profile_table <- clinical_profile_table[!duplicated(clinical_profile_table$patient_id), ]

# recoding fields in catalan or 0/1 being used
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

# Nicer display names for the clinical table's columns
clinical_column_labels <- c(patient_id = "Patient ID", sexe = "Sex", estadi2 = "Stage", MSS = "MSI/MSS Status", BRAF = "BRAF", KRAS = "KRAS", TP53 = "TP53", Recaiguda = "Relapse",
                            edat_IQ = "Age",
                            MeanMethTumor = "Mean Methylation (Tumor)", MeanMethNormal = "Mean Methylation (Normal)",MethShift = "Methylation Shift (Tumor \u2212 Normal)",
                            nTumorSamples = "Tumor Samples", nNormalSamples = "Normal Samples")
clinical_label_for <- function(col_id) {
  if (col_id %in% names(clinical_column_labels)) return(unname(clinical_column_labels[col_id]))
  tools::toTitleCase(gsub("_", " ", col_id))
}

# Clinical Explorer filter panel: only the headline demographic/clinical fields the user asked for
# (sex, age, stage, relapse) get a filter control, even though the table itself still shows every column
clinical_numeric_cols <- intersect("edat_IQ", names(clinical_profile_table))
clinical_categorical_cols <- intersect(c("sexe", "estadi2", "Recaiguda"), names(clinical_profile_table))

# choices for a categorical filter, with missing values grouped into their own selectable "(Missing)" option
clinical_categorical_choices <- function(col_vals) {
  vals <- as.character(col_vals)
  vals[is.na(col_vals)] <- "(Missing)"
  sort(unique(vals))
}

# one checkbox group per categorical clinical column, all options selected by default
build_clinical_categorical_filter_inputs <- function(cols, df, compact = FALSE) {
  lapply(cols, function(col) {
    choices <- clinical_categorical_choices(df[[col]])
    if (compact) {
      div(style = "min-width:150px;",
          checkboxGroupInput(inputId = paste0("clinfilter_", col), label = clinical_label_for(col), choices = choices, selected = choices, inline = TRUE))
    } else {
      div(style = "flex: 1 1 220px; min-width: 200px; max-width: 320px;",
          tags$div(style = "font-weight:700; font-size:11.5px; text-transform:uppercase; letter-spacing:0.6px; color:#7d92a3; margin-bottom:6px;", clinical_label_for(col)),
          checkboxGroupInput(inputId = paste0("clinfilter_", col), label = NULL, choices = choices, selected = choices, inline = TRUE))
    }
  })
}

# one range slider per numeric clinical column (age, CNV counts, survival times, methylation summaries, etc.)
build_clinical_numeric_filter_inputs <- function(cols, df, compact = FALSE) {
  lapply(cols, function(col) {
    vals <- df[[col]]
    vals <- vals[is.finite(vals)]
    rng <- if (length(vals) > 0) range(vals) else c(0, 1)
    span <- rng[2] - rng[1]
    step <- if (span > 20 || all(vals == round(vals))) 1 else 0.01
    lo <- if (step >= 1) floor(rng[1]) else floor(rng[1] * 100) / 100
    hi <- if (step >= 1) ceiling(rng[2]) else ceiling(rng[2] * 100) / 100
    if (hi <= lo) hi <- lo + step
    width_style <- if (compact) "min-width:160px; max-width:200px;" else "flex: 1 1 220px; min-width: 200px; max-width: 280px;"
    div(style = width_style,
        sliderInput(inputId = paste0("clinfilter_", col), label = clinical_label_for(col),
                    min = lo, max = hi, value = c(lo, hi), step = step, width = "100%"))
  })
}

# applies every active Clinical Explorer filter to a clinical_profile_table subset
apply_clinical_filters <- function(df, numeric_cols, categorical_cols, input) {
  for (col in categorical_cols) {
    sel <- input[[paste0("clinfilter_", col)]]
    if (is.null(sel)) next
    vals <- as.character(df[[col]])
    vals[is.na(df[[col]])] <- "(Missing)"
    df <- df[vals %in% sel, , drop = FALSE]
  }
  for (col in numeric_cols) {
    rng <- input[[paste0("clinfilter_", col)]]
    if (is.null(rng) || length(rng) != 2) next
    col_vals <- df[[col]]
    keep <- is.na(col_vals) | (col_vals >= rng[1] & col_vals <= rng[2])
    df <- df[keep, , drop = FALSE]
  }
  df
}

# turns the long shift table into a patient x bin matrix
build_shift_matrix <- function(df, bin_ids) {
  m <- tapply(df$shift, list(df$patient_id, as.character(df$bin_id)), FUN = identity)
  keep_cols <- bin_ids[bin_ids %in% colnames(m)]
  m[, keep_cols, drop = FALSE]
}

# generic long-to-wide reshape (id, key, value) -> matrix, used for both PCA/UMAP and the network plot
build_wide_matrix <- function(df, id_col, key_col, value_col) {
  df <- df[is.finite(df[[value_col]]), c(id_col, key_col, value_col)]
  if (nrow(df) == 0) return(NULL)
  mat <- tapply(df[[value_col]], list(df[[id_col]], as.character(df[[key_col]])), FUN = identity)
  mat <- mat[, colSums(!is.finite(mat)) == 0, drop = FALSE]
  mat <- mat[rowSums(!is.finite(mat)) == 0, , drop = FALSE]
  if (nrow(mat) == 0 || ncol(mat) == 0) return(NULL)
  mat
}

# picks colors for however many categories a clinical variable has, recycling the palette
categorical_palette <- function(values) {
  base_palette <- c("#0e7c86", "#d1495b", "#3aa9c9", "#e0a339", "#2fae66", "#16324f", "#7d92a3")
  levels_now <- sort(unique(values))
  setNames(base_palette[((seq_along(levels_now) - 1) %% length(base_palette)) + 1], levels_now)
}

# one ggplot theme used everywhere so all the plots look the same instead of each one being styled separately
theme_app <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(plot.subtitle = element_text(size = rel(0.82), color = "#7d92a3", margin = margin(b = 10)),
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

# Converts an hclust object into a data.frame of line segments
hclust_to_segments <- function(hc) {
  n <- length(hc$order)
  leaf_pos <- integer(n)
  leaf_pos[hc$order] <- seq_len(n)
  
  merge_x <- numeric(nrow(hc$merge))
  merge_y <- hc$height
  
  # x/y of either a leaf (negative index into the original items) or an earlier merge
  get_xy <- function(idx) {
    if (idx < 0) list(x = leaf_pos[-idx], y = 0) else list(x = merge_x[idx], y = merge_y[idx])
  }
  
  segments <- vector("list", nrow(hc$merge) * 3)
  k <- 0
  for (i in seq_len(nrow(hc$merge))) {
    left  <- get_xy(hc$merge[i, 1])
    right <- get_xy(hc$merge[i, 2])
    merge_x[i] <- (left$x + right$x) / 2
    k <- k + 1; segments[[k]] <- data.frame(x = left$x,  y = left$y,     xend = left$x,  yend = merge_y[i])
    k <- k + 1; segments[[k]] <- data.frame(x = right$x, y = right$y,    xend = right$x, yend = merge_y[i])
    k <- k + 1; segments[[k]] <- data.frame(x = left$x,  y = merge_y[i], xend = right$x, yend = merge_y[i])
  }
  do.call(rbind, segments)
}

# builds a minimal ggplot dendrogram from a segment table produced by hclust_to_segments()
plot_dendrogram <- function(segments, orientation = c("top", "left")) {
  orientation <- match.arg(orientation)
  p <- ggplot(segments) + theme_void() + theme(plot.margin = margin(1, 1, 1, 1))
  if (orientation == "top") {
    p + geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "#7d92a3", linewidth = 0.3) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  } else {
    p + geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend), color = "#7d92a3", linewidth = 0.3) +
      scale_x_continuous(expand = expansion(mult = c(0.05, 0)))
  }
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
# gene set expanded (82 genes) using https://www.intogen.org/search and UCSC Genome Browser (manually done)

gene_lookup_table <- data.frame(
  gene = c("APC", "TP53", "KRAS", "PIK3CA", "SMAD4", "FBXW7", "SOX9", "AMER1",
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
  
  # Exemple 1: chr:start[-end] coordinate format
  coord_pattern <- "^(chr|Chr|CHR)?([0-9]{1,2}|[XxYy])\\s*:\\s*([0-9]+)\\s*(-\\s*([0-9]+))?$" # regular expression, learned in ap3 lecture
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
                       "-", format(end, big.mark = ",", scientific = FALSE))))
  }
  
  # Exemple 2: bin ID format "chr_position"
  bin_pattern <- "^([0-9]{1,2}|[XxYy])_([0-9]+)$" # regular expression, learned in ap3 lecture
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
  
  # Exemple 3: gene symbol already annotated on one or more bins
  gene_query <- toupper(q)
  has_annotation <- !is.na(bin_table$gene_names) &
    grepl(paste0("(^|[,;])\\s*", gene_query, "\\s*($|[,;])"), toupper(bin_table$gene_names))
  in_annotation <- bin_table[has_annotation,]
  if (nrow(in_annotation) >= 1) {
    return(list(status = "ok", chr = in_annotation$chr[1],
                start = min(in_annotation$bin_start), end = max(in_annotation$bin_end),
                message = paste0("Jumped to gene ", gene_query)))
  }
  
  # Exemple 4: gene symbol not in bin_table, so fall back to the static gene_lookup_table
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
      # single unlabeled column, could be a header row that's actually the first ID
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

# builds https URLs to external genomic resources for a given region
build_external_links <- function(chr, start, end, gene = NA_character_, genome_build = "hg19") {
  chr_lab <- paste0("chr", chr)
  ensembl_sub <- if (identical(genome_build, "hg19")) "grch37." else ""
  
  ucsc_url <- paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=", genome_build,"&position=", chr_lab, "%3A", format(round(start), scientific = FALSE),
                     "-", format(round(end), scientific = FALSE))
  ensembl_url <- paste0("https://", ensembl_sub, "ensembl.org/Homo_sapiens/Location/View?r=",chr, "%3A", format(round(start), scientific = FALSE), "-", format(round(end), scientific = FALSE) )
  # TCGA colorectal cohorts (COAD = colon, READ = rectal adenocarcinoma) via the GDC portal
  tcga_url <- "https://portal.gdc.cancer.gov/exploration?filters=%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-COAD%22%2C%22TCGA-READ%22%5D%7D%7D"
  
  links <- list(UCSC = ucsc_url, Ensembl = ensembl_url, TCGA = tcga_url)
  
  first_gene <- trimws(strsplit(as.character(gene)[1], "[,;]")[[1]])[1]
  if (!is.null(first_gene) && !is.na(first_gene) && nzchar(first_gene)) {
    links$Gene_info <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", utils::URLencode(first_gene), "%5Bsym%5D+AND+human%5Borgn%5D")
  }
  links
}

# Compact per-bin "Quick Links" HTML (UCSC, Ensembl, TCGA)
build_bin_quick_links_html <- function(chr, start, end, gene) {
  links <- build_external_links(chr, start, end, gene)
  tag_list <- list(
    tags$a(bs_icon("box-arrow-up-right", size = "0.8em"), " UCSC",
           href = links$UCSC, target = "_blank", rel = "noopener noreferrer",
           style = "margin-right:10px; font-size:12px; color:#2b6cb0; text-decoration:none; white-space:nowrap;"),
    tags$a(bs_icon("box-arrow-up-right", size = "0.8em"), " Ensembl",
           href = links$Ensembl, target = "_blank", rel = "noopener noreferrer",
           style = "margin-right:10px; font-size:12px; color:#0e7c86; text-decoration:none; white-space:nowrap;"),
    tags$a(bs_icon("box-arrow-up-right", size = "0.8em"), " TCGA",
           href = links$TCGA, target = "_blank", rel = "noopener noreferrer",
           style = "margin-right:10px; font-size:12px; color:#d1495b; text-decoration:none; white-space:nowrap;")
  )
  if (!is.null(links$Gene_info)) {
    tag_list <- c(tag_list, list(
      tags$a(bs_icon("box-arrow-up-right", size = "0.8em"), " NCBI Gene",
             href = links$Gene_info, target = "_blank", rel = "noopener noreferrer",
             style = "font-size:12px; color:#8a5a00; text-decoration:none; white-space:nowrap;")
    ))
  }
  as.character(tagList(tag_list))
}
bin_table$quick_links_html <- mapply(build_bin_quick_links_html,
                                     bin_table$chr, bin_table$bin_start, bin_table$bin_end, bin_table$gene_names,
                                     SIMPLIFY = TRUE
)

# plain-text equivalent of the column above (no HTML), used only for the CSV download
build_bin_quick_links_text <- function(chr, start, end, gene) {
  links <- build_external_links(chr, start, end, gene)
  parts <- c(paste0("UCSC: ", links$UCSC),
             paste0("Ensembl: ", links$Ensembl),
             paste0("TCGA: ", links$TCGA)
  )
  if (!is.null(links$Gene_info)) parts <- c(parts, paste0("NCBI Gene: ", links$Gene_info))
  paste(parts, collapse = " | ")
}
bin_table$quick_links_text <- mapply(
  build_bin_quick_links_text,
  bin_table$chr, bin_table$bin_start, bin_table$bin_end, bin_table$gene_names,
  SIMPLIFY = TRUE
)

# Plot description function
plot_desc <- function(...) {
  tags$div(style = "font-size:12px; color:#5c7182; margin-top:8px; padding-top:8px; border-top:1px solid #e3e8ec; line-height:1.45;", ...)
}

# Export helpers (PNG + CSV for every plot)

# small style shared by every export button so they all look the same
export_btn_style <- "font-size:11.5px; padding:4px 12px; border-radius:6px;"

# Row of export buttons placed under a ggplot/base-R plot (plotOutput)
plot_export_ui <- function(id) {
  div(style = "display:flex; gap:8px; margin-top:8px;",
      downloadButton(paste0(id, "_dl_png"), "Export plot (.png)", icon = icon("image"),
                     class = "btn-sm btn-outline-secondary", style = export_btn_style),
      downloadButton(paste0(id, "_dl_csv"), "Export data (.csv)", icon = icon("table"),
                     class = "btn-sm btn-outline-secondary", style = export_btn_style))
}

# Row of export buttons placed under a plotly plot (plotlyOutput)
plotly_export_ui <- function(id) {
  div(style = "display:flex; gap:8px; margin-top:8px;",
      tags$button(type = "button", class = "btn btn-sm btn-outline-secondary", style = export_btn_style,
                  onclick = sprintf("Plotly.downloadImage(document.getElementById('%s'), {format:'png', filename:'%s'});", id, id),
                  icon("image"), " Export plot (.png)"),
      downloadButton(paste0(id, "_dl_csv"), "Export data (.csv)", icon = icon("table"),
                     class = "btn-sm btn-outline-secondary", style = export_btn_style))
}

# Registers the two downloadHandlers behind plot_export_ui() for a ggplot-based plot
register_plot_export <- function(output, id, plot_fn, data_fn, width = 9, height = 6.5, dpi = 150) {
  output[[paste0(id, "_dl_png")]] <- downloadHandler(
    filename = function() paste0(id, ".png"),
    content = function(file) {
      ggsave(file, plot = plot_fn(), width = width, height = height, dpi = dpi, bg = "white")
    })
  output[[paste0(id, "_dl_csv")]] <- downloadHandler(
    filename = function() paste0(id, ".csv"),
    content = function(file) {
      write.csv(data_fn(), file, row.names = FALSE)
    })
}

# Registers just the CSV downloadHandler behind plotly_export_ui() (PNG is handled client-side).
register_plotly_export <- function(output, id, data_fn) {
  output[[paste0(id, "_dl_csv")]] <- downloadHandler(
    filename = function() paste0(id, ".csv"),
    content = function(file) {
      write.csv(data_fn(), file, row.names = FALSE)
    })
}

# Small bundle of responsive CSS rules for the Genome Browser toolbar + plots
gb_responsive_css <- "
.gb-toolbar { display:flex; flex-wrap:wrap; gap:14px; align-items:center; width:100%; }
.gb-toolbar-group { display:flex; gap:6px; align-items:center; }
.gb-search-group { flex:1 1 320px; min-width:220px; }
.gb-nav-group, .gb-zoom-group { flex:0 0 auto; }
.gb-plot-wrap { width:100%; }
.gb-plot-wrap .html-widget.plotly { width:100% !important; }
@media (max-width: 767.98px) {
  .gb-toolbar { gap:10px; }
  .gb-toolbar-group { width:100%; }
  .gb-nav-group, .gb-zoom-group { justify-content:center; }
  .gb-search-group { min-width:0; }
}
@media (max-width: 991.98px) {
  .gb-zoom-group actionButton, .gb-zoom-group .btn { flex: 1 1 auto; }
}

/* App-wide responsive rules: plots/tables/cards adapt to viewport instead of staying at a fixed desktop size */
.shiny-plot-output, .html-widget, .html-widget.plotly { max-width: 100%; }
.card-body { min-width: 0; } /* lets flex/grid children shrink instead of forcing horizontal scroll on the whole card */
.dataTables_wrapper { width: 100%; overflow-x: auto; -webkit-overflow-scrolling: touch; }
@media (max-width: 767.98px) {
  .card-header { font-size: 13.5px; padding: 10px 12px; }
  .card-body { padding: 12px; }
  .dataTables_wrapper .dataTables_filter input,
  .dataTables_wrapper .dataTables_length select { font-size: 12px; }
  table.dataTable { font-size: 12px; }
  .bslib-value-box .value-box-value { font-size: 20px; }
}
@media (max-width: 575.98px) {
  table.dataTable { font-size: 11px; }
}
"

# CSS for the Home / Overview landing page
home_css <- "
.home-page { padding: 4px 2px 6px 2px; }
.home-hero {
  position: relative; width: calc(100% + 3rem); margin: -1.5rem -1.5rem 30px -1.5rem;
  padding: 44px 48px 38px 48px; overflow: hidden; color: #eaf2f6;
  display: flex; flex-wrap: wrap; gap: 36px; align-items: center; justify-content: space-between;
  background: radial-gradient(1100px 460px at 12% -15%, rgba(58,169,201,0.30), transparent 60%),
              linear-gradient(135deg, #0b2436 0%, #16324f 55%, #0e5a63 130%);
  min-height: 300px;
}
.home-hero::after { content:''; position:absolute; right:-70px; top:-90px; width:300px; height:300px;
  border-radius:50%; background:rgba(58,169,201,0.16); }
.home-eyebrow { font-size:12px; font-weight:700; letter-spacing:1.6px; color:#7fd8cf; text-transform:uppercase; margin-bottom:14px; }
.home-hero-text { max-width:600px; position:relative; z-index:1; }
.home-hero-text h1 { font-size:clamp(24px,2.6vw,36px); font-weight:700; line-height:1.18; margin-bottom:14px; color:#ffffff; }
.home-hero-lede { font-size:14.5px; line-height:1.65; color:#c7d8e3; margin-bottom:22px; }
.home-hero-actions { display:flex; gap:12px; flex-wrap:wrap; }
.home-cta-btn { border:none; border-radius:8px; padding:10px 20px; font-weight:600; font-size:13.5px; cursor:pointer;
  transition:transform .15s ease, box-shadow .15s ease; display:inline-flex; align-items:center; gap:8px; }
.home-cta-primary { background:#0e7c86; color:#fff; box-shadow:0 6px 18px rgba(14,124,134,0.45); }
.home-cta-primary:hover { transform:translateY(-2px); box-shadow:0 10px 22px rgba(14,124,134,0.55); }
.home-cta-ghost { background:rgba(255,255,255,0.06); color:#eaf2f6; border:1px solid rgba(255,255,255,0.35); }
.home-cta-ghost:hover { background:rgba(255,255,255,0.14); transform:translateY(-2px); }
.home-hero-stats { display:grid; grid-template-columns:repeat(2,minmax(128px,1fr)); gap:12px; position:relative; z-index:1; }
.home-stat-card { background:rgba(255,255,255,0.08); border:1px solid rgba(255,255,255,0.16); border-radius:12px; padding:14px 18px; min-width:136px; }
.home-stat-value { font-size:25px; font-weight:700; color:#ffffff; line-height:1; }
.home-stat-label { font-size:11px; color:#a9c1d1; text-transform:uppercase; letter-spacing:0.5px; margin-top:6px; }
.home-section { margin-bottom:30px; }
.home-section-title { font-size:19px; font-weight:700; color:#16324f; margin-bottom:6px; display:flex; align-items:center; gap:10px; }
.home-section-lede { font-size:13px; color:#6c8093; margin-bottom:16px; max-width:840px; line-height:1.55; }
.home-grid-2 { display:grid; grid-template-columns:1.05fr 1fr; gap:20px; }
.home-card { background:#ffffff; border:1px solid #e3e8ec; border-radius:14px; padding:22px 24px; box-shadow:0 2px 10px rgba(15,40,60,0.04); }
.home-card h3 { font-size:15px; font-weight:700; color:#16324f; display:flex; align-items:center; gap:9px; margin-bottom:12px; }
.home-card p { font-size:13px; color:#3a5470; line-height:1.65; margin-bottom:10px; }
.home-card p:last-child { margin-bottom:0; }
.home-source-list { list-style:none; padding:0; margin:0; display:flex; flex-direction:column; gap:13px; }
.home-source-item { display:flex; gap:12px; align-items:flex-start; }
.home-source-icon { flex:0 0 auto; width:32px; height:32px; border-radius:9px; background:#eaf6f4; color:#0e7c86;
  display:flex; align-items:center; justify-content:center; font-size:15px; }
.home-source-text b { display:block; font-size:12.5px; color:#16324f; margin-bottom:2px; }
.home-source-text span { font-size:12px; color:#6c8093; line-height:1.5; }
.home-func-grid { display:grid; grid-template-columns:repeat(auto-fit, minmax(225px,1fr)); gap:14px; }
.home-func-card { background:#ffffff; border:1px solid #e3e8ec; border-radius:14px; padding:17px 19px; cursor:pointer;
  transition:all .15s ease; text-align:left; display:block; width:100%; }
.home-func-card:hover { border-color:#0e7c86; box-shadow:0 8px 20px rgba(14,124,134,0.14); transform:translateY(-3px); }
.home-func-icon { width:36px; height:36px; border-radius:10px; background:linear-gradient(135deg,#0e7c86,#3aa9c9); color:#fff;
  display:flex; align-items:center; justify-content:center; font-size:16px; margin-bottom:11px; }
.home-func-card h4 { font-size:13.5px; font-weight:700; color:#16324f; margin-bottom:5px; }
.home-func-card p { font-size:12px; color:#6c8093; line-height:1.55; margin-bottom:0; }
.home-func-arrow { font-size:11.5px; color:#0e7c86; font-weight:600; margin-top:9px; display:inline-flex; align-items:center; gap:4px; }
.home-flow { display:flex; align-items:stretch; gap:0; flex-wrap:wrap; }
.home-flow-stage { flex:1 1 200px; background:#ffffff; border:1px solid #e3e8ec; border-radius:14px; padding:19px 18px; min-width:200px; }
.home-flow-stage.home-flow-data { border-top:4px solid #3aa9c9; }
.home-flow-stage.home-flow-process { border-top:4px solid #0e7c86; }
.home-flow-stage.home-flow-modules { border-top:4px solid #16324f; }
.home-flow-stage.home-flow-insight { border-top:4px solid #2fae66; }
.home-flow-stage h5 { font-size:12px; font-weight:700; text-transform:uppercase; letter-spacing:0.5px; color:#7d92a3; margin-bottom:10px; }
.home-flow-stage ul { padding-left:17px; margin-bottom:0; font-size:12px; color:#3a5470; line-height:1.65; }
.home-flow-connector { flex:0 0 44px; display:flex; align-items:center; justify-content:center; color:#0e7c86; font-size:19px; }
.home-footer { width:calc(100% + 3rem); margin:6px -1.5rem -1.5rem -1.5rem; padding:14px 48px; background:#f7f9fa;
  border-top:1px solid #e7ecef; font-size:11px; color:#8092a3; text-align:center; }
@media (max-width: 991.98px) {
  .home-hero { flex-direction:column; align-items:flex-start; }
  .home-hero-stats { width:100%; grid-template-columns:repeat(2,1fr); }
  .home-grid-2 { grid-template-columns:1fr; }
  .home-flow { flex-direction:column; }
  .home-flow-connector { transform:rotate(90deg); padding:6px 0; }
}
"

# Home page helper builders (small functions so the markup below stays readable)

# One stat chip in the hero (e.g. "128 Patients"). `value_ui` can be a plain number/string
home_stat <- function(value_ui, label) {
  div(class = "home-stat-card",
      div(class = "home-stat-value", value_ui),
      div(class = "home-stat-label", label))
}

# One row in the "Data Sources" card.
home_source_item <- function(title, desc) {
  div(class = "home-source-item",
      div(class = "home-source-text", tags$b(title), tags$span(desc)))
}

# One clickable card in the "Key Functionalities" grid
home_func_card <- function(title, desc, nav_target) {
  tags$button(type = "button", class = "home-func-card",
              onclick = sprintf("Shiny.setInputValue('home_nav_target', '%s', {priority:'event'});", nav_target),
              h4(title), p(desc),
              div(class = "home-func-arrow", "Open tab \u2192"))
}

# UI Definition

ui <- page_sidebar(title = div(style = "display:flex; justify-content:space-between; align-items:center; width:100%; padding:18px 10px; border-bottom:1px solid rgba(255,255,255,0.10);",
                               
                               # Left side
                               div(h2("Exploration of Epigenomic Data in Colorectal Cancer",style = "margin:0; font-weight:600; font-size:28px; line-height:1.2;"),
                                   tags$div("Institut Germans Trias i Pujol (IGTP) · Universitat Politècnica de Catalunya (UPC)", style = "font-size:13px; color:#9fb3c8; margin-top:4px;"))),
                   
                   theme = bslib::bs_add_rules(bs_theme(version = 5, base_font = font_google("IBM Plex Sans"), heading_font = font_google("Libre Franklin"),
                                                        bg = "#f4f7f9", fg = "#0b2436", primary = "#0e7c86", secondary = "#16324f",success  = "#2fae66", info = "#3aa9c9", warning = "#e0a339", danger = "#d1495b", "navbar-bg"  = "#0b2436", base_font_size_scale = 0.98), c(gb_responsive_css, home_css)),
                   
                   sidebar = sidebar(width = 300,
                                     # Landing message shown only on the Home tab, in place of the data filters
                                     # (which don't apply until the person has picked an exploration tab).
                                     conditionalPanel("input.main_nav == 'Home'",
                                                      div(style = "text-align:center; padding:14px 4px;",
                                                          tags$p(style = "font-size:13px; color:#0b2436; margin-top:12px; line-height:1.6;",
                                                                 "Pick a tab above to start exploring the chromosome, bin, and view-level controls will appear here."))),
                                     conditionalPanel("input.main_nav != 'Home'",
                                                      tags$div("DATA SELECTION", style = "font-size:14px; font-weight:700; letter-spacing:0.8px; color:#7d92a3; margin-bottom:0px;"),
                                                      selectInput("chr", "Chromosome:", choices = chrom_list), selectizeInput("bins", "Selected bins:", choices = NULL, multiple = TRUE, options = list(placeholder = "Find and select bins (ex: 1_1000000)...", plugins = list("remove_button"), maxOptions = length(bin_table$bin_id) + 100)),
                                                      actionButton("add_chr_bins", "Add all bins on the chromosome", icon = bs_icon("plus-circle"), class = "btn-sm class=btn-outline-light w-100"),
                                                      actionButton("clear_bins", "Clear bin selection", icon = bs_icon("x-circle"), class = "btn-sm class=btn-outline-light w-100"),
                                                      fileInput("bins_file", "Or upload bins to select:", accept = c(".csv", ".tsv", ".txt", "text/csv", "text/tab-separated-values", "text/plain"), placeholder = "No file selected", buttonLabel = "Browse..."),
                                                      tags$div("Accepts .csv/.tsv/.txt or a plain list of bin IDs like 1_1000000, one per line or comma-separated).", style = "font-size:11px; color:#7d92a3; margin-top:-10px; margin-bottom:8px;"),
                                                      hr(style = "margin:1px 0; border-color:#3a5470;"),
                                                      tags$div("VISUALIZATION LEVEL", style = "font-size:14px; font-weight:700; letter-spacing:0.8px; color:#7d92a3; margin-bottom:0px;"),
                                                      radioButtons("view_level", NULL,
                                                                   choices = c("By Bin (aggregated)" = "bins", "By Sample (raw)" = "samples"),
                                                                   selected = "bins", inline = TRUE),
                                                      tags$div("Applies to plots where both views make sense: switches between pre-aggregated per-bin means and unaggregated per-sample methylation values.",
                                                               style = "font-size:11px; color:#7d92a3; margin-top:-10px; margin-bottom:8px;"),
                                                      hr(style = "margin:1px 0; border-color:#3a5470;"))),
                   
                   navset_card_underline(title = "Exploration and Visualization", id = "main_nav",
                                         # 0. Home
                                         nav_panel(title = "Home", value = "Home",
                                                   div(class = "home-page",
                                                       
                                                       # Hero
                                                       div(class = "home-hero",
                                                           div(class = "home-hero-text",
                                                               tags$div(class = "home-eyebrow", "Research Platform \u00b7 IGTP & UPC"),
                                                               tags$h1("Mapping the Epigenomic Landscape of Colorectal Cancer"),
                                                               tags$p(class = "home-hero-lede",
                                                                      "An interactive platform for exploring genome-wide DNA methylation differences between tumor and matched normal colorectal tissue, and for relating those differences to each patient's mutation status and clinical profile."),
                                                               div(class = "home-hero-actions",
                                                                   tags$button(type = "button", class = "home-cta-btn home-cta-primary",
                                                                               onclick = "Shiny.setInputValue('home_nav_target', 'Genome Browser', {priority:'event'});",
                                                                               "Launch Genome Browser"),
                                                                   tags$button(type = "button", class = "home-cta-btn home-cta-ghost",
                                                                               onclick = "Shiny.setInputValue('home_nav_target', 'Overview', {priority:'event'});",
                                                                               "See Data Overview"))),
                                                           div(class = "home-hero-stats",
                                                               home_stat(textOutput("n_patients", inline = TRUE), "Patients"),
                                                               home_stat(textOutput("n_samples", inline = TRUE), "Samples (Tumor + Normal)"),
                                                               home_stat(format(nrow(bin_table), big.mark = ","), "Genomic Bins (1 Mb)"),
                                                               home_stat(length(chrom_list), "Chromosomes Covered"))),
                                                       
                                                       # Objective + Data sources ----
                                                       div(class = "home-section home-grid-2",
                                                           div(class = "home-card",
                                                               h3("Project Objective"),
                                                               p("This platform was built to make genome-wide DNA methylation differences between colorectal tumor and matched normal tissue easy to explore, without writing code. The genome is split into 1\u00a0Mb bins, and methylation is compared bin-by-bin between tumor and normal samples across every chromosome."),
                                                               p("Beyond the tumor/normal contrast, the platform lets you relate methylation shifts to each patient's clinical and mutation profile, sex, stage, MSI/MSS status, and BRAF/KRAS/TP53 mutation status \u2014 and to overlay genomic context such as CpG islands, Alu repeats, and COSMIC cancer-gene annotation, in order to help surface candidate regions for further study.")),
                                                           div(class = "home-card",
                                                               h3("Data Sources"),
                                                               tags$ul(class = "home-source-list",
                                                                       home_source_item("Clinical & Mutation Metadata",
                                                                                        "One record per sample: patient ID, sample type (Tumor/Normal), sex, age, tumor stage, MSI/MSS status, relapse, and BRAF/KRAS/TP53 mutation status."),
                                                                       home_source_item("Genomic Bin Table",
                                                                                        "The genome split into 1\u00a0Mb bins, each with mean/SD methylation (Tumor and Normal), CpG and Alu content, and COSMIC Cancer Gene Census + promoter annotation."),
                                                                       home_source_item("Per-sample Methylation Values",
                                                                                        "Raw, unaggregated methylation measurements per sample and bin, used for sample-level views such as density plots, PCA/UMAP, and the similarity network.")))),
                                                       
                                                       # Key functionalities
                                                       div(class = "home-section",
                                                           div(class = "home-section-title", "Key Functionalities"),
                                                           p(class = "home-section-lede", "Every card below opens the matching tab to start exploring."),
                                                           div(class = "home-func-grid",
                                                               home_func_card("Genome Browser", "Zoom into any chromosome, search by coordinate, bin ID, or gene, and compare Tumor vs Normal methylation base-pair by base-pair.", "Genome Browser"),
                                                               home_func_card("Genome-wide Profile", "A Manhattan-style view across all chromosomes flags statistically significant and outlier bins in one glance.", "Genome-wide Profile"),
                                                               home_func_card("Feature \u00d7 Chromosome Heatmap", "Cluster patients and bins by methylation shift, and overlay a clinical feature as a colour strip.", "Feature \u00d7 Chromosome Heatmap"),
                                                               home_func_card("Overview", "A summary dashboard of the Tumor vs Normal shift, split by mutation status, stage/MSI, and sex.", "Overview"),
                                                               home_func_card("Clinical Explorer", "Cross-reference methylation with each patient's full clinical and mutation profile, filterable by subgroup.", "Clinical Explorer"),
                                                               home_func_card("Bin Table", "The full annotated bin-level table \u2014 methylation stats, CpG/Alu content, COSMIC genes and promoters \u2014 exportable to CSV.", "Bin Table"),
                                                               home_func_card("Tumor vs. Normal", "Density plots, PCA/UMAP projections, and a patient similarity network reveal global Tumor/Normal separation.", "Tumor vs. Normal"))),
                                                       
                                                       div(class = "home-footer",
                                                           "Exploration of Epigenomic Data in Colorectal Cancer \u00b7 Institut Germans Trias i Pujol (IGTP) & Universitat Polit\u00e8cnica de Catalunya (UPC) \u00b7 Built with R Shiny, bslib, plotly, ggplot2 and patchwork."))),
                                         
                                         # 1. Genome Browser
                                         nav_panel("Genome Browser",
                                                   div(style = "font-size:13px; color:#6c757d; margin-bottom:10px;",
                                                       "Browse methylation across the whole genome: pick a chromosome, zoom in/out, or search by coordinates, bin ID, or gene name."),
                                                   card(full_screen = TRUE,
                                                        card_header("All chromosomes, click a point to jump to that chromosome"),
                                                        div(class = "gb-plot-wrap", plotlyOutput("genome_overview_plot", height = "calc(max(320px, min(40vh, 420px)))", width = "100%")),
                                                        plotly_export_ui("genome_overview_plot")),
                                                   card(full_screen = TRUE,
                                                        card_header(textOutput("browser_position_header", inline = TRUE)),
                                                        div(class = "gb-toolbar",
                                                            div(class = "gb-toolbar-group gb-search-group",
                                                                textInput("browser_search", label = NULL, placeholder = "e.g. 12:25000000-26000000, 12_25000000, or KRAS", width = "100%"),
                                                                actionButton("browser_search_go", NULL, icon = bs_icon("search"), class = "btn-sm btn-outline-secondary")),
                                                            div(class = "gb-toolbar-group gb-nav-group",
                                                                actionButton("prev_chr", NULL, icon = bs_icon("chevron-left"), class = "btn-sm btn-outline-secondary"),
                                                                selectInput("browser_chr", NULL, choices = chrom_list, selected = "1", width = "90px"),
                                                                actionButton("next_chr", NULL, icon = bs_icon("chevron-right"), class = "btn-sm btn-outline-secondary")),
                                                            div(class = "gb-toolbar-group gb-zoom-group",
                                                                actionButton("zoom_in", NULL, icon = bs_icon("zoom-in"), class = "btn-sm btn-outline-secondary"),
                                                                actionButton("zoom_out", NULL, icon = bs_icon("zoom-out"), class = "btn-sm btn-outline-secondary"),
                                                                actionButton("zoom_reset", "Whole chromosome", icon = bs_icon("arrow-counterclockwise"), class = "btn-sm btn-outline-secondary"))),
                                                        div(class = "gb-plot-wrap", plotlyOutput("browser_chr_plot", height = "calc(max(440px, min(72vh, 760px)))", width = "100%")),
                                                        plotly_export_ui("browser_chr_plot"),
                                                        plot_desc("Tumor vs Normal methylation line-by-line across one chromosome."),
                                                        tags$div("Tip: bins currently in your sidebar selection are outlined on the plot above. Drag on the plot to zoom in, double-click to reset, or use the range slider below the plot to pan across the chromosome. Use the expand icon in the card header to view the plot full-screen.",
                                                                 style = "font-size:11px; color:#7d92a3; margin-top:6px;")),
                                                   card(card_header("External Resources for This Region"), uiOutput("browser_external_links"))),
                                         
                                         # 2. Genome-wide Profile
                                         nav_panel("Genome-wide Profile", div(style = "font-size:13px; color:#6c757d; margin-bottom:10px;","Genome-wide Tumor vs. Normal methylation across 1 Mb bins. Significant bins (q < 0.05) and top 5% extreme differences are flagged. Hover for details; drag to zoom, double-click to reset, and use the camera icon to export."),
                                                   card(card_header("Manhattan-style plot: Tumor \u2212 Normal methylation difference"), plotlyOutput("manhattan_plot", height = "calc(max(360px, min(58vh, 620px)))"), plotly_export_ui("manhattan_plot"), plot_desc("Every bin's Tumor\u2212Normal difference across the whole genome, colours mark statistically significant and outlier bins."),
                                                        div(style = "display:flex; gap:20px; flex-wrap:wrap; margin-top:12px; padding-top:10px; border-top:1px solid #e3e8ec; font-size:12px; color:#3a5470;",
                                                            tags$span(tags$span(style = "display:inline-block; width:10px; height:10px; border-radius:50%; background:#8D99AE; margin-right:6px;"), "Significant hypermethylation (q < 0.05)"),
                                                            tags$span(tags$span(style = "display:inline-block; width:10px; height:10px; border-radius:50%; background:#76C893; margin-right:6px;"), "Significant hypomethylation (q < 0.05)"),
                                                            tags$span(tags$span(style = "display:inline-block; width:10px; height:10px; border-radius:50%; background:#F4A261; margin-right:6px;"), "Extreme outlier (top 5% |\u0394|, not significant)"),
                                                            tags$span(tags$span(style = "display:inline-block; width:10px; height:10px; border-radius:50%; background:#168AAD; margin-right:6px;"), "Not significant")))),
                                         
                                         # 3. Feature x Chromosome Heatmap
                                         nav_panel("Feature × Chromosome Heatmap", card(full_screen = TRUE, card_header("Diverging heatmap of methylation shift (feature × chromosome)"), div(style = "font-size:13px; color:#6c757d; margin-bottom:10px;", "Each row is a patient, each column a bin, and the colour is the Tumor \u2212 Normal methylation shift for that patient/bin. Both axes are ordered by hierarchical clustering, dendrograms along the top and left, so patients and bins with similar methylation-shift patterns sit next to each other. Select a feature below to add a colour strip showing that clinical variable alongside the heatmap (this only changes the colour strip, not the row order). Use the zoom slider to zoom in for detail or out to see the complete graph at once, click \"Fit to screen\" to size it automatically, or use the expand icon (top-right of this card) for a full-screen view."),
                                                                                        div(style = "display:flex; flex-wrap:wrap; gap:28px; align-items:flex-end; margin-bottom:12px;",
                                                                                            div(style = "flex:1 1 220px; max-width:340px;", selectInput("heatmap_feature", "Select Feature:", choices = heatmap_feature_choices, selected = "none")),
                                                                                            div(style = "flex:1 1 220px; min-width:220px; max-width:300px;", sliderInput("heatmap_zoom", "Zoom:", min = 10, max = 400, value = 100, step = 5, post = "%", width = "100%")),
                                                                                            div(style = "display:flex; align-items:center; gap:8px; padding-left:20px; padding-bottom:2px; border-left:1px solid #e3e8ec;",
                                                                                                actionButton("heatmap_zoom_fit", "Fit to screen", icon = bs_icon("aspect-ratio"), class = "btn-sm btn-outline-secondary"),
                                                                                                actionButton("heatmap_zoom_reset", "Reset (100%)", icon = bs_icon("arrow-counterclockwise"), class = "btn-sm btn-outline-secondary"))),
                                                                                        tags$div(id = "heatmap_scroll_container", style = "overflow:auto; width:100%; height:92vh; border:1px solid #e3e8ec; border-radius:8px; background:#ffffff;",
                                                                                                 plotOutput("feature_heatmap")),
                                                                                        plot_export_ui("feature_heatmap"),
                                                                                        tags$script(HTML("
                                                                                          $(document).on('click', '#heatmap_zoom_fit', function() {
                                                                                            var el = document.getElementById('heatmap_scroll_container');
                                                                                            if (el) {
                                                                                              Shiny.setInputValue('heatmap_container_dims', {w: el.clientWidth, h: el.clientHeight, nonce: Math.random()}, {priority: 'event'});
                                                                                            }
                                                                                          });
                                                                                        ")),
                                                                                        plot_desc("Red/blue tiles show methylation shifts by patient and bin; both axes are ordered by dendrogram clustering on similarity, and the optional strip shows the chosen clinical feature. Zoom out (or click \"Fit to screen\") to see the whole heatmap without scrolling; zoom in to read fine detail, then scroll or drag inside the plot area to pan around."))),
                                         
                                         # 4. Overview
                                         nav_panel("Overview", uiOutput("overview_content")),
                                         
                                         # 5. Clinical Explorer
                                         nav_panel("Clinical Explorer", div(style = "font-size:13px; color:#6c757d; margin-bottom:10px;", "Compare tumor methylation by mutation status, and cross-reference patients against their full clinical profile. Filter the patient set below to restrict both the boxplot and the table to a specific subgroup. Click a row in the table to highlight that patient in the plot; use the search box or column filters to narrow the table down further."),
                                                   uiOutput("selected_patient_banner"), layout_columns(col_widths = breakpoints(sm = 12, md = 12, lg = c(5, 7)), card(card_header("Methylation by mutation status"),
                                                                                                                                                                      div(style = "display:flex; flex-wrap:wrap; gap:18px; align-items:flex-end; margin-bottom:6px;",
                                                                                                                                                                          div(style = "min-width:160px;", selectInput("mutation_gene", "Gene:", choices = c("KRAS", "BRAF", "TP53"))),
                                                                                                                                                                          div(style = "min-width:220px;", radioButtons("mutation_stat_test", "Statistical test:", choices = c("Wilcoxon rank-sum" = "wilcox", "Welch's t-test" = "ttest"), selected = "wilcox", inline = TRUE))),
                                                                                                                                                                      plotOutput("clinical_boxplot", height = "calc(max(300px, min(46vh, 460px)))"), plot_export_ui("clinical_boxplot"), plot_desc("Compares mutant vs wild-type methylation, with a p-value (Wilcoxon or Welch's t-test), among patients currently matching the filters above.")),
                                                                                                       card(card_header(div(style = "display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:8px;",
                                                                                                                            div(style = "display:flex; align-items:center; gap:8px;", bs_icon("table", size = "1.05em"), "Clinical metadata explorer"),
                                                                                                                            uiOutput("clinical_table_count_badge"))),
                                                                                                            div(style = "display:flex; flex-wrap:wrap; gap:10px 16px; align-items:center; background:#f7f9fa; border:1px solid #e7ecef; border-radius:8px; padding:8px 12px; margin-bottom:10px;",
                                                                                                                tags$span(style = "font-size:11px; font-weight:700; text-transform:uppercase; letter-spacing:0.4px; color:#7d92a3; display:flex; align-items:center; gap:4px;", bs_icon("sliders", size = "0.9em"), "Filter:"),
                                                                                                                build_clinical_categorical_filter_inputs(clinical_categorical_cols, clinical_profile_table, compact = TRUE),
                                                                                                                build_clinical_numeric_filter_inputs(clinical_numeric_cols, clinical_profile_table, compact = TRUE),
                                                                                                                actionButton("clinfilter_reset", "Reset", icon = bs_icon("arrow-counterclockwise"), class = "btn-sm btn-outline-danger", style = "font-size:11px; padding:3px 10px; border-radius:5px; margin-left:auto;")),
                                                                                                            div(style = "font-size:12px; color:#7d92a3; margin-bottom:8px;", "One row per patient. Use the column filters below the headers to narrow results, click any column header to sort, and click a row to select that patient."),
                                                                                                            DTOutput("clinical_table"), plot_desc("A searchable table of patient info, click a row to highlight them.")))),
                                         
                                         # 6. Bin Table
                                         nav_panel("Bin Table", card(card_header(div(style = "display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:8px;",
                                                                                     div(style = "display:flex; align-items:center; gap:8px;", bs_icon("table", size = "1.05em"), "Filterable Bin-level Data"),
                                                                                     uiOutput("bintable_count_badge"))),
                                                                     tags$details(open = "open", style = "background:#ffffff; border:1px solid #d8e0e6; border-radius:10px; padding:0; margin-bottom:20px; box-shadow:0 1px 4px rgba(22,50,79,0.08); overflow:hidden;",
                                                                                  tags$summary(style = "font-weight:700; font-size:14.5px; color:#ffffff; background:#16324f; cursor:pointer; display:flex; align-items:center; gap:8px; padding:12px 18px; list-style:none;", bs_icon("sliders", size = "1.1em"), "Filter Bin Metrics",
                                                                                               tags$span(style = "margin-left:auto; font-weight:400; font-size:11px; color:#c7d2da;", "Click to expand / collapse")), div(style = "padding:18px 20px 20px 20px;",
                                                                                                                                                                                                                          
                                                                                                                                                                                                                          # Chromosome selection
                                                                                                                                                                                                                          div(style = "margin-bottom:16px;",
                                                                                                                                                                                                                              div(style = "display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:8px; margin-bottom:8px;",
                                                                                                                                                                                                                                  tags$span(style = "font-weight:700; font-size:11.5px; text-transform:uppercase; letter-spacing:0.6px; color:#7d92a3; display:flex; align-items:center; gap:6px;",
                                                                                                                                                                                                                                            bs_icon("bar-chart-steps"), "Chromosome"),
                                                                                                                                                                                                                                  div(style = "display:flex; gap:6px;", actionButton("binfilter_chr_all", "Select all", class = "btn-sm btn-outline-secondary", style = "font-size:11px; padding:2px 10px; border-radius:5px;"), actionButton("binfilter_chr_clear", "Clear", class = "btn-sm btn-outline-secondary", style = "font-size:11px; padding:2px 10px; border-radius:5px;"))),
                                                                                                                                                                                                                              div(style = "background:#f7f9fa; border:1px solid #e7ecef; border-radius:8px; padding:10px 12px;", checkboxGroupInput("binfilter_chr", label = NULL, choices = chrom_list, selected = chrom_list, inline = TRUE))),
                                                                                                                                                                                                                          tags$hr(style = "border-top:1px solid #eef2f5; margin:14px 0;"),
                                                                                                                                                                                                                          
                                                                                                                                                                                                                          # Gene search
                                                                                                                                                                                                                          div(style = "margin-bottom:16px;",
                                                                                                                                                                                                                              div(style = "display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:8px; margin-bottom:8px;",
                                                                                                                                                                                                                                  tags$span(style = "font-weight:700; font-size:11.5px; text-transform:uppercase; letter-spacing:0.6px; color:#7d92a3; display:flex; align-items:center; gap:6px;", bs_icon("search"), "Gene"), actionButton("binfilter_gene_clear", "Clear", class = "btn-sm btn-outline-secondary", style = "font-size:11px; padding:2px 10px; border-radius:5px;")),
                                                                                                                                                                                                                              div(style = "background:#f7f9fa; border:1px solid #e7ecef; border-radius:8px; padding:10px 12px;", selectizeInput("binfilter_gene", label = NULL, choices = NULL, selected = character(0), options = list(placeholder = "Search a gene, e.g. TP53...", maxOptions = 200), width = "100%"),
                                                                                                                                                                                                                                  tags$div(style = "font-size:11px; color:#7d92a3; margin-top:4px;", "Shows every bin annotated with this exact gene symbol (COSMIC Cancer Gene Census annotation)."))),
                                                                                                                                                                                                                          tags$hr(style = "border-top:1px solid #eef2f5; margin:14px 0;"),
                                                                                                                                                                                                                          
                                                                                                                                                                                                                          # Numeric metric filters
                                                                                                                                                                                                                          tags$div(style = "font-weight:700; font-size:11.5px; text-transform:uppercase; letter-spacing:0.6px; color:#7d92a3; margin-bottom:10px; display:flex; align-items:center; gap:6px;", bs_icon("sliders2"), "Methylation & Content Metrics"),
                                                                                                                                                                                                                          div(style = "display:flex; flex-wrap:wrap; gap:16px; align-items:flex-start;",
                                                                                                                                                                                                                              build_bin_table_filter_inputs(bin_table_filters, bin_table)),
                                                                                                                                                                                                                          tags$hr(style = "border-top:1px solid #eef2f5; margin:18px 0 14px 0;"),
                                                                                                                                                                                                                          div(style = "display:flex; justify-content:flex-end;", actionButton("binfilter_reset", "Reset all filters", icon = bs_icon("arrow-counterclockwise"), class = "btn-sm btn-outline-danger", style = "font-size:12px; border-radius:6px; padding:5px 14px;"))) ),
                                                                     div(style = "display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:12px; margin-bottom:15px;", downloadButton( "bintable_download", "Download CSV Data", icon = icon("download"), class = "btn-lg btn-primary", style = "font-weight: 600; padding: 10px 22px; font-size: 15px; border-radius: 6px;"),
                                                                         uiOutput("bintable_gene_note"),
                                                                         uiOutput("bintable_promoter_note")),
                                                                     div(style = "margin-top:4px;", DTOutput("bintable")))),
                                         
                                         # 7. Tumor vs Normal
                                         nav_panel("Tumor vs. Normal", navset_tab(
                                           nav_panel("Methylation density", card(full_screen = TRUE, card_header("Methylation density (Tumor vs Normal)"), plotOutput("tn_density", height = "58vh"), plot_export_ui("tn_density"), plot_desc("Shape of the methylation distribution, Tumor vs Normal overlapped."))),
                                           nav_panel("PCA / UMAP", card(full_screen = TRUE, card_header("PCA / UMAP (interactive)"),
                                                                        div(style = "display:flex; flex-wrap:wrap; gap:18px; align-items:flex-end; margin-bottom:10px;",
                                                                            div(style = "min-width:160px;", radioButtons("proj_method", "Method:", choices = c("PCA", "UMAP"), inline = TRUE)),
                                                                            div(style = "min-width:150px;", selectInput("proj_xdim", "X axis:", choices = setNames(1:10, paste0("Component ", 1:10)), selected = 1)),
                                                                            div(style = "min-width:150px;", selectInput("proj_ydim", "Y axis:", choices = setNames(1:10, paste0("Component ", 1:10)), selected = 2)),
                                                                            div(style = "min-width:260px; flex:1 1 260px;", selectizeInput("proj_samples", "Samples (optional filter):", choices = NULL, multiple = TRUE, options = list(placeholder = "All samples (leave empty to include everyone)...", plugins = list("remove_button"))))),
                                                                        plotlyOutput("tn_projection", height = "54vh"),
                                                                        plotly_export_ui("tn_projection"),
                                                                        plot_desc("PCA finds the directions where samples vary the most, UMAP groups similar samples together. Pick any of the first 10 components/dimensions for each axis (PCA1\u2013PCA10 / UMAP1\u2013UMAP10), and optionally restrict the projection to specific samples."))),
                                           nav_panel("Patient similarity network", card(full_screen = TRUE, card_header("Patient similarity network"),
                                                                                        div(style = "display:flex; flex-wrap:wrap; gap:18px; align-items:flex-end; margin-bottom:10px;",
                                                                                            div(style = "min-width:260px; max-width:420px; flex:1 1 300px;",
                                                                                                selectizeInput("network_color_feature", "Colour nodes by (pick one or more):", choices = network_color_choices, selected = "Type", multiple = TRUE,
                                                                                                               options = list(plugins = list("remove_button")))),
                                                                                            div(style = "min-width:240px; max-width:340px; flex:1 1 260px;",
                                                                                                sliderInput("network_corr_threshold", "Minimum correlation for an edge (corr \u2265):", min = -1, max = 1, value = 0.7, step = 0.01, width = "100%"))),
                                                                                        plotOutput("tn_network", height = "54vh"),
                                                                                        plot_export_ui("tn_network"),
                                                                                        plot_desc("Each dot is one sample (Tumor or Normal, shown separately, as in the metadata table): a circle for Tumor and a triangle for Normal. Lines connect sample pairs whose methylation profiles correlate at or above the chosen threshold, so you can see at a glance whether the samples that resemble each other are Tumor or Normal. Pick more than one colour variable to compare side by side.")))))),
                   div(style = "display:flex; align-items:center; justify-content:center; gap:18px; margin-top:24px; padding:16px 0; border-top:1px solid #dde3e8; color:#6c757d; font-size:13px;", tags$span("Paula Artiz Dueñas, Bioinformatics Student", style = "margin-left:10px;")))

# Server Definition

server <- function(input, output, session) {
  
  # "Selected bins:" always offers only the currently chosen chromosome's bins
  bins_choices_for_chr <- function(chr, extra = character(0)) {
    union(bin_choices_by_chr[[chr]], extra)
  }
  
  updateSelectizeInput(session, "bins", choices = bins_choices_for_chr(isolate(input$chr)), selected = character(0), server = TRUE)
  updateSelectizeInput(session, "binfilter_gene", choices = all_gene_names, selected = character(0), server = TRUE)
  updateSelectizeInput(session, "proj_samples", choices = sort(unique(metadata$sample_id)), selected = character(0), server = TRUE)
  
  # picking a different chromosome re-scopes the "Selected bins" search to that chromosome
  observeEvent(input$chr, {
    updateSelectizeInput(session, "bins", choices = bins_choices_for_chr(input$chr, input$bins), selected = input$bins, server = TRUE)
  }, ignoreInit = TRUE)
  
  observeEvent(input$add_chr_bins, {
    chr_bins <- bin_choices_by_chr[[input$chr]]
    updated <- union(input$bins, chr_bins)
    updateSelectizeInput(session, "bins", choices = bins_choices_for_chr(input$chr, updated), selected = updated, server = TRUE)
  })
  
  observeEvent(input$clear_bins, {
    updateSelectizeInput(session, "bins", choices = bins_choices_for_chr(input$chr), selected = character(0), server = TRUE)
  })
  
  observeEvent(input$bins_file, {
    req(input$bins_file)
    res <- parse_bin_file(input$bins_file$datapath, input$bins_file$name, bin_table, chrom_list)
    if (length(res$matched) == 0) {
      showNotification("Couldn't match any bins in file.", type = "error", duration = 7)
      return()
    }
    updated <- union(input$bins, res$matched)
    updateSelectizeInput(session, "bins", choices = bins_choices_for_chr(input$chr, updated), selected = updated, server = TRUE)
    n_unmatched <- length(res$unmatched)
    if (n_unmatched == 0) {
      showNotification(
        sprintf("All %d bin IDs in the file matched and were added to your selection.", length(res$matched)),
        type = "message", duration = 7)
    } else {
      preview <- paste(utils::head(res$unmatched, 10), collapse = ", ")
      more <- if (n_unmatched > 10) sprintf(" (+%d more)", n_unmatched - 10) else ""
      showNotification(
        sprintf("Selected %d matching bins. %d ID(s) from the file don't exist in this dataset and were skipped: %s%s",
                length(res$matched), n_unmatched, preview, more),
        type = "warning", duration = 15)
    }
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
    span(style = "font-size:12px; font-weight:600; color:#16324f; background:#eef2f5; border:1px solid #d8e0e6; border-radius:20px; padding:4px 12px; white-space:nowrap;",
         paste0(format(n, big.mark = ","), " of ", format(total, big.mark = ","), " bins"))
  })
  
  output$n_samples <- renderText({ nrow(metadata) })
  output$n_patients <- renderText({ length(unique(metadata$patient_id)) })
  
  # 0. Home: clicking a hero button or a "Key Functionalities" card sets home_nav_target
  observeEvent(input$home_nav_target, {
    nav_select("main_nav", selected = input$home_nav_target, session = session)
  })
  
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
  
  # small label reused on every Overview / Tumor vs. Normal plot so it's clear which region is shown
  region_label <- reactive({
    df <- selected_bin_table()
    n <- nrow(df)
    chrs <- unique(as.character(df$chr))
    chr_txt <- if (length(chrs) > 0) paste0(" (chr", paste(chrs, collapse = ", chr"), ")") else ""
    paste0("Selected ", n, " bin", if (n != 1) "s" else "", chr_txt)
  })
  
  overview_plot_meta <- list(
    list(id = "overview_boxplot", title = "Tumor vs. Normal Methylation", desc = "Shows if tumors are more or less methylated than normals."),
    list(id = "overview_mutations", title = "Mutation Status vs. Methylation Shift", desc = "Shows if having a KRAS/BRAF/TP53 mutation changes methylation shift."),
    list(id = "overview_composition", title = "Stage & MSI/MSS vs. Methylation Shift", desc = "Shows if shift differs by stage or MSS status."),
    list(id = "overview_sex_comparison", title = "Sex vs. Methylation Shift", desc = "Shows if methylation shift differs between males and females (relevant on chrX/chrY)."))
  
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
        card(card_header(
          div(style = "display:flex; justify-content:space-between; align-items:center; width:100%;",
              tags$span(pm$title),
              actionButton(paste0("expand_", pm$id), NULL, icon = bs_icon("arrows-fullscreen"),class = "btn-sm btn-outline-secondary", title = "Expand", style = "padding:2px 7px; margin-left:auto; margin-right:-5px;"))),
          plotOutput(pm$id, height = "calc(max(230px, min(34vh, 320px)))", click = paste0(pm$id, "_click")),
          plot_export_ui(pm$id),
          plot_desc(pm$desc)) # Show graphs descriptions
      })
      # stacks to a single column on phones/tablets (< md), 2x2 grid from tablet-landscape/desktop (md) up
      do.call(layout_columns, c(list(col_widths = breakpoints(sm = 12, md = 6), row_heights = "auto"), card_list))
      
    } else {
      # enlarged single-plot view, with a button to go back to the 2x2 grid
      meta <- Filter(function(pm) identical(pm$id, zoom), overview_plot_meta)[[1]]
      tagList(
        div(style = "margin-bottom:10px;",
            actionButton("overview_back", "Go back to Overview", icon = bs_icon("arrow-left-circle"), class = "btn-sm btn-outline-primary")),
        card(card_header(meta$title), plotOutput(meta$id, height = "calc(max(360px, min(64vh, 640px)))"), plot_export_ui(meta$id)))
    }
  })
  
  overview_sex_comparison_build <- reactive({
    # Counted at the patient level (via overview_patient_shift), not per bin-row
    df <- overview_patient_shift()
    validate(need(!is.null(df) && nrow(df) > 0, "No data available for the current selection."))
    
    # normalize case/whitespace before the lookup so stray formatting in the source data
    sex_labels <- c("DONA" = "Female", "HOME" = "Male")
    df$Sex <- sex_labels[toupper(trimws(as.character(df$sexe)))]
    validate(need(any(!is.na(df$Sex)), "No sex information available for the current selection."))
    df <- df[!is.na(df$Sex), ]
    validate(need(any(is.finite(df$shift)), "Methylation shift could not be computed for the current selection (missing Tumor or Normal values)."))
    df <- df[is.finite(df$shift), ]
    df$Sex <- factor(df$Sex, levels = intersect(c("Female", "Male"), unique(df$Sex)))
    
    n_by_sex <- table(df$Sex)
    df$SexLabel <- factor(paste0(as.character(df$Sex), " (n=", n_by_sex[as.character(df$Sex)], ")"),
                          levels = paste0(names(n_by_sex), " (n=", n_by_sex, ")"))
    
    p <- ggplot(df, aes(x = SexLabel, y = shift, fill = Sex)) +
      geom_hline(yintercept = 0, color = "#7d92a3", linewidth = 0.4, linetype = "dashed") +
      geom_jitter(width = 0.08, size = 0.9, alpha = 0.3, color = "#16324f") +
      geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.92, color = "#16324f", linewidth = 0.4) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = "white", color = "#16324f", stroke = 0.8) +
      scale_fill_manual(values = c("Female" = "#d1495b", "Male" = "#3aa9c9")) +
      labs(x = NULL, y = "Methylation shift (Tumor \u2212 Normal)", subtitle = region_label(), caption = "\u25c7 mean") +
      theme_app() +
      theme(legend.position = "none")
    
    # only annotate a test statistic when both sexes are actually present
    if (length(levels(df$Sex)) == 2) {
      test_res <- run_group_test(df$shift[df$Sex == "Female"], df$shift[df$Sex == "Male"], method = "wilcox")
      y_max <- max(df$shift, na.rm = TRUE)
      p <- p + annotate("text", x = 1.5, y = y_max * 1.12, label = test_res$label, size = 3.2, color = "#16324f", fontface = "bold") +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.14)))
    }
    list(plot = p, data = df[, c("patient_id", "Sex", "shift")])
  })
  output$overview_sex_comparison <- renderPlot({ overview_sex_comparison_build()$plot })
  register_plot_export(output, "overview_sex_comparison",
                       function() overview_sex_comparison_build()$plot,
                       function() overview_sex_comparison_build()$data)
  
  # TRUE when the user wants raw per-sample values instead of pre-aggregated per-bin means
  view_by_samples <- reactive({ identical(input$view_level, "samples") })
  
  overview_boxplot_build <- reactive({
    if (view_by_samples()) {
      df <- overview_ml_selected()
      validate(need(nrow(df) > 0, "No bins selected. Choose a chromosome or bins in the sidebar."))
      plot_df <- data.frame(Type = df$Type, Methylation = df$methylation)
      plot_df <- plot_df[is.finite(plot_df$Methylation) & plot_df$Type %in% c("Normal", "Tumor"), ]
      validate(need(nrow(plot_df) > 0, "No valid methylation values for the current selection."))
      y_lab <- "Methylation (per sample \u00d7 bin)"
    } else {
      df <- selected_bin_table()
      validate(need(nrow(df) > 0, "No bins selected. Choose a chromosome or bins in the sidebar."))
      plot_df <- data.frame( Type = rep(c("Normal", "Tumor"), each = nrow(df)), Methylation = c(df$mean_methylation_normal, df$mean_methylation_tumor))
      plot_df <- plot_df[is.finite(plot_df$Methylation), ]
      validate(need(nrow(plot_df) > 0, "No valid methylation values for the current selection."))
      y_lab <- "Mean methylation (per bin)"
    }
    plot_df$Type <- factor(plot_df$Type, levels = c("Normal", "Tumor"))
    
    p <- ggplot(plot_df, aes(x = Type, y = Methylation, fill = Type)) +
      geom_jitter(width = 0.07, size = 0.9, alpha = 0.25, color = "#16324f") +
      geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.92, color = "#16324f", linewidth = 0.4) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = "white", color = "#16324f", stroke = 0.8) +
      scale_fill_manual(values = c("Tumor" = "#d1495b", "Normal" = "#3aa9c9")) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0.06))) +
      labs(x = NULL, y = y_lab, subtitle = region_label(), caption = "\u25c7 mean") +
      theme_app() +
      theme(legend.position = "none")
    list(plot = p, data = plot_df)
  })
  output$overview_boxplot <- renderPlot({ overview_boxplot_build()$plot })
  register_plot_export(output, "overview_boxplot",
                       function() overview_boxplot_build()$plot,
                       function() overview_boxplot_build()$data)
  
  overview_mutations_build <- reactive({
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
      levels = c("Wild-type", "Mutant", "Unknown"))
    n_by_gene <- table(long_df$Gene)
    long_df$GeneLabel <- factor(
      paste0(long_df$Gene, " (n=", n_by_gene[long_df$Gene], ")"),
      levels = paste0(genes, " (n=", n_by_gene[genes], ")"))
    
    p <- ggplot(long_df, aes(x = Status, y = Shift, fill = Status)) +
      geom_hline(yintercept = 0, color = "#7d92a3", linewidth = 0.4, linetype = "dashed") +
      geom_jitter(width = 0.1, size = 0.9, alpha = 0.3, color = "#16324f") +
      geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.92, color = "#16324f", linewidth = 0.4) +
      scale_fill_manual(values = c("Wild-type" = "#3aa9c9", "Mutant" = "#d1495b", "Unknown" = "darkgray")) +
      facet_wrap(~ GeneLabel, nrow = 1, scales = "free_x") +
      labs(x = NULL, y = "Methylation shift (Tumor \u2212 Normal)", subtitle = region_label()) +
      theme_app() + # Predefined before
      theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
    list(plot = p, data = long_df[, c("Gene", "Status", "Shift")])
  })
  output$overview_mutations <- renderPlot({ overview_mutations_build()$plot })
  register_plot_export(output, "overview_mutations",
                       function() overview_mutations_build()$plot,
                       function() overview_mutations_build()$data)
  
  overview_composition_build <- reactive({
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
    
    p <- ggplot(agg, aes(x = Level, y = Shift, fill = Feature)) +
      geom_hline(yintercept = 0, color = "#7d92a3", linewidth = 0.4) +
      geom_col(width = 0.6, show.legend = FALSE) +
      geom_text(aes(y = Shift, label = paste0("n=", n), vjust = ifelse(Shift >= 0, -0.5, 1.3)), size = 2.7, color = "#16324f") +
      scale_fill_manual(values = c("estadi2" = "#d1495b", "MSS" = "#3aa9c9")) +
      scale_y_continuous(expand = expansion(mult = c(0.15, 0.18))) +
      facet_wrap(~ FeatureLabel, scales = "free_x", nrow = 1) +
      labs(x = NULL, y = "Mean methylation shift (Tumor \u2212 Normal)", subtitle = region_label()) +
      theme_app()
    list(plot = p, data = agg[, c("Feature", "FeatureLabel", "Level", "Shift", "n")])
  })
  output$overview_composition <- renderPlot({ overview_composition_build()$plot })
  register_plot_export(output, "overview_composition",
                       function() overview_composition_build()$plot,
                       function() overview_composition_build()$data)
  
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
    # keep the sidebar "Chromosome:" selector 
    if (!identical(input$chr, input$browser_chr)) {
      updateSelectInput(session, "chr", selected = input$browser_chr)
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
      nav$chr<- res$chr # set this first so the browser_chr observer above sees it and skips its reset
      nav$start <- max(1, floor(res$start))
      nav$end <- min(chrom_lengths[[res$chr]], ceiling(res$end))
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
    paste0("Genome Browser - chr", nav$chr, ": ",
           format(round(nav$start), big.mark = ",", scientific = FALSE), "-",
           format(round(nav$end), big.mark = ",", scientific = FALSE),
           " (", format_mb_label(nav$end - nav$start), ")")
  })
  
  # quick access to UCSC / Ensembl / TCGA for whatever region is currently shown in the genome browser
  output$browser_external_links <- renderUI({
    req(nav$chr, nav$start, nav$end)
    region_links <- build_external_links(nav$chr, nav$start, nav$end)
    
    div(style = "display:flex; gap:16px; flex-wrap:wrap; font-size:13px;",
        tags$a(bs_icon("box-arrow-up-right"), " UCSC Genome Browser", href = region_links$UCSC, target = "_blank", rel = "noopener noreferrer"),
        tags$a(bs_icon("box-arrow-up-right"), " Ensembl", href = region_links$Ensembl, target = "_blank", rel = "noopener noreferrer"),
        tags$a(bs_icon("box-arrow-up-right"), " TCGA", href = region_links$TCGA, target = "_blank", rel = "noopener noreferrer")
    )
  })
  
  browser_view_bins <- reactive({
    req(nav$chr, nav$start, nav$end)
    df <- bin_table[bin_table$chr == nav$chr & bin_table$bin_end >= nav$start & bin_table$bin_start <= nav$end, ]
    df[order(df$bin_position), ]
  })
  
  output$genome_overview_plot <- renderPlotly({
    df <- genome_wide_bins
    p <- plot_ly(df, x = ~genome_x, y = ~overall_meth, color = ~chr_parity, colors = c("#16324f", "#0e7c86"), customdata = ~as.character(chr),
                 type = "scatter", mode = "markers", marker = list(size = 4),text = ~bin_id, hoverinfo = "text", source = "genome_overview") |>
      layout(showlegend = FALSE, xaxis = list(title = "", tickvals = as.numeric(chr_mid), ticktext = chrom_list, tickangle = 45), yaxis = list(title = "Mean methylation", range = c(0, 1), automargin = TRUE), shapes = chr_boundary_shapes, margin = list(t = 20, r = 20, b = 40, l = 50), autosize = TRUE) |>
      config(displayModeBar = FALSE, responsive = TRUE) |>
      event_register("plotly_click")
  })
  register_plotly_export(output, "genome_overview_plot",
                         function() genome_wide_bins[, c("bin_id", "chr", "bin_position", "genome_x", "overall_meth")])
  
  # per-sample-per-bin methylation values within the current genome browser window
  browser_view_samples <- reactive({
    bins_now <- browser_view_bins()
    req(nrow(bins_now) > 0)
    df <- ml_annot[ml_annot$bin_id %in% bins_now$bin_id, c("sample_id", "patient_id", "Type", "bin_id", "bin_position", "methylation")]
    df[df$Type %in% c("Tumor", "Normal"), ]
  })
  
  output$browser_chr_plot <- renderPlotly({
    if (view_by_samples()) {
      df <- browser_view_samples()
      if (nrow(df) == 0) {
        return(
          plotly_empty(type = "scatter", mode = "markers") |>
            layout(title = list(text = "No bins in this region", font = list(size = 14))) |>
            config(displayModeBar = FALSE, responsive = TRUE))
      }
      # highlight only bins explicitly picked in the sidebar's "Selected bins" list -- deliberately
      # NOT selected_bins()'s whole-chromosome fallback, so an empty pick shows no highlighting here
      hl <- df$bin_id %in% input$bins
      base_size <- max(3, min(8, 4000 / nrow(df)))
      hl_size <- base_size + 4
      p <- plot_ly(df, x = ~bin_position / 1e6)
      for (ty in c("Tumor", "Normal")) {
        sub <- df[df$Type == ty, ]
        hl_sub <- sub$bin_id %in% input$bins
        p <- add_trace(p, data = sub, x = ~bin_position / 1e6, y = ~methylation, type = "scatter", mode = "markers", name = ty,
                       marker = list(color = if (ty == "Tumor") "#d1495b" else "#3aa9c9",
                                     size = ifelse(hl_sub, hl_size, base_size),
                                     opacity = 0.55,
                                     line = list(color = "#0b2436", width = ifelse(hl_sub, 1.2, 0))),
                       text = ~paste0(bin_id, "<br>Sample: ", sample_id, "<br>Patient: ", patient_id, "<br>", ty, ": ", round(methylation, 3)),
                       hoverinfo = "text")
      }
      p |>
        layout(xaxis = list(title = paste0("Position on chr", nav$chr, " (Mb)"), automargin = TRUE,
                            rangeslider = list(visible = TRUE, thickness = 0.08, bgcolor = "#eef2f5", bordercolor = "#d8e0e6", borderwidth = 1)),
               yaxis = list(title = "Methylation (per sample)", range = c(0, 1), automargin = TRUE),
               legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1, yanchor = "bottom"),
               hovermode = "closest",
               margin = list(t = 50, r = 20, b = 50, l = 55),
               autosize = TRUE) |>
        config(displayModeBar = FALSE, responsive = TRUE)
    } else {
      df <- browser_view_bins()
      if (nrow(df) == 0) {
        return(
          plotly_empty(type = "scatter", mode = "markers") |>
            layout(title = list(text = "No bins in this region", font = list(size = 14))) |>
            config(displayModeBar = FALSE, responsive = TRUE))
      }
      # highlight only bins explicitly picked in the sidebar's "Selected bins" list
      hl <- df$bin_id %in% input$bins
      base_size <- max(4, min(9, 700 / nrow(df)))
      hl_size <- base_size + 5
      plot_ly(df, x = ~bin_position / 1e6) |>
        add_trace(y = ~mean_methylation_tumor, type = "scatter", mode = "lines+markers", name = "Tumor",line = list(color = "#d1495b", width = 1.6),
                  marker = list(color = "#d1495b", size = ifelse(hl, hl_size, base_size), line = list(color = "#0b2436", width = ifelse(hl, 1.5, 0))),
                  text = ~paste0(bin_id, "<br>Tumor mean: ", round(mean_methylation_tumor, 3)), hoverinfo = "text") |>
        add_trace(y = ~mean_methylation_normal, type = "scatter", mode = "lines+markers", name = "Normal",
                  line = list(color = "#3aa9c9", width = 1.6), marker = list(color = "#3aa9c9", size = ifelse(hl, hl_size, base_size),line = list(color = "#0b2436", width = ifelse(hl, 1.5, 0))),
                  text = ~paste0(bin_id, "<br>Normal mean: ", round(mean_methylation_normal, 3)), hoverinfo = "text") |>
        layout(xaxis = list(title = paste0("Position on chr", nav$chr, " (Mb)"), automargin = TRUE,
                            rangeslider = list(visible = TRUE, thickness = 0.08, bgcolor = "#eef2f5", bordercolor = "#d8e0e6", borderwidth = 1)),
               yaxis = list(title = "Mean methylation", range = c(0, 1), automargin = TRUE),
               legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1, yanchor = "bottom"),
               hovermode = "x unified",
               margin = list(t = 50, r = 20, b = 50, l = 55),
               autosize = TRUE) |>
        config(displayModeBar = FALSE, responsive = TRUE)
    }
  })
  register_plotly_export(output, "browser_chr_plot", function() {
    if (view_by_samples()) browser_view_samples() else browser_view_bins()[, c("bin_id", "chr", "bin_position", "mean_methylation_tumor", "mean_methylation_normal")]
  })
  
  # 3. Tumor vs. Normal
  
  # sample x bin matrix of raw methylation for the currently selected bins (and, optionally, selected samples), used by PCA/UMAP
  projection_matrix <- reactive({
    bins_now <- selected_bins()
    req(length(bins_now) > 0)
    df <- ml_annot[ml_annot$bin_id %in% bins_now, c("sample_id", "bin_id", "methylation")]
    samples_now <- input$proj_samples
    if (!is.null(samples_now) && length(samples_now) > 0) {
      df <- df[df$sample_id %in% samples_now, , drop = FALSE]
    }
    mat <- build_wide_matrix(df, "sample_id", "bin_id", "methylation")
    req_ok <- !is.null(mat) && ncol(mat) >= 2
    if (!req_ok) return(NULL)
    keep_cols <- apply(mat, 2, function(x) stats::sd(x) > 0)
    mat <- mat[, keep_cols, drop = FALSE]
    if (nrow(mat) < 3 || ncol(mat) < 2) return(NULL)
    mat
  })
  
  # sample x bin matrix of raw methylation (Tumor and Normal kept as separate rows) for the selected bins
  network_sample_matrix <- reactive({
    bins_now <- selected_bins()
    req(length(bins_now) > 0)
    df <- ml_annot[ml_annot$bin_id %in% bins_now, c("sample_id", "bin_id", "methylation")]
    mat <- build_wide_matrix(df, "sample_id", "bin_id", "methylation")
    if (is.null(mat)) return(NULL)
    keep_cols <- apply(mat, 2, function(x) stats::sd(x) > 0)
    mat <- mat[, keep_cols, drop = FALSE]
    if (nrow(mat) < 3 || ncol(mat) < 2) return(NULL)
    mat
  })
  
  tn_density_build <- reactive({
    if (view_by_samples()) {
      df <- ml_annot[ml_annot$bin_id %in% selected_bins() & ml_annot$Type %in% c("Tumor", "Normal"), c("Type", "methylation")]
      validate(need(nrow(df) > 0, "No bins selected. Choose a chromosome or bins in the sidebar."))
      plot_df <- data.frame(Type = df$Type, Methylation = df$methylation)
      plot_df <- plot_df[is.finite(plot_df$Methylation), ]
      validate(need(nrow(plot_df) > 0, "No valid methylation values for the current selection."))
      x_lab <- "Methylation (per sample \u00d7 bin)"
    } else {
      df <- selected_bin_table()
      validate(need(nrow(df) > 0, "No bins selected. Choose a chromosome or bins in the sidebar."))
      plot_df <- data.frame(Type = rep(c("Tumor", "Normal"), each = nrow(df)),
                            Methylation = c(df$mean_methylation_tumor, df$mean_methylation_normal))
      plot_df <- plot_df[is.finite(plot_df$Methylation), ]
      validate(need(nrow(plot_df) > 0, "No valid methylation values for the current selection."))
      x_lab <- "Mean methylation (per bin)"
    }
    
    p <- ggplot(plot_df, aes(x = Methylation, fill = Type, color = Type)) +
      geom_density(alpha = 0.35, linewidth = 0.9, na.rm = TRUE) +
      scale_fill_manual(values = c("Tumor" = "#d1495b", "Normal" = "#3aa9c9")) +
      scale_color_manual(values = c("Tumor" = "#d1495b", "Normal" = "#3aa9c9")) +
      labs(x = x_lab, y = "Density", title = region_label(), fill = NULL, color = NULL) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "top")
    list(plot = p, data = plot_df)
  })
  output$tn_density <- renderPlot({ tn_density_build()$plot })
  register_plot_export(output, "tn_density",
                       function() tn_density_build()$plot,
                       function() tn_density_build()$data)
  
  tn_projection_build <- reactive({
    method <- input$proj_method
    if (is.null(method)) method <- "PCA"
    mat <- projection_matrix()
    validate(need(!is.null(mat),"Not enough samples with complete data for the selected bins (and sample filter) to compute a projection."))
    
    xdim <- as.integer(input$proj_xdim)
    ydim <- as.integer(input$proj_ydim)
    if (is.na(xdim)) xdim <- 1
    if (is.na(ydim)) ydim <- 2
    
    if (identical(method, "PCA")) {
      pca <- prcomp(mat, center = TRUE, scale. = TRUE)
      n_comp <- min(10, ncol(pca$x))
      validate(need(n_comp >= 2, "Not enough bins selected to compute at least two principal components."))
      validate(need(xdim <= n_comp && ydim <= n_comp,
                    paste0("Only ", n_comp, " principal component(s) can be computed for the current selection. Pick a component up to PC", n_comp, ".")))
      var_exp <- round(100 * summary(pca)$importance[2, c(xdim, ydim)], 1)
      coords <- as.data.frame(pca$x[, c(xdim, ydim), drop = FALSE])
      names(coords) <- c("Dim1", "Dim2")
      xlab <- paste0("PC", xdim, " (", var_exp[1], "%)")
      ylab <- paste0("PC", ydim, " (", var_exp[2], "%)")
    } else {
      validate(need(requireNamespace("uwot", quietly = TRUE), "UMAP requires the 'uwot' package, which isn't installed."))
      n_comp <- min(10, nrow(mat) - 1, ncol(mat))
      validate(need(n_comp >= 2, "Not enough samples/bins in the current selection to compute at least two UMAP dimensions."))
      validate(need(xdim <= n_comp && ydim <= n_comp,
                    paste0("Only ", n_comp, " UMAP dimension(s) can be computed for the current selection. Pick a dimension up to UMAP", n_comp, ".")))
      set.seed(42)
      n_neighbors <- max(2, min(15, nrow(mat) - 1))
      um <- uwot::umap(mat, n_neighbors = n_neighbors, n_components = n_comp)
      coords <- as.data.frame(um[, c(xdim, ydim), drop = FALSE])
      names(coords) <- c("Dim1", "Dim2")
      xlab <- paste0("UMAP ", xdim)
      ylab <- paste0("UMAP ", ydim)
    }
    
    coords$sample_id <- rownames(mat)
    coords <- merge(coords, metadata[, c("sample_id", "patient_id", "Type", "sexe")], by = "sample_id")
    
    p <- plot_ly(coords, x = ~Dim1, y = ~Dim2, color = ~Type,
                 colors = c("Tumor" = "#d1495b", "Normal" = "#3aa9c9"),
                 type = "scatter", mode = "markers",
                 marker = list(size = 10, line = list(color = "#0b2436", width = 1)),
                 text = ~paste0("Sample: ", sample_id, "<br>Patient: ", patient_id, "<br>Type: ", Type),
                 hoverinfo = "text") |>
      layout(xaxis = list(title = xlab), yaxis = list(title = ylab), legend = list(orientation = "h", x = 0, y = 1.08))
    list(plot = p, data = coords)
  })
  output$tn_projection <- renderPlotly({ tn_projection_build()$plot })
  register_plotly_export(output, "tn_projection", function() tn_projection_build()$data)
  
  tn_network_build <- reactive({
    mat <- network_sample_matrix()
    validate(need(!is.null(mat), "Not enough overlapping data across samples for the selected bins to build a similarity network."))
    
    # one or more clinical/annotation variables to colour by (small multiples if more than one is picked)
    features <- input$network_color_feature
    features <- intersect(features, names(sample_annotation))
    if (length(features) == 0) features <- "Type"
    
    # user-configurable similarity cutoff for drawing an edge (e.g. corr >= 0.7)
    corr_threshold <- input$network_corr_threshold
    if (is.null(corr_threshold) || !is.finite(corr_threshold)) corr_threshold <- 0.7
    
    cor_mat <- stats::cor(t(mat)) # Pearson correlation for similarity between samples
    layout_xy <- tryCatch(stats::cmdscale(stats::as.dist(1 - cor_mat), k = 2), error = function(e) NULL)
    validate(need(!is.null(layout_xy), "Couldn't compute a layout for the current selection."))
    
    nodes_base <- data.frame(sample_id = rownames(mat), x = layout_xy[, 1], y = layout_xy[, 2], stringsAsFactors = FALSE)
    nodes_base <- merge(nodes_base, sample_annotation, by = "sample_id", all.x = TRUE)
    nodes_base$Type[is.na(nodes_base$Type)] <- "unknown"
    
    # only draw an edge between sample pairs whose correlation meets the chosen threshold
    ids <- rownames(mat)
    edge_df <- NULL
    if (length(ids) >= 2) {
      pairs <- utils::combn(ids, 2)
      edge_df <- data.frame(
        from = pairs[1, ], to = pairs[2, ],
        corr = cor_mat[cbind(pairs[1, ], pairs[2, ])],
        stringsAsFactors = FALSE)
      edge_df <- edge_df[is.finite(edge_df$corr) & edge_df$corr >= corr_threshold, ]
      if (nrow(edge_df) > 0) {
        edge_df <- merge(edge_df, setNames(nodes_base[, c("sample_id", "x", "y")], c("from", "x_from", "y_from")), by = "from")
        edge_df <- merge(edge_df, setNames(nodes_base[, c("sample_id", "x", "y")], c("to", "x_to", "y_to")), by = "to")
      }
    }
    validate(need(nrow(nodes_base) > 0, "No samples available for the current bin selection."))
    
    type_shapes <- c("Tumor" = 21, "Normal" = 24, "unknown" = 22)
    
    # builds one coloured network panel for a single variable; reused per feature so several can be shown side by side
    make_network_panel <- function(feature, show_subtitle) {
      nodes <- nodes_base
      nodes$color_label <- as.character(nodes[[feature]])
      nodes$color_label[is.na(nodes$color_label)] <- "unknown"
      color_palette <- categorical_palette(nodes$color_label)
      feature_label <- names(network_color_choices)[network_color_choices == feature][1]
      
      p <- ggplot()
      if (!is.null(edge_df) && nrow(edge_df) > 0) {
        p <- p + geom_segment(
          data = edge_df, aes(x = x_from, y = y_from, xend = x_to, yend = y_to, alpha = corr),
          color = "#7d92a3", linewidth = 0.4, show.legend = FALSE)
      }
      p + geom_point(data = nodes, aes(x = x, y = y, fill = color_label, shape = Type), size = 5, color = "#0b2436", stroke = 0.6) +
        geom_text(data = nodes, aes(x = x, y = y, label = patient_id), vjust = -1.3, size = 2.8, color = "#16324f") +
        scale_fill_manual(values = color_palette, name = feature_label) +
        scale_shape_manual(values = type_shapes, name = "Type", drop = TRUE) +
        guides(fill = guide_legend(override.aes = list(shape = 21))) +
        labs(x = "Dimension 1", y = "Dimension 2",
             title = if (length(features) > 1) feature_label else "Patient similarity network (methylation correlation)",
             subtitle = if (show_subtitle) paste0(region_label(), " \u2014 edges at corr \u2265 ", corr_threshold) else NULL) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(t = 0, b = 0),
              panel.grid.minor = element_blank())
    }
    
    p <- if (length(features) == 1) {
      make_network_panel(features[1], show_subtitle = TRUE)
    } else {
      panels <- lapply(features, make_network_panel, show_subtitle = FALSE)
      n_col <- if (length(panels) >= 3) 2 else length(panels)
      patchwork::wrap_plots(panels, ncol = n_col, guides = "collect") +
        patchwork::plot_annotation(
          title = "Patient similarity network (methylation correlation)",
          subtitle = paste0(region_label(), " \u2014 edges at corr \u2265 ", corr_threshold, "; circle = Tumor, triangle = Normal"),
          theme = theme(plot.title = element_text(face = "bold", size = 13, color = "#16324f"),
                        plot.subtitle = element_text(size = 10, color = "#7d92a3"))) &
        theme(legend.position = "bottom")
    }
    node_export <- nodes_base[, intersect(c("sample_id", "patient_id", "Type", "x", "y", features), names(nodes_base))]
    list(plot = p, data = node_export)
  })
  output$tn_network <- renderPlot({ tn_network_build()$plot })
  register_plot_export(output, "tn_network",
                       function() tn_network_build()$plot,
                       function() tn_network_build()$data,
                       width = 10, height = 7)
  
  # 4. Genome-wide Profile
  output$manhattan_plot <- renderPlotly({
    df <- genome_wide_bins
    validate(need(nrow(df) > 0, "No genome-wide methylation data available."))
    
    cat_colors <- c("Not significant" = "#168AAD","Extreme outlier" = "#F4A261","Significant hypomethylation" = "#76C893","Significant hypermethylation" = "#8D99AE")
    
    plot_ly(df, x = ~genome_x, y = ~mean_diff_tumor_normal, color = ~manhattan_category,
            colors = cat_colors, type = "scattergl", mode = "markers",
            marker = list(size = 6, opacity = 0.9, line = list(width = 0)),
            text = ~paste0(
              "Bin: ", bin_id,
              "<br>Chr: ", as.character(chr),
              "<br>Diff (Tumor \u2212 Normal): ", round(mean_diff_tumor_normal, 3),
              "<br>q-value: ", ifelse(is.na(q_value), "NA", format(signif(q_value, 3), scientific = TRUE)),
              "<br>Patients tested: ", ifelse(is.na(n), 0, n)),
            hoverinfo = "text") |>
      layout(xaxis = list(title = "", tickvals = as.numeric(chr_mid), ticktext = chrom_list, tickangle = 45),
             yaxis = list(title = "Methylation difference (Tumor \u2212 Normal)", zeroline = TRUE, zerolinewidth = 1.4, zerolinecolor = "#16324f"),
             shapes = chr_boundary_shapes,
             legend = list(orientation = "h", x = 0, y = 1.12),
             hovermode = "closest") |>
      config(displayModeBar = TRUE,
             toImageButtonOptions = list(format = "png", filename = "genome_wide_methylation_profile"))
  })
  register_plotly_export(output, "manhattan_plot",
                         function() genome_wide_bins[, c("bin_id", "chr", "bin_position", "mean_diff_tumor_normal", "q_value", "n", "manhattan_category")])
  
  # 5. Feature × Chromosome Heatmap
  max_dendro_bins <- 300
  
  # current zoom level (%), falling back to 100 before the slider has initialized
  heatmap_zoom_pct <- reactive({
    z <- input$heatmap_zoom
    if (is.null(z) || !is.finite(z)) 100 else z
  })
  
  observeEvent(input$heatmap_zoom_reset, {
    updateSliderInput(session, "heatmap_zoom", value = 100)
  })
  
  observeEvent(input$heatmap_container_dims, {
    d <- heatmap_build()
    natural_w <- max(400, round(d$n_bins * 14 + 160))
    natural_h <- max(300, round(d$n_patients * 20 + 160))
    avail_w <- input$heatmap_container_dims$w
    avail_h <- input$heatmap_container_dims$h
    if (!is.null(avail_w) && !is.null(avail_h) && is.finite(avail_w) && is.finite(avail_h) &&
        natural_w > 0 && natural_h > 0) {
      fit_pct <- floor(min(avail_w / natural_w, avail_h / natural_h) * 100)
      fit_pct <- max(10, min(400, fit_pct))
      updateSliderInput(session, "heatmap_zoom", value = fit_pct)
    }
  })
  
  # all the data prep (matrix, clustering, dendrogram segments) as one reactive
  heatmap_build <- reactive({
    feature <- input$heatmap_feature
    if (is.null(feature)) feature <- "none"
    
    bins_now <- selected_bins()
    validate(need(length(bins_now) > 0, "Select at least one bin in the sidebar."))
    
    shift_df <- patient_bin_shift[patient_bin_shift$bin_id %in% bins_now, ]
    validate(need(nrow(shift_df) > 0, "No methylation shift data for the selected bins."))
    max_abs <- max(abs(shift_df$shift), na.rm = TRUE)
    
    # patient (row) x bin (column) matrix of shifts, used to cluster BOTH axes by similarity
    bins_present <- bin_table$bin_id[bin_table$bin_id %in% unique(shift_df$bin_id)]
    mat <- build_shift_matrix(shift_df, bins_present)
    mat[is.na(mat)] <- 0
    validate(need(nrow(mat) >= 2, "Need at least two patients with data to build a clustered heatmap."))
    
    # row dendrogram: cluster patients by their shift pattern across the selected bins
    row_hc <- hclust(dist(mat))
    row_order <- rownames(mat)[row_hc$order]
    row_segments <- hclust_to_segments(row_hc)
    n_patients <- length(row_order)
    
    # column dendrogram: cluster bins by their shift pattern across patients
    show_col_dendro <- ncol(mat) >= 2 && ncol(mat) <= max_dendro_bins
    if (show_col_dendro) {
      col_hc <- hclust(dist(t(mat)))
      col_order <- colnames(mat)[col_hc$order]
      col_segments <- hclust_to_segments(col_hc)
    } else {
      col_order <- bin_table$bin_id[bin_table$bin_id %in% colnames(mat)]
      col_segments <- NULL
    }
    n_bins <- length(col_order)
    
    shift_df$bin_id <- factor(shift_df$bin_id, levels = col_order)
    shift_df$patient_id <- factor(shift_df$patient_id, levels = row_order)
    shift_df$x <- as.integer(shift_df$bin_id)
    shift_df$y <- as.integer(shift_df$patient_id)
    
    list(feature = feature, shift_df = shift_df, max_abs = max_abs, row_order = row_order,
         col_order = col_order, n_patients = n_patients, n_bins = n_bins,
         row_segments = row_segments, col_segments = col_segments, show_col_dendro = show_col_dendro)
  })
  
  feature_heatmap_build <- reactive({
    d <- heatmap_build()
    
    # main heatmap uses numeric x/y (rather than discrete scales) so its tile positions line up
    p_main <- ggplot(d$shift_df, aes(x = x, y = y, fill = shift)) +
      geom_tile(color = "white", linewidth = 0.15) +
      scale_fill_gradient2(low = "#3aa9c9", mid = "white", high = "#d1495b", midpoint = 0, limits = c(-d$max_abs, d$max_abs), name = "Methylation\nshift\n(Tumor \u2212 Normal)") +
      scale_x_continuous(breaks = seq_len(d$n_bins), labels = d$col_order, expand = c(0, 0)) +
      scale_y_continuous(breaks = seq_len(d$n_patients), labels = d$row_order, expand = c(0, 0)) +
      labs(x = if (d$show_col_dendro) "Bin (dendrogram-ordered)" else "Bin (genome order; too many bins to cluster)",
           y = "Patient (dendrogram-ordered)") +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))
    
    row_dendro <- plot_dendrogram(d$row_segments, orientation = "left") +
      scale_y_continuous(limits = c(0.5, d$n_patients + 0.5), expand = c(0, 0))
    
    col_dendro <- if (d$show_col_dendro) {
      plot_dendrogram(d$col_segments, orientation = "top") +
        scale_x_continuous(limits = c(0.5, d$n_bins + 0.5), expand = c(0, 0))
    } else {
      plot_spacer()
    }
    
    heatmap_export_data <- d$shift_df[, c("patient_id", "bin_id", "chr", "meth_tumor", "meth_normal", "shift")]
    
    if (identical(d$feature, "none")) {
      design <- "
      AB
      CD
      "
      p_none <- plot_spacer() + col_dendro + row_dendro + p_main +
        patchwork::plot_layout(design = design, widths = c(1, 8), heights = c(1, 8))
      return(list(plot = p_none, data = heatmap_export_data))
    }
    
    # colour strip for the chosen clinical feature, row-aligned with the main heatmap's row order
    ann_df <- patient_annotation[patient_annotation$patient_id %in% d$row_order, ]
    ann_df$patient_id <- factor(ann_df$patient_id, levels = d$row_order)
    ann_df$y <- as.integer(ann_df$patient_id)
    ann_df$x <- "Feature"
    ann_df$value <- as.character(ann_df[[d$feature]])
    ann_df$value[is.na(ann_df$value)] <- "unknown"
    
    levels_cat <- sort(unique(ann_df$value))
    base_palette <- c("#0e7c86", "#d1495b", "#3aa9c9", "#e0a339", "#2fae66", "#16324f", "#7d92a3")
    pal <- setNames(base_palette[seq_along(levels_cat)], levels_cat)
    
    feature_label <- names(heatmap_feature_choices)[heatmap_feature_choices == d$feature]
    
    p_ann <- ggplot(ann_df, aes(x = x, y = y, fill = value)) +
      geom_tile(color = "white", linewidth = 0.15) +
      scale_fill_manual(values = pal, name = feature_label) +
      scale_y_continuous(limits = c(0.5, d$n_patients + 0.5), expand = c(0, 0)) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 12) +
      theme(axis.text.y = element_blank(), axis.text.x = element_blank(), panel.grid = element_blank())
    
    design <- "
    AAB
    CDE
    "
    p_feat <- plot_spacer() + col_dendro + row_dendro + p_ann + p_main +
      patchwork::plot_layout(design = design, widths = c(1, 0.6, 8), heights = c(1, 8)) # plot_layout() from patchwork package
    list(plot = p_feat, data = heatmap_export_data)
  })
  
  output$feature_heatmap <- renderPlot({
    feature_heatmap_build()$plot
  },
  # the plot's rendered pixel size is the "natural" data-driven size
  height = function() {
    d <- heatmap_build()
    natural_h <- max(300, round(d$n_patients * 20 + 160))
    max(150, round(natural_h * heatmap_zoom_pct() / 100))
  },
  width = function() {
    d <- heatmap_build()
    natural_w <- max(400, round(d$n_bins * 14 + 160))
    max(200, round(natural_w * heatmap_zoom_pct() / 100))
  })
  register_plot_export(output, "feature_heatmap",
                       function() feature_heatmap_build()$plot,
                       function() feature_heatmap_build()$data,
                       width = 12, height = 9)
  
  # 6. Clinical Explorer
  
  # filtered clinical profile, driven by the "Filter Patients" panel; used by both the table and the boxplot
  filtered_clinical_profile <- reactive({
    apply_clinical_filters(clinical_profile_table, clinical_numeric_cols, clinical_categorical_cols, input)
  })
  
  observeEvent(input$clinfilter_reset, {
    for (col in clinical_categorical_cols) {
      updateCheckboxGroupInput(session, paste0("clinfilter_", col), selected = clinical_categorical_choices(clinical_profile_table[[col]]))
    }
    for (col in clinical_numeric_cols) {
      vals <- clinical_profile_table[[col]]
      vals <- vals[is.finite(vals)]
      rng <- if (length(vals) > 0) range(vals) else c(0, 1)
      updateSliderInput(session, paste0("clinfilter_", col), value = rng)
    }
  })
  
  # live "N of TOTAL patients" badge shown in the Clinical metadata explorer card header
  output$clinical_table_count_badge <- renderUI({
    n <- nrow(filtered_clinical_profile())
    total <- nrow(clinical_profile_table)
    span(style = "font-size:12px; font-weight:600; color:#16324f; background:#eef2f5; border:1px solid #d8e0e6; border-radius:20px; padding:4px 12px; white-space:nowrap;",
         paste0(format(n, big.mark = ","), " of ", format(total, big.mark = ","), " patients"))
  })
  
  # table data as its own reactive so DT row indices (input$clinical_table_rows_selected)
  clinical_table_data <- reactive({
    filtered_clinical_profile()
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
  
  clinical_boxplot_build <- reactive({
    gene <- input$mutation_gene
    req(gene)
    test_method <- input$mutation_stat_test
    if (is.null(test_method)) test_method <- "wilcox"
    
    status_labels <- c("0" = "Wild-type", "1" = "Mutant")
    filtered_ids <- filtered_clinical_profile()$patient_id
    
    if (view_by_samples()) {
      # one row per Tumor sample x selected bin (unaggregated), instead of one row per patient
      validate(need(gene %in% names(patient_annotation), paste0(gene, " status is not available in this dataset.")))
      bins_now <- selected_bins()
      df <- ml_annot[ml_annot$bin_id %in% bins_now & ml_annot$Type == "Tumor", c("sample_id", "patient_id", "bin_id", "methylation")]
      df <- df[df$patient_id %in% filtered_ids, ]
      df <- merge(df, patient_annotation[, c("patient_id", gene)], by = "patient_id")
      names(df)[names(df) == gene] <- "GeneStatus"
      status_raw <- as.character(df$GeneStatus)
      df$Status <- ifelse(status_raw %in% names(status_labels), status_labels[status_raw], NA_character_)
      df$Value <- df$methylation
      df <- df[!is.na(df$Status) & is.finite(df$Value), ]
      validate(need(nrow(df) > 0, "No patients match the current filters. Adjust the filters above."))
      validate(need(length(unique(df$Status)) == 2,
                    paste0("Need both Wild-type and Mutant tumor samples for ", gene, " in the current bin selection.")))
      y_lab <- "Tumor methylation (per sample \u00d7 bin)"
    } else {
      df <- overview_patient_shift()
      validate(need(!is.null(df) && nrow(df) > 0, "No data available for the current bin selection."))
      validate(need(gene %in% names(df), paste0(gene, " status is not available in this dataset.")))
      
      df <- df[df$patient_id %in% filtered_ids, ]
      validate(need(nrow(df) > 0, "No patients match the current filters. Adjust the filters above."))
      
      status_raw <- as.character(df[[gene]])
      df$Status <- ifelse(status_raw %in% names(status_labels), status_labels[status_raw], NA_character_)
      df$Value <- df$Tumor
      df <- df[!is.na(df$Status) & is.finite(df$Value), ]
      validate(need(length(unique(df$Status)) == 2,
                    paste0("Need both Wild-type and Mutant tumor samples for ", gene, " in the current bin selection.")))
      y_lab <- "Mean tumor methylation (per bin)"
    }
    df$Status <- factor(df$Status, levels = c("Wild-type", "Mutant"))
    
    n_by_status <- table(df$Status)
    test_res <- run_group_test(df$Value[df$Status == "Mutant"], df$Value[df$Status == "Wild-type"], method = test_method)
    y_max <- max(df$Value, na.rm = TRUE)
    
    p <- ggplot(df, aes(x = Status, y = Value, fill = Status)) +
      geom_jitter(width = 0.08, size = 1.1, alpha = 0.35, color = "#16324f") +
      geom_boxplot(width = 0.45, outlier.shape = NA, alpha = 0.92, color = "#16324f", linewidth = 0.4) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = "white", color = "#16324f", stroke = 0.8) +
      scale_fill_manual(values = c("Wild-type" = "#3aa9c9", "Mutant" = "#d1495b")) +
      scale_x_discrete(labels = paste0(levels(df$Status), "\n(n=", as.numeric(n_by_status[levels(df$Status)]), ")")) +
      annotate("text", x = 1.5, y = y_max * 1.1, label = test_res$label, size = 3.5, color = "#16324f", fontface = "bold") +
      scale_y_continuous(limits = c(0, max(1, y_max * 1.18)), expand = expansion(mult = c(0.02, 0.02))) +
      labs( x = NULL, y = y_lab, title = paste0(gene, " mutation status vs. tumor methylation"), subtitle = region_label()) +
      theme_app() +
      theme(legend.position = "none")
    
    sel_id <- selected_patient_id()
    if (!is.null(sel_id) && sel_id %in% df$patient_id) {
      hl <- df[df$patient_id == sel_id, ]
      p <- p + geom_point(data = hl, aes(x = Status, y = Value), shape = 21, size = 4.8, fill = "#e0a339", color = "#16324f", stroke = 1.1) +
        geom_text(data = hl, aes(x = Status, y = Value, label = paste0("Patient ", patient_id)), vjust = -1.3, size = 3, color = "#16324f", fontface = "bold")
    }
    export_cols <- intersect(c("patient_id", "sample_id", "bin_id", "Status", "Value", "Normal", "shift"), names(df))
    list(plot = p, data = df[, export_cols])
  })
  output$clinical_boxplot <- renderPlot({ clinical_boxplot_build()$plot })
  register_plot_export(output, "clinical_boxplot",
                       function() clinical_boxplot_build()$plot,
                       function() clinical_boxplot_build()$data)
  
  output$clinical_table <- renderDT({
    df <- clinical_table_data()
    display_df <- df
    names(display_df) <- vapply(names(df), clinical_label_for, character(1))
    num_cols <- vapply(display_df, is.numeric, logical(1))
    display_df[num_cols] <- lapply(display_df[num_cols], round, 3)
    datatable(display_df,filter = "top",  selection = "single",rownames = FALSE,  options = list(scrollX = TRUE, pageLength = 10, autoWidth = TRUE),class = "stripe hover compact")
  })
  
  output$selected_patient_banner <- renderUI({
    sel_id <- selected_patient_id()
    if (is.null(sel_id)) {
      return(div(style = "font-size:12px; color:#7d92a3; background:#eef2f5; border:1px dashed #c3ccd3; border-radius:8px; padding:10px 16px; margin-bottom:14px;",
                 bs_icon("info-circle"), " No patient selected, click a row in the clinical metadata table to highlight that patient in the plot."))
    }
    row <- clinical_profile_table[clinical_profile_table$patient_id == sel_id, ]
    if (nrow(row) == 0) return(NULL)
    
    chip <- function(col_id, fallback = "\u2014") {
      raw <- if (col_id %in% names(row) && length(row[[col_id]]) == 1) row[[col_id]] else NA
      val <- if (is.na(raw)) fallback else if (is.numeric(raw)) format(round(raw, 3), nsmall = 3) else as.character(raw)
      div(style = "display:flex; flex-direction:column; gap:2px;",
          tags$span(clinical_label_for(col_id), style = "font-size:10px; color:#7d92a3; text-transform:uppercase; letter-spacing:0.4px;"),
          tags$span(val, style = "font-size:14px; font-weight:600; color:#16324f;"))
    }
    
    chip_cols <- intersect(c("sexe", "estadi2", "MSS", "BRAF", "KRAS", "TP53", "Recaiguda", "MethShift"), names(row))
    div(style = "display:flex; align-items:center; gap:26px; flex-wrap:wrap; background:#eef2f5; border:1px solid #d8e0e6; border-radius:8px; padding:12px 18px; margin-bottom:14px;",
        div(style = "display:flex; align-items:center; gap:8px; font-weight:700; color:#16324f;",
            bs_icon("person-check-fill", size = "1.3em"), paste0("Patient ", sel_id)),
        lapply(chip_cols, chip),
        actionButton("clear_patient_selection", "Clear selection", icon = bs_icon("x-circle"),
                     class = "btn-sm btn-outline-secondary", style = "margin-left:auto;"))
  })
  
  # 7. Bin Table
  
  # note shown if gene annotation couldn't be computed at all
  output$bintable_gene_note <- renderUI({
    if (!gene_annotation_available) {
      return(div(style = "font-size:12px; color:#a34; background:#fdf1ef; border:1px solid #f0d6d1; border-radius:8px; padding:8px 14px; max-width:520px;",
                 bs_icon("exclamation-triangle")," Gene annotation source files were not found when the app started"))
    }
  })
  
  # note shown if promoter annotation couldn't be computed at all
  output$bintable_promoter_note <- renderUI({
    if (!promoter_annotation_available) {
      return(div(style = "font-size:12px; color:#a34; background:#fdf1ef; border:1px solid #f0d6d1; border-radius:8px; padding:8px 14px; max-width:520px;",
                 bs_icon("exclamation-triangle")," Promoter annotation source file was not found when the app started"))
    }
  })
  
  output$bintable <- renderDT({
    display_df <- build_display_bin_table(filtered_bin_table())
    meth_cols  <- intersect(c("Tumor Mean", "Normal Mean", "Tumor SD", "Normal SD", "\u0394 (Tumor \u2212 Normal)"), names(display_df))
    pct_cols   <- intersect(c("Alu Meth %", "CpG Meth %"), names(display_df))
    count_cols <- intersect(c("Alu Count", "CpG Count", "Gene Count"), names(display_df))
    delta_col  <- "\u0394 (Tumor \u2212 Normal)"
    links_col  <- which(names(display_df) == "Quick Links")
    
    dt <- datatable(display_df,
                    container = bin_table_header_sketch,
                    rownames = FALSE,
                    class = "stripe hover compact",
                    escape = if (length(links_col) > 0) -links_col else TRUE,
                    options = list(scrollX = TRUE, pageLength = 15, lengthMenu = list(c(10, 15, 25, 50, -1), c("10", "15", "25", "50", "All")), dom = "ltip", language = list(search = "Quick search:", lengthMenu = "Show _MENU_ bins per page"), columnDefs = list(list(className = "dt-center", targets = "_all"))))
    
    if (length(meth_cols) > 0) dt <- formatRound(dt, meth_cols, 3)
    if (length(pct_cols) > 0) dt <- formatRound(dt, pct_cols, 2)
    if (length(count_cols) > 0) dt <- formatRound(dt, count_cols, 0)
    
    # colour the tumor/normal shift so a hyper- vs hypo-methylated bin is visible at a glance
    if (delta_col %in% names(display_df)) {
      dt <- formatStyle(dt, delta_col, fontWeight = "600", color = styleInterval(0, c("#0e7c86", "#d1495b")))
    }
    # make bins that carry one or more COSMIC cancer census genes stand out from the majority with none
    if ("Gene Count" %in% names(display_df)) {
      dt <- formatStyle(dt, "Gene Count", fontWeight = styleInterval(0, c("normal", "700")), color = styleInterval(0, c("#adb8c0", "#16324f")), backgroundColor = styleInterval(0, c("transparent", "#e8f5f3")))
    }
    dt
  })
  
  output$bintable_download <- downloadHandler(filename = function() {"bintable.csv"},content = function(file) {
    write.csv(build_display_bin_table(filtered_bin_table(), plain = TRUE), file, row.names = FALSE)
  })
}

shinyApp(ui = ui, server = server)