# Exploration of Epigenomic Data in Colorectal Cancer
# Paula Artiz Dueñas - UPC

library(shiny)
library(bslib)
library(bsicons) # icons used for value boxes
library(DT)
library(plotly)
library(ggplot2)

data <- readRDS("/Users/paulaartizduenas/Desktop/Project/Dataset/Data Processed/data_app.rds")

metadata <- data$metadata
bin_table <- data$bin_table
methylation_long <- data$methylation_long
chrom_list <- c(as.character(1:22), "X", "Y")

# Ordered bin choices for the bin-selector (grouped by chromosome, sorted by position)
bin_table <- bin_table[order(match(bin_table$chr, chrom_list), bin_table$bin_position), ]
bin_choices_by_chr <- split(bin_table$bin_id, factor(bin_table$chr, levels = chrom_list))

# Bin Table tab: derived columns + registries

# Human-readable "start–end" label for a bin, e.g. "1 MB-2 MB"
format_mb_label <- function(bp) {
  mb <- bp / 1e6
  if (abs(mb - round(mb)) < 1e-9) paste0(round(mb)," MB")
  else paste0(format(round(mb, 1), nsmall = 1)," MB")
}
format_bin_coordinates <- function(start, end) {
  paste0(format_mb_label(start), "\u2013", format_mb_label(end))
}
bin_table$bin_coordinates <- mapply(format_bin_coordinates, bin_table$bin_start, bin_table$bin_end)

# Per-sample tumor/normal methylation columns used to compute per-bin SD
tumor_sample_pattern  <- "^tumor_sample"
normal_sample_pattern <- "^normal_sample"
tumor_sample_cols  <- grep(tumor_sample_pattern, names(bin_table), value = TRUE)
normal_sample_cols <- grep(normal_sample_pattern, names(bin_table), value = TRUE)

# Standard deviation not found in data_app.rds
if (!"sd_methylation_tumor" %in% names(bin_table)) bin_table$sd_methylation_tumor <- NA_real_
if (!"sd_methylation_normal" %in% names(bin_table)) bin_table$sd_methylation_normal <- NA_real_

# Mean Difference Metric (Tumor - Normal)
bin_table$mean_diff_tumor_normal <- bin_table$mean_methylation_tumor - bin_table$mean_methylation_normal

# Counts of Alus and CpGs per 1 Mb bin, computed upstream in Week_2_With_Prevalence.R
if (!"n_alu" %in% names(bin_table)) bin_table$n_alu <- NA_real_
if (!"n_cpg" %in% names(bin_table)) bin_table$n_cpg <- NA_real_

# Columns shown in the Bin Table tab, in display order
bin_table_columns <- list(
  list(id = "bin_id", label = "Bin ID"),
  list(id = "bin_coordinates", label = "Bin Coordinates"),
  list(id = "mean_methylation_tumor", label = "Mean Meth (Tumor)", digits = 3),
  list(id = "mean_methylation_normal", label = "Mean Meth (Normal)", digits = 3),
  list(id = "sd_methylation_tumor", label = "SD Meth (Tumor)", digits = 3),
  list(id = "sd_methylation_normal",label = "SD Meth (Normal)", digits = 3),
  list(id = "mean_diff_tumor_normal",  label = "Mean Diff (Tumor - Normal)", digits = 3),
  list(id = "n_alu", label = "Alu Count"),
  list(id = "n_cpg", label = "CpG Count"))

# Numeric metrics exposed as filters
bin_table_filters <- list(
  list(id = "mean_methylation_tumor", label = "Mean Meth (Tumor)", op = ">="),
  list(id = "mean_methylation_normal", label = "Mean Meth (Normal)", op = ">="),
  list(id = "sd_methylation_tumor", label = "SD Meth (Tumor)", op = ">="),
  list(id = "sd_methylation_normal",label = "SD Meth (Normal)", op = ">="),
  list(id = "mean_diff_tumor_normal", label = "Mean Diff (Tumor - Normal)", op = ">=")
)

# Builds visual grid of slider inputs for the Bin Table filters
build_bin_table_filter_inputs <- function(filters, df) {
  lapply(filters, function(f) {
    vals <- df[[f$id]]
    vals <- vals[is.finite(vals)]
    rng <- if (length(vals) > 0) range(vals) else c(0, 1)
    
    lo <- floor(rng[1] * 100) / 100
    hi <- ceiling(rng[2] * 100) / 100
    if (hi <= lo) hi <- lo + 0.01
    
    div(
      style = "flex: 1 1 200px; min-width: 180px; max-width: 250px;",
      sliderInput(
        inputId = paste0("binfilter_", f$id),
        label = paste0(f$label, " ", f$op),
        min = lo, max = hi, value = lo, step = 0.01,
        width = "100%"
      )
    )
  })
}

# Applies all configured filters to a bin_table subset.
apply_bin_table_filters <- function(df, filters, input) {
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
  df
}

# Builds the display data frame from bin_table_columns.
build_display_bin_table <- function(df, columns = bin_table_columns) {
  col_ids <- vapply(columns, function(c) c$id, character(1))
  out <- df[, col_ids, drop = FALSE]
  for (i in seq_along(columns)) {
    if (!is.null(columns[[i]]$digits)) out[[i]] <- round(out[[i]], columns[[i]]$digits)
  }
  names(out) <- vapply(columns, function(c) c$label, character(1))
  out
}

# Length of each chromosome
chrom_lengths <- tapply(bin_table$bin_end, bin_table$chr, max)
chrom_lengths <- setNames(as.numeric(chrom_lengths[chrom_list]), chrom_list)

# Cumulative genome-wide offset
chr_offset <- setNames(numeric(length(chrom_list)), chrom_list)
running <- 0
for (chr in chrom_list) {
  chr_offset[chr] <- running
  running <- running + chrom_lengths[[chr]]
}
total_genome_len <- running
chr_mid <- chr_offset + chrom_lengths / 2

# Genome-wide bin table
genome_wide_bins <- bin_table
genome_wide_bins$chr <- factor(genome_wide_bins$chr, levels = chrom_list)
genome_wide_bins <- genome_wide_bins[order(genome_wide_bins$chr, genome_wide_bins$bin_position), ]
genome_wide_bins$genome_x <- chr_offset[as.character(genome_wide_bins$chr)] + genome_wide_bins$bin_position
genome_wide_bins$overall_meth <- rowMeans(
  genome_wide_bins[, c("mean_methylation_tumor", "mean_methylation_normal")], na.rm = TRUE)
genome_wide_bins$chr_parity <- factor(as.integer(genome_wide_bins$chr) %% 2)

# Dashed vertical lines marking chromosome boundaries, reused on both plots
chr_boundary_shapes <- lapply(as.numeric(chr_offset)[-1], function(x) {
  list(type = "line", x0 = x, x1 = x, y0 = 0, y1 = 1, yref = "paper",
       line = list(color = "rgba(150,150,150,0.4)", width = 1, dash = "dot"))
})

# Static lookup table of hg19/GRCh37 coordinates for colorectal-cancer driver genes.
# Used only as a fallback when a searched gene isn't already annotated in bin_table$genes.
# Gene set expanded (82 genes) with https://www.intogen.org/search and UCSC Genome Browser 

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

# Parses a free-text genome browser search box into a (chr, start, end) target
parse_genomic_search <- function(query, bin_table, chrom_list, gene_lookup_table, pad = 2e6) {
  q <- trimws(query)
  if (!nzchar(q)) return(list(status = "error", message = "Please enter a coordinate, bin ID, or gene name."))
  q_nospace <- gsub(",", "", q)
  
  # Case 1: chr:start[-end] coordinate format
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
  
  # Case 2: bin ID format "chr_position"
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
  
  # Case 3: gene symbol already annotated on one or more bins
  gene_query <- toupper(q)
  has_annotation <- !is.na(bin_table$genes) &
    grepl(paste0("(^|[,;])\\s*", gene_query, "\\s*($|[,;])"), toupper(bin_table$genes))
  in_annotation <- bin_table[has_annotation, ]
  if (nrow(in_annotation) >= 1) {
    return(list(status = "ok", chr = in_annotation$chr[1],
                start = min(in_annotation$bin_start), end = max(in_annotation$bin_end),
                message = paste0("Jumped to gene ", gene_query)))
  }
  
  # Case 4: gene symbol not in bin_table, fall back to the static gene_lookup_table
  in_lookup <- gene_lookup_table[toupper(gene_lookup_table$gene) == gene_query, ]
  if (nrow(in_lookup) == 1) {
    return(list(
      status = "ok", chr = in_lookup$chr,
      start = max(1, in_lookup$start - pad), end = in_lookup$end + pad,
      message = paste0("Jumped to ", gene_query, " (", in_lookup$chr, ":",
                       format(in_lookup$start, big.mark = ",", scientific = FALSE), "-",
                       format(in_lookup$end, big.mark = ",", scientific = FALSE), ")")
    ))
  }
  list(status = "error", message = paste0("'", q, "' was not recognized as a coordinate, bin ID, or gene name."))
}

# Parses an uploaded file of bins (.csv, .tsv, .txt, or unrecognized extension) into a
# vector of bin_id values that exist in bin_table

parse_bin_file <- function(filepath, filename, bin_table, chrom_list) {
  ext <- tolower(tools::file_ext(filename))
  raw_ids <- character(0)
  tbl <- NULL
  
  # Only attempt structured table
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
  
  # Fall back to treating the file as plain, unstructured text
  if (length(raw_ids) == 0) {
    txt <- tryCatch(readLines(filepath, warn = FALSE), error = function(e) character(0))
    txt <- paste(txt, collapse = "\n")
    raw_ids <- strsplit(txt, "[,;\\s]+", perl = TRUE)[[1]]
  }
  
  raw_ids <- trimws(raw_ids)
  raw_ids <- raw_ids[nzchar(raw_ids)]
  n_entries <- length(unique(raw_ids))
  
  # Normalize formatting differences: drop "chr" prefix, unify ":"/"-" separators to "_"
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

# Zooms in and out
zoom_view <- function(start, end, chr_len, factor, min_width = 2e6) {
  width <- end - start
  center <- (start + end) / 2
  new_width <- max(min_width, min(chr_len, width * factor)) # factor < 1 zooms in, > 1 zooms out
  new_start <- center - new_width / 2
  new_end <- center + new_width / 2
  if (new_start < 1) {new_end <- new_end + (1 - new_start); new_start <- 1}
  if (new_end > chr_len) {new_start <- new_start - (new_end - chr_len); new_end <- chr_len}
  new_start <- max(1, new_start)
  list(start = round(new_start), end = round(new_end))}

# UI Definition

ui <- page_sidebar(
  title = div(
    style = "display: flex; justify-content: space-between; align-items: center; width: 100%; gap: 1px; padding: 22px 5px; border-bottom: 1px solid rgba(255,255,255,0.08);",
    div(
      style = "flex:1 1 auto; min-width:0;",
      h2("Exploration of Epigenomic Data in Colorectal Cancer", style = "margin:0; font-weight:600; letter-spacing:0.2px; white-space:normal;")
    ),
    div(
      style = "display:flex; gap:14px; align-items:center; margin-left:auto;",
      div(
        style = "display:flex; align-items:center; gap:8px; border:1px solid #3a5470; border-radius:8px; padding:8px 14px;",
        bs_icon("person", size = "2em"),
        textOutput("n_patients", inline = TRUE), " patients"
      ),
      div(
        style = "display:flex; align-items:center; gap:8px; border:1px solid #3a5470; border-radius:8px; padding:8px 14px;",
        bs_icon("clipboard-data", size = "2em"),
        textOutput("n_samples", inline = TRUE), " samples"
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
    hr(style = "margin:1px 0; border-color:#3a5470;")
    ),
  
  navset_card_underline(
    title = "Exploration and Visualization",
    
    # 1. Overview
    nav_panel("Overview", layout_columns(col_widths = c(4, 4, 4), card(card_header("Tumor vs Normal Methylation"), plotOutput("overview_boxplot", height = "260px"), height = "310px"), card(card_header("Mutation status (KRAS / BRAF / TP53)"), plotOutput("overview_mutations", height = "260px"), height = "310px"), card(card_header("Sex / Stage / MSS composition"), plotOutput("overview_composition", height = "260px"), height = "310px"))),
    
    # 2. Genome Browser
    nav_panel("Genome Browser", div(style = "font-size:13px; color:#6c757d; margin-bottom:10px;", "Browse methylation across the whole genome: pick a chromosome, zoom in/out, or search by coordinates, bin ID, or gene name."), card(card_header("All chromosomes, click a point to jump to that chromosome"), plotlyOutput("genome_overview_plot", height = "420px")), card(card_header(textOutput("browser_position_header", inline = TRUE)), div(style = "display:flex; flex-wrap:wrap; gap:0px; align-items:center", div(style = "display:flex; gap:6px; align-items:center; flex:1 1 320px; min-width:280px;", textInput("browser_search", label = NULL, placeholder = "e.g. 12:25000000-26000000, 12_25000000, or KRAS", width = "100%"), actionButton("browser_search_go", NULL, icon = bs_icon("search"), class = "btn-sm btn-outline-secondary")), div(style = "display:flex; gap:4px; align-items:center;", actionButton("prev_chr", NULL, icon = bs_icon("chevron-left"), class = "btn-sm btn-outline-secondary"), selectInput("browser_chr", NULL, choices = chrom_list, selected = "1", width = "90px"), actionButton("next_chr", NULL, icon = bs_icon("chevron-right"), class = "btn-sm btn-outline-secondary")), div(style = "display:flex; gap:4px;", actionButton("zoom_in", NULL, icon = bs_icon("zoom-in"), class = "btn-sm btn-outline-secondary"), actionButton("zoom_out", NULL, icon = bs_icon("zoom-out"), class = "btn-sm btn-outline-secondary"), actionButton("zoom_reset", "Whole chromosome", icon = bs_icon("arrow-counterclockwise"), class = "btn-sm btn-outline-secondary"))), plotlyOutput("browser_chr_plot", height = "560px"), tags$div("Tip: bins currently in your sidebar selection are outlined on the plot above. You can also drag directly on the plot to zoom, and double-click it to reset.", style = "font-size:11px; color:#7d92a3; margin-top:6px;"))),
    
    # 3. Tumor vs Normal
    nav_panel("Tumor vs. Normal", layout_columns(col_widths = c(6, 6, 12), card(card_header("Methylation density (Tumor vs Normal)"), plotOutput("tn_density")), card(card_header("PCA / UMAP (interactive)"), radioButtons("proj_method", NULL, choices = c("PCA", "UMAP"), inline = TRUE), plotlyOutput("tn_projection", height = "400px")), card(card_header("Patient similarity network"), plotOutput("tn_network", height = "500px")))),
    
    # 4. Genome-wide Profile
    nav_panel("Genome-wide Profile", card(card_header("Manhattan-style plot: Tumor - Normal methylation difference"), plotlyOutput("manhattan_plot", height = "550px"))),
    
    # 5. Feature x Chromosome Heatmap
    nav_panel("Feature × Chromosome Heatmap", card(card_header("Diverging heatmap of selected features from metadata"), plotOutput("feature_heatmap", height = "600px"))),
    
    # 6. Clinical Explorer
    nav_panel("Clinical Explorer", layout_columns(col_widths = c(5, 7), card(card_header("Methylation by mutation status"), selectInput("mutation_gene", "Gene:", choices = c("KRAS", "BRAF", "TP53")), plotOutput("clinical_boxplot")), card(card_header("Clinical metadata table"), DTOutput("clinical_table")))),
    
    # 7. Bin Table
    nav_panel(
      "Bin Table",
      card(
        card_header("Filterable bin-level data"),
        downloadButton("bintable_download", "Download CSV"),
        DTOutput("bintable")
      )
    )
  ),
  
  div(
    style = "
    display:flex;
    align-items:center;
    justify-content:center;
    gap:18px;
    margin-top:24px;
    padding:16px 0;
    border-top:1px solid #dde3e8;
    color:#6c757d;
    font-size:13px;
  ", tags$span(
    "Paula Artiz Dueñas, Institut Germans Trias i Pujol (IGTP) - Universitat Politècnica de Catalunya (UPC)",
    style = "margin-left:10px;"
  )
  )
)

# Server

server <- function(input, output, session) {
  
  updateSelectizeInput(session, "bins", choices = bin_table$bin_id, server = TRUE)
  
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
  
  filtered_bin_table <- reactive({
    apply_bin_table_filters(selected_bin_table(), bin_table_filters, input)
  })
  
  output$n_samples <- renderText({ nrow(metadata) })
  output$n_patients <- renderText({ length(unique(metadata$patient_id)) })
  
  # 1. Overview
  
  output$overview_sex_comparison <- renderPlot({
  })
  
  output$overview_boxplot <- renderPlot({
  })
  
  output$overview_mutations <- renderPlot({
  })
  
  output$overview_composition <- renderPlot({
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
      nav$chr   <- res$chr # set first so the browser_chr observer above sees it and skips its reset
      nav$start <- max(1, floor(res$start))
      nav$end   <- min(chrom_lengths[[res$chr]], ceiling(res$end))
      if (nav$end <= nav$start) nav$end <- min(chrom_lengths[[res$chr]], nav$start + 2e6)
      if (!identical(input$browser_chr, res$chr)) updateSelectInput(session, "browser_chr", selected = res$chr)
      showNotification(res$message, type = "message", duration = 4)
    } else {
      showNotification(res$message, type = "error", duration = 5)
    }
  })
  
  # Clicking a point on the all-chromosomes overview jumps to that chromosome
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
    # Only highlight bins the user has explicitly picked in the sidebar
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
  output$tn_density <- renderPlot({
  })
  
  output$tn_projection <- renderPlotly({
  })
  
  output$tn_network <- renderPlot({
  })
  
  # 4. Genome-wide Profile
  output$manhattan_plot <- renderPlotly({
  })
  
  # 5. Feature × Chromosome Heatmap
  output$feature_heatmap <- renderPlot({
    
  })
  
  # 6. Clinical Explorer
  output$clinical_boxplot <- renderPlot({
  })
  
  output$clinical_table <- renderDT({
  })
  
  # 7. Bin Table
  output$bintable <- renderDT({
    datatable(selected_bin_table(),
              options = list(scrollX = TRUE, pageLength = 15),
              rownames = FALSE
    )
  })
  
  output$bintable_download <- downloadHandler(
    filename = function() {
      "bintable.csv"},
    content = function(file) {
      write.csv(selected_bin_table(), file, row.names = FALSE)
    })
  
}

shinyApp(ui = ui, server = server)
