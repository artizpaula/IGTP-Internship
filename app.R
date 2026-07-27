# Exploration of Epigenomic Data in Colorectal Cancer
# Paula Artiz Dueñas - UPC

library(shiny)
library(bslib)
library(bsicons) # icons for value boxes
library(DT)
library(plotly)

data <- readRDS("/Users/paulaartizduenas/Desktop/Internship - IGTP/Dataset/Data Processed/data_app.rds")

bin_table <- data$bin_table
chrom_list <- c(as.character(1:22), "X", "Y")

# Ordered bin choices for the bin-selector (grouped by chromosome, sorted by position)
bin_table <- bin_table[order(match(bin_table$chr, chrom_list),bin_table$bin_position), ]
bin_choices_by_chr <- split(bin_table$bin_id, factor(bin_table$chr, levels = chrom_list))

# Length of each chromosome, inferred from the last bin on it (bin_end)
chrom_lengths <- tapply(bin_table$bin_end, bin_table$chr, max)
chrom_lengths <- setNames(as.numeric(chrom_lengths[chrom_list]), chrom_list)

# Cumulative genome-wide offset so all chromosomes can sit end-to-end on one axis
chr_offset <- setNames(numeric(length(chrom_list)), chrom_list)
running <- 0
for (chr in chrom_list) {
  chr_offset[chr] <- running
  running <- running + chrom_lengths[[chr]]}
total_genome_len <- running
chr_mid <- chr_offset + chrom_lengths / 2

# Genome-wide bin table
genome_wide_bins <- bin_table
genome_wide_bins$chr <- factor(genome_wide_bins$chr, levels = chrom_list)
genome_wide_bins <- genome_wide_bins[order(genome_wide_bins$chr, genome_wide_bins$bin_position), ] # bin's within-chromosome position -> genome-wide x-coordinate
genome_wide_bins$genome_x <- chr_offset[as.character(genome_wide_bins$chr)] + genome_wide_bins$bin_position
genome_wide_bins$overall_meth <- rowMeans(
  genome_wide_bins[, c("mean_methylation_tumor", "mean_methylation_normal")], na.rm = TRUE)

# Dashed vertical lines marking chromosome boundaries, reused on both plots
chr_boundary_shapes <- lapply(as.numeric(chr_offset)[-1], function(x) {
  list(type = "line", x0 = x, x1 = x, y0 = 0, y1 = 1, yref = "paper",
       line = list(color = "rgba(150,150,150,0.4)", width = 1, dash = "dot"))
})

# Print base-pair distance
format_bp <- function(bp) {
  bp <- abs(bp)
  if (bp >= 1e6) sprintf("%.2f Mb", bp / 1e6)
  else if (bp >= 1e3) sprintf("%.0f kb", bp / 1e3)
  else sprintf("%.0f bp", bp)
}

# Static lookup table of hg19/UCSC coordinates for a few well-known cancer genes.
# Used only as a fallback when a searched gene isn't already annotated in bin_table$genes.
gene_not_found <- data.frame(
  gene  = c("APC", "PIK3CA", "SMAD4", "PTEN",
            "MLH1", "MSH2", "MSH6", "FBXW7"),
  chr   = c("5", "3", "18", "10", "3", "2", "2", "4"),
  start = c(112073582, 178866145, 48556583, 89623382,
            37035009, 47630295, 48010284, 153241696),
  end   = c(112181936, 178957881, 48611412, 89731687,
            37092337, 47710362, 48034092, 153457244),
  stringsAsFactors = FALSE
) 

# Parses a free-text genome browser search box into a (chr, start, end) target
parse_genomic_search <- function(query, bin_table, chrom_list, gene_not_found, pad = 2e6) {
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
  in_annotation <- bin_table[has_annotation, ]
  if (nrow(in_annotation) >= 1) {
    return(list(status = "ok", chr = in_annotation$chr[1],
                start = min(in_annotation$bin_start), end = max(in_annotation$bin_end),
                message = paste0("Jumped to gene ", gene_query)))
  }
  
  # case 4: gene symbol not in bin_table, fall back to the static gene_not_found lookup
  in_lookup <- gene_not_found[toupper(gene_not_found$gene) == gene_query, ]
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

# Define UI
ui <- page_sidebar(
  title = div(
    style = "
    display: flex;
    justify-content: space-between;
    align-items: center;
    width: 100%;
    gap: 1px;
    padding: 22px 5px;
    border-bottom: 1px solid rgba(255,255,255,0.08);
  ",
    # Left block: title
    div(
      style = "
      flex:1 1 auto;
      min-width:0;
    ",
      h2("Exploration of Epigenomic Data in Colorectal Cancer", style = "margin:0; font-weight:600; letter-spacing:0.2px; white-space:normal;")
    ),
    div(
      style = "
      display:flex;
      gap:14px;
      align-items:center;
      margin-left:auto;
    ",
      div(
        style = "
        display:flex;
        align-items:center;
        gap:8px;
        border:1px solid #3a5470;
        border-radius:8px;
        padding:8px 14px;
      ",
        bs_icon("person", size = "2em"),
        textOutput("n_patients", inline = TRUE),
        " patients"
      ),
      
      div(
        style = "
        display:flex;
        align-items:center;
        gap:8px;
        border:1px solid #3a5470;
        border-radius:8px;
        padding:8px 14px;
      ",
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
    
    tags$div(
      "DATA SELECTION",
      style = "font-size:14px; font-weight:700; letter-spacing:0.8px; color:#7d92a3; margin-bottom:0px;"
    ),
    selectInput("chr", "Chromosome:", choices = chrom_list),
    
    selectizeInput(
      "bins",
      "Selected bins:",
      choices = NULL,
      multiple = TRUE,
      options = list(
        placeholder = "Find and select bins (ex: 1_1000000)...",
        plugins = list("remove_button")
      )
    ),
    
    actionButton(
      "add_chr_bins",
      "Add all bins on the chromosome",
      icon = bs_icon("plus-circle"),
      class = "btn-sm class=btn-outline-light w-100"),
    
    actionButton(
      "clear_bins",
      "Clear bin selection",
      icon = bs_icon("x-circle"),
      class = "btn-sm class=btn-outline-light w-100"
    ),
    
    hr(style = "margin:1px 0; border-color:#3a5470;"),
    
    tags$div(
      "DISPLAY OPTIONS",
      style = "font-size:14px; font-weight:700; letter-spacing:0.8px; color:#7d92a3; margin-bottom:0px;"
    ),
    
    selectInput("clinical_var", "Colour/group by:", choices = c(
      "Sample type"    = "Type",
      "Sex"            = "sexe_label",
      "Stage"          = "estadi2",
      "MSI/MSS status" = "MSS",
      "KRAS status"    = "KRAS",
      "BRAF status"    = "BRAF",
      "TP53 status"    = "TP53"
    )),
    
    selectInput(
      "feature_type",
      "Genomic feature:",
      choices = c(
        "All", "Promoter", "Gene body",
        "Intergenic", "Enhancer"
      )
    )
  ),
  
  navset_card_underline(
    title = "Exploration and Visualization",
    
    # 1. Overview
    nav_panel(
      "Overview",
      layout_columns(
        col_widths = c(4, 4, 4),
        card(
          card_header("Tumor vs Normal Methylation"),
          plotOutput("overview_boxplot", height = "260px"),
          height = "310px"
        ),
        card(
          card_header("Mutation status (KRAS / BRAF / TP53)"),
          plotOutput("overview_mutations", height = "260px"),
          height = "310px"
        ),
        card(
          card_header("Sex / Stage / MSS composition"),
          plotOutput("overview_composition", height = "260px"),
          height = "310px"
        )
      )
    ),
    
    # 2. Genome Browser
    nav_panel(
      "Genome Browser",
      div(
        style = "font-size:13px; color:#6c757d; margin-bottom:10px;",
        "Browse methylation across the whole genome: pick a chromosome, zoom in/out, ",
        "or search by coordinates, bin ID, or gene name."
      ),
      
      card(
        card_header(textOutput("browser_position_header", inline = TRUE)),
        
        div(
          style = "display:flex; flex-wrap:wrap; gap:10px; align-items:center; margin-bottom:12px;",
          
          div(
            style = "display:flex; gap:6px; align-items:center; flex:1 1 320px; min-width:280px;",
            textInput("browser_search", label = NULL,
                      placeholder = "e.g. 12:25000000-26000000, 12_25000000, or KRAS", width = "100%"),
            actionButton("browser_search_go", NULL, icon = bs_icon("search"), class = "btn-sm btn-outline-secondary")
          ),
          
          div(
            style = "display:flex; gap:4px; align-items:center;",
            actionButton("prev_chr", NULL, icon = bs_icon("chevron-left"), class = "btn-sm btn-outline-secondary"),
            selectInput("browser_chr", NULL, choices = chrom_list, selected = "1", width = "90px"),
            actionButton("next_chr", NULL, icon = bs_icon("chevron-right"), class = "btn-sm btn-outline-secondary")
          ),
          
          div(
            style = "display:flex; gap:4px;",
            actionButton("zoom_in", NULL, icon = bs_icon("zoom-in"), class = "btn-sm btn-outline-secondary"),
            actionButton("zoom_out", NULL, icon = bs_icon("zoom-out"), class = "btn-sm btn-outline-secondary"),
            actionButton("zoom_reset", "Whole chromosome", icon = bs_icon("arrow-counterclockwise"),
                         class = "btn-sm btn-outline-secondary")
          )
        ),
        
        plotlyOutput("browser_chr_plot", height = "460px"),
        
        tags$div(
          "Tip: bins currently in your sidebar selection are outlined on the plot above. ",
          "You can also drag directly on the plot to zoom, and double-click it to reset.",
          style = "font-size:11px; color:#7d92a3; margin-top:6px;"
        )
      ),
      
      # Press Enter in the search box to trigger the same action as the Go button
      tags$script(HTML(
        "$(document).on('keydown', '#browser_search', function(e) {
           if (e.key === 'Enter') { e.preventDefault(); $('#browser_search_go').trigger('click'); }
         });"
      ))
    ),
    
    # 3. Tumor vs. Normal
    nav_panel(
      "Tumor vs. Normal",
      layout_columns(
        col_widths = c(6, 6, 12),
        card(
          card_header("Methylation density (Tumor vs Normal)"),
          plotOutput("tn_density")
        ),
        card(
          card_header("PCA / UMAP (interactive)"),
          radioButtons(
            "proj_method",
            NULL,
            choices = c("PCA", "UMAP"),
            inline = TRUE
          ),
          plotlyOutput("tn_projection", height = "400px")
        ),
        card(
          card_header("Patient similarity network"),
          plotOutput("tn_network", height = "500px")
        )
      )
    ),
    
    # 4. Genome-wide Profile
    nav_panel(
      "Genome-wide Profile",
      card(
        card_header("Manhattan-style plot: Tumor - Normal methylation difference"),
        plotlyOutput("manhattan_plot", height = "550px")
      )
    ),
    
    # 5. Feature × Chromosome Heatmap
    nav_panel(
      "Feature × Chromosome Heatmap",
      card(
        card_header("Diverging heatmap of methylation shift (feature × chromosome)"),
        plotOutput("feature_heatmap", height = "600px")
      )
    ),
    
    # 6. Clinical Explorer
    nav_panel(
      "Clinical Explorer",
      layout_columns(
        col_widths = c(5, 7),
        card(
          card_header("Methylation by mutation status"),
          selectInput(
            "mutation_gene",
            "Gene:",
            choices = c("KRAS", "BRAF", "TP53")
          ),
          plotOutput("clinical_boxplot")
        ),
        card(
          card_header("Clinical metadata table"),
          DTOutput("clinical_table")
        )
      )
    ),
    
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
  
  # Bin selection
  
  # Populate the bin selection (server-side)
  updateSelectizeInput(
    session, "bins",
    choices = bin_table$bin_id,
    server = TRUE
  )
  
  # Button: add all bins of the currently selected chromosome to the selection
  observeEvent(input$add_chr_bins, {
    chr_bins <- bin_choices_by_chr[[input$chr]]
    updated <- union(input$bins, chr_bins)
    updateSelectizeInput(
      session, "bins",
      choices = bin_table$bin_id,
      selected = updated,
      server = TRUE
    )
  })
  
  # Button: clear bin selection
  observeEvent(input$clear_bins, {
    updateSelectizeInput(
      session, "bins",
      choices = bin_table$bin_id,
      selected = character(0),
      server = TRUE
    )
  })
  
  # Bins currently in scope for the analyses: whatever the user picked,or (if none picked) every bin in the selected chromosome
  selected_bins <- reactive({ # use the bins explicitly selected by the user
    if (length(input$bins) > 0) {
      input$bins
    } else {
      bin_choices_by_chr[[input$chr]]
    }
  })
  
  # Subset of bin_table matching the current bin selection
  selected_bin_table <- reactive({
    bin_table[bin_table$bin_id %in% selected_bins(), ]
  })
  
  # 1. Overview
  output$n_samples <- renderText({
    nrow(metadata)
  })
  
  output$n_patients <- renderText({
    length(unique(metadata$patient_id))
  })
  
  output$overview_sex_comparison <- renderPlot({
  })
  
  output$overview_boxplot <- renderPlot({
  })
  
  output$overview_mutations <- renderPlot({
  })
  
  output$overview_composition <- renderPlot({
  })
  
  
  # 2. Genome Browser

  observeEvent(input$chr, {
    updateSelectInput(session, "browser_chr", selected = input$chr)
  })
  
  output$genome_browser <- renderPlot({
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
    }
  )
  
}

# Run the application
shinyApp(ui = ui, server = server)
