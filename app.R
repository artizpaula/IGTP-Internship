# Exploration of Epigenomic Data in Colorectal Cancer
# Paula Artiz Dueñas - UPC

library(shiny)
library(bslib)
library(bsicons) # icons for value boxes
library(DT)
library(plotly)

data <- readRDS("/Users/paulaartizduenas/Desktop/Internship - IGTP/Dataset/Data Processed/data_app.rds")

metadata <- data$metadata # 124 samples x 25 columns
methylation_long <- data$methylation_long
bin_table <- data$bin_table

chrom_list <- c(as.character(1:22), "X", "Y")

# Ordered bin choices for the bin-selector (grouped by chromosome, sorted by position)
bin_table <- bin_table[order(match(bin_table$chr, chrom_list),bin_table$bin_position), ]
bin_choices_by_chr <- split(bin_table$bin_id, factor(bin_table$chr, levels = chrom_list))

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
    bg = "#f4f7f9",
    fg = "#0b2436",
    primary = "#0e7c86",
    secondary = "#16324f",
    success  = "#2fae66",
    info = "#3aa9c9",
    warning = "#e0a339",
    danger = "#d1495b",
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
    
    selectInput(
      "clinical_var",
      "Colour/group by:",
      choices = c(
        "sample_type", "sex", "stage", "mss_status",
        "KRAS", "BRAF", "TP53"
      )
    ),
    
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
      card(
        card_header("Methylation across selected bins (Tumor vs Normal)"),
        plotlyOutput("browser_line", height = "500px")
      )
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
  
  # Bins currently in scope for the analyses: whatever the user picked,
  # or (if none picked) every bin in the selected chromosome
  selected_bins <- reactive({
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
  output$browser_line <- renderPlotly({
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
