# Exploration of Epigenomic Data in Colorectal Cancer
# Paula Artiz Dueñas - UPC

library(shiny)
library(bslib)
library(bsicons) # icons for value boxes
library(DT)
library(plotly)

data <- readRDS(
  "/Users/paulaartizduenas/Desktop/Internship - IGTP/Dataset/Data Processed/data_app.rds")

metadata <- data$metadata # 124 samples x 25 columns
methylation_long <- data$methylation_long
bin_table <- data$bin_table

chrom_list <- c(as.character(1:22), "X", "Y")

# Define UI
ui <- page_sidebar(
  title = div(
    style = "
    display: flex;
    justify-content: space-between;
    align-items: center;
    width: 100%;
  ",
    div(
      br(),
      h2("Exploration of Epigenomic Data in Colorectal Cancer"),
      p(
        "Paula Artiz Dueñas – Universitat Politècnica de Catalunya (UPC)",
        style = "margin-top:8px; font-size:16px; color:#6c757d;"
      )
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
    version      = 5,
    base_font    = font_google("Inter"),
    heading_font = font_google("Space Grotesk"),
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
    selectInput("chr", "Chromosome:", choices = chrom_list),
    
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
      layout_columns(
        col_widths = c(9, 3),
        card(
          card_header("Methylation along selected chromosome (Tumor vs Normal)"),
          plotlyOutput("browser_line", height = "500px")
        ),
        card(
          card_header("External links"), # configurar més endavant
          p("Configurar més endavant"),
          uiOutput("browser_links")
        )
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
        br(),
        br(),
        DTOutput("bintable")
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
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
  
  output$browser_links <- renderUI({
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
  })
  
  output$bintable_download <- downloadHandler(
    filename = function() {
      "bintable.csv"
    },
    content = function(file) {
      # Download logic to be implemented
    }
  )
  
}

# Run the application
shinyApp(ui = ui, server = server)
