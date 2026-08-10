# Exploration of Epigenomic Data in Colorectal Cancer

##### Paula Artiz Dueñas - UPC

### Goal of the project

The main objective of this project is to develop an interactive web application using R Shiny to explore and visualize genomic bins. The application will enable users to analyze methylation patterns, compare tumor and normal samples, filter genomic regions, and explore genomic annotations through interactive visualizations.

### Functioning of the App

The app is organized into a sidebar for data selection and a set of tabs for exploration and visualization.

***Sidebar: Data Selection***
- Choose a chromosome and search/select specific genomic bins (e.g. `1_1000000`).
- Add all bins from the selected chromosome or clear the current selection.
- Upload a list of bins from a `.csv`, `.tsv`, or `.txt` file.

***Tabs: Exploration and Visualization***
1. **Overview**: Summary panels comparing tumor vs. normal methylation, mutation status, clinical stage/MSI-MSS composition, and sex distribution for the current bin selection.
2. **Genome Browser**: Browse methylation across the whole genome; jump to a chromosome, zoom in/out, or search by coordinates, bin ID, or gene name.
3. **Tumor vs. Normal**: Density plots, PCA/UMAP projection, and a patient similarity network based on methylation.
4. **Genome-wide Profile**: A Manhattan-style plot of the tumor–normal methylation difference across all bins, flagging statistically significant and outlier bins.
5. **Feature × Chromosome Heatmap**: A heatmap of methylation shift per patient and chromosome, with an optional colour strip for a selected clinical feature.
6. **Clinical Explorer**: Compare tumor methylation by mutation status (KRAS, BRAF, TP53) and browse/filter the full clinical metadata table, with patient selection linked to the plot.
7. **Bin Table**: A searchable, filterable table of all genomic bins with their methylation statistics and annotations, downloadable as a CSV.

### Necessary Datasets

The app loads a single pre-processed data object (`data_app.rds`) containing:
- Sample and patient metadata (clinical variables, mutation status, etc.).
- A per-bin table with methylation summary statistics and annotation columns.
- Long-format methylation values used for the exploratory plots.

This data object is generated upstream from the project's raw data, using the following processing scripts, which must be run **in this order** before the app is run:

1. `Alus_and_CpGs.R`
2. `Week_2_With_Prevalence.R`
3. `Gene_and_Functional_Annotation.R` -> Not developed yet

Running them in sequence produces `data_app.rds`, which is then loaded directly by `app_final_.R`.

### Necessary Files

- `app_final_.R`: Main Shiny application (UI + server logic).
- `data_app.rds`: Pre-processed dataset consumed by the app (path set at the top of the script and should be updated to match your local environment).
- Required R packages: `shiny`, `bslib`, `bsicons`, `DT`, `plotly`, `ggplot2`, `patchwork`.

### Running the App

1. Ensure the required R packages are installed.
2. Update the path to `data_app.rds` in `app_final_.R` if needed.
3. Run the app from R/RStudio with:
   
   ```r
   shiny::runApp("app_final_.R")
   ```
