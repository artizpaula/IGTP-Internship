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

---

### Processing Pipeline - Scripts and the Order They Must Run In

The app does **not** read raw data directly. It reads a single pre-processed object, `data_app.rds`, which is built by running three scripts, in this order:

| Order | Script | Role |
|---|---|---|
| 1 | `Alus_and_CpGs.R` | Defines helper functions that read raw per-CpG methylation tables and count CpGs/Alu elements per 1 Mb bin. Not run directly, it is `source()`-d by `Week_2_With_Prevalence.R`. |
| 2 | `Gene_Annotation.R` | Defines helper functions that overlap bin coordinates with the COSMIC Cancer Gene Census to annotate bins with gene names/IDs. Also `source()`-d by `Week_2_With_Prevalence.R`, not run directly. |
| 3 | `Week_2_With_Prevalence.R` | The actual driver script. Loads metadata, reshapes the per-sample bin JSON files, sources the two scripts above, computes CpG/Alu counts, gene annotation, prevalence statistics, and the final bin table, then writes out `data_app.rds`. |
| 4 | `app_final_.R` | The Shiny app itself. Loads `data_app.rds` (produced by step 3) and does not need to be "run" as a pipeline step, launch it with `shiny::runApp()`. |

You only ever execute `Week_2_With_Prevalence.R` yourself; it internally `source()`s `Alus_and_CpGs.R` and `Gene_Annotation.R`, so all three must simply be present in the same folder.

---

### Necessary Datasets (raw inputs you must supply)

None of the raw data is included in this repository, you must provide the following files/folders yourself, matching the exact formats below, before running `Week_2_With_Prevalence.R`.

1. **`Metadata_all_runs_combined.csv`** : clinical/sample metadata table.
   - Despite the `.csv` extension, it must be **tab-delimited** (read with `read.delim(..., sep = "\t")`).
   - Must contain at least a `sample_id` column and a `Sample2` column (used to derive `patient_id`).

2. **`counts_bins_norm_mean/`**: a folder containing one file per sample named `counts_<sample_id>.txt`.
   - Each file is a **JSON** object nested as `{chromosome: {bin_position: methylation_value}}`.
   - These are parsed and reshaped into the long-format methylation table (`Methylation_long.csv`).

3. **Raw per-CpG methylation tables** (used only for CpG/Alu counting), either:
   - loose `.tsv`/`.txt` files, or
   - `.tar.gz` / `_tar.gz` / `.tgz` archives containing such files,
   - each with (at minimum) tab-separated columns `chr`, `end`, `CpG`, `alu_region`, and optionally `meth`/`meth_pct`.
   - These live in a dedicated folder (referred to as `dataset_alus_cpgs` in the script) that is **separate** from the bins folder above.

4. **COSMIC Cancer Gene Census reference**, used for gene annotation:
   - File: `Cosmic_CancerGeneCensus_v101_GRCh37.tsv`
   - Inside a subfolder named `Cosmic_CancerGeneCensus_Tsv_v101_GRCh37/`
   - Must contain the columns `GENE_SYMBOL`, `COSMIC_GENE_ID`, `CHROMOSOME`, `GENOME_START`, `GENOME_STOP`.
   - This TSV must be downloaded separately from COSMIC (https://cancer.sanger.ac.uk/census), it is **not** redistributed here due to licensing.

5. **`CRC_curated_genes.txt`** (optional), a plain text list of curated colorectal-cancer gene symbols, one per line.
   - The path is defined (`crc_gene_list_path`) but **is currently not passed** into `read_gene_reference()` in `Week_2_With_Prevalence.R` (line ~192 calls `read_gene_reference(cosmic_tsv_path)` only). As shipped, the pipeline annotates bins against the **full** COSMIC census, not a curated subset. You can safely omit this file unless you edit the script to pass it in.

#### Suggested local folder layout

Any layout works as long as the path variables (see below) are updated to match it, but a layout mirroring the original project looks like:

```
Project/
├── R Scripts/
│   ├── Alus_and_CpGs.R
│   ├── Gene_Annotation.R
│   ├── Week_2_With_Prevalence.R
│   └── app_final_.R
└── Dataset/
    ├── Metadata/
    │   └── Metadata_all_runs_combined.csv
    ├── Bins/
    │   └── counts_bins_norm_mean/
    │       ├── counts_<sample1>.txt
    │       └── ...
    ├── Alus_CpGs/
    │   ├── <raw methylation files or .tar.gz archives>
    └── Archivo/
        ├── Cosmic_CancerGeneCensus_Tsv_v101_GRCh37/
        │   └── Cosmic_CancerGeneCensus_v101_GRCh37.tsv
        └── CRC_curated_genes.txt   # optional
```

Running the pipeline will also create a `Data Processed/` folder (a sibling of `Dataset/Metadata/`, i.e. `Project/Dataset/Data Processed/`) containing all intermediate CSVs and the final `data_app.rds`.

---

### Paths That Must Be Edited Before Running

All of the scripts as provided contain **hardcoded, machine-specific paths** from the original author's computer. These will not work on any other machine and must be changed.

#### `Week_2_With_Prevalence.R`
- **Line 7**: `setwd("/Users/paulaartizduenas/Desktop/Project/R Scripts")`
  → Change to the folder on your machine containing `Alus_and_CpGs.R`, `Gene_Annotation.R`, and `Week_2_With_Prevalence.R` itself (they must all be in the same folder, since they're `source()`-d by relative path on lines 9 and 11).
- **Lines 15–18**: the four dataset root paths:
  ```r
  dataset_metadata <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Metadata"
  dataset_bins <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Bins"
  dataset_alus_cpgs <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Alus_CpGs"
  dataset_gene_annotation <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Archivo"
  ```
  → All four **must** be updated to point to your local copies of the datasets described above. Everything else in the script (`metadata_dir`, `bins_dir`, `output_dir`, `cosmic_tsv_path`, `crc_gene_list_path`) is built automatically from these four with `file.path()`, so you do not need to edit anything past line 28, as long as your folder/file names inside each of the four roots match exactly what's listed in the "Necessary Datasets" section above (e.g. the metadata file must literally be named `Metadata_all_runs_combined.csv` inside `dataset_metadata`).
- No other paths need editing in this file. `output_dir` is auto-created (`dir.create(..., recursive = TRUE)`), so it doesn't need to pre-exist.

#### `Alus_and_CpGs.R` and `Gene_Annotation.R`
- No hardcoded paths, every function takes a `filepath`/`cosmic_tsv_path` argument. Nothing to edit here; they just need to sit next to `Week_2_With_Prevalence.R` so the `source()` calls in step above can find them.

#### `app_final_.R`
- **Line 12**: `data <- readRDS("data_app.rds")`
  → This is a relative path, so `data_app.rds` (produced by `Week_2_With_Prevalence.R`, found in the `Data Processed/` output folder) must either be **copied into the same folder as `app_final_.R`**, or you must change this line to the correct full/relative path where you keep it, e.g.:
  ```r
  data <- readRDS("/path/to/Data Processed/data_app.rds")
  ```
- **Lines 58–66**: fallback lookup for the COSMIC census TSV and `Gene_Annotation.R`, used **only** if `data_app.rds` is somehow missing gene annotation columns (it won't be, if you ran the full pipeline above successfully):
  ```r
  cosmic_tsv_path <- find_first_existing(c(
    "Archivo/Cosmic_CancerGeneCensus_Tsv_v101_GRCh37/Cosmic_CancerGeneCensus_v101_GRCh37.tsv",
    "../Archivo/Cosmic_CancerGeneCensus_Tsv_v101_GRCh37/Cosmic_CancerGeneCensus_v101_GRCh37.tsv",
    "Data/Archivo/Cosmic_CancerGeneCensus_Tsv_v101_GRCh37/Cosmic_CancerGeneCensus_v101_GRCh37.tsv"
  ))
  gene_annotation_script <- find_first_existing(c(
    "Gene_Annotation.R",
    "../Gene_Annotation.R"
  ))
  ```
  → This is a safety net, not required for normal use. You only need to touch it if `data_app.rds` lacks the `gene_count`/`gene_ids`/`gene_names` columns (e.g. you ran the pipeline without the COSMIC TSV present) and want the app to recompute gene annotation on load instead. In that case, place a copy of the `Archivo/` folder and `Gene_Annotation.R` in one of the three candidate relative locations above, next to `app_final_.R`.
- No other paths in `app_final_.R` need editing. The `fileInput` for uploading a custom bin list (`.csv`/`.tsv`/`.txt`) is handled entirely through Shiny's upload widget and does not require a path edit.

---

### Necessary Files (summary)

**Scripts (included in this repo):**
- `Alus_and_CpGs.R`: CpG/Alu counting helper functions (sourced by `Week_2_With_Prevalence.R`).
- `Gene_Annotation.R`: gene annotation helper functions (sourced by `Week_2_With_Prevalence.R`).
- `Week_2_With_Prevalence.R`: pipeline driver script; produces `data_app.rds`.
- `app_final_.R`: main Shiny application (UI + server logic).

**Data you must supply (not included, see "Necessary Datasets" above):**
- `Metadata_all_runs_combined.csv`
- `counts_bins_norm_mean/` (per-sample JSON bin files)
- Raw per-CpG methylation files/archives (for CpG/Alu counting)
- `Cosmic_CancerGeneCensus_v101_GRCh37.tsv` (COSMIC Cancer Gene Census, downloaded separately)
- `CRC_curated_genes.txt` (optional, currently unused by default)

**Generated by running the pipeline:**
- `data_app.rds`: pre-processed dataset consumed by the app (also several intermediate CSVs: `Metadata_clean.csv`, `Methylation_long.csv`, `CpG_Alu_bin_annotation.csv`, `Gene_bin_overlaps.csv`, `Bin_annotation_template.csv`, `Bin_prevalence_detection.csv`, `Bin_table.csv`).

**Required R packages:**
- For the pipeline (`Week_2_With_Prevalence.R`, which sources `Alus_and_CpGs.R` and `Gene_Annotation.R`): `jsonlite`, `data.table`.
- For the app (`app_final_.R`): `shiny`, `bslib`, `bsicons`, `DT`, `plotly`, `ggplot2`, `patchwork`.

---

### Running the Pipeline and App

1. Install the required R packages listed above.
2. Place `Alus_and_CpGs.R`, `Gene_Annotation.R`, and `Week_2_With_Prevalence.R` in the same folder.
3. Edit the paths in `Week_2_With_Prevalence.R` as described above (line 7 and lines 15–18) to point to your local script folder and dataset folders.
4. Run the pipeline:
   ```r
   source("Week_2_With_Prevalence.R")
   ```
   This produces `data_app.rds` inside a `Data Processed/` folder next to your metadata folder.
5. Copy (or symlink) `data_app.rds` into the same folder as `app_final_.R`, or edit line 12 of `app_final_.R` to point directly to it.
6. Run the app from R/RStudio with:
   ```r
   shiny::runApp("app_final_.R")
   ```
