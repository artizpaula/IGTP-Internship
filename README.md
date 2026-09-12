# Exploration of Epigenomic Data in Colorectal Cancer

##### Paula Artiz Dueñas, UPC

- **App Link: https://paulaartiz.shinyapps.io/app_igtp/**
- **App Alus Link: https://paulaartiz.shinyapps.io/app_igtp_alu/**

### What This Project Is About

This project is an interactive **R Shiny web app** for exploring and visualizing DNA methylation data in colorectal cancer. It lets users look at genomic "bins" (1 Mb chunks of the genome), compare methylation between tumor and normal samples, filter by region, and explore gene and promoter annotations, all through interactive plots and tables.

There are actually **two apps** built on top of the same pipeline:

| App | Reads | What it adds |
|---|---|---|
| `app.R` | `data_app.rds` | The core bin-level exploration described below. |
| `app_alus.R` | `data_app.rds` **and** `data_app_alus.rds` | Everything in `app.R`, plus an **Alu Table** tab (one row per individual Alu element instead of per 1 Mb bin) and a Genome Browser that can jump straight to, and highlight, whichever Alus you've selected. |

`app_alus.R` needs an extra preprocessing step (`Alu_Table_Preprocessing_alus.R`) to build `data_app_alus.rds` — see the pipeline table below.

### What the App Can Do

The app has a sidebar for picking data, plus several tabs for exploring it.

**Sidebar, choosing your data:**
- Pick a chromosome and search/select specific bins (e.g. `1_1000000`).
- Add every bin from a chromosome at once, or clear your selection.
- Upload your own list of bins from a `.csv`, `.tsv`, or `.txt` file.
- *(`app_alus.R` only)* Do the same thing for individual Alu elements (e.g. `1:51584-51880`): select them, add a whole chromosome's worth, upload a list, or jump the Genome Browser straight to your current Alu selection with one click.

**Tabs, exploring and visualizing:**
1. **Overview**: summary panels comparing tumor vs. normal methylation, mutation status, clinical stage/MSI-MSS, and sex distribution for the bins you selected.
2. **Genome Browser**: scroll through methylation across the whole genome; jump to a chromosome, zoom, or search by coordinates, bin ID, or gene name. *(`app_alus.R` only)* also draws a track of individual Alu elements along the bottom of the plot, highlighting whichever ones are in your sidebar selection, and offers a toggle to hide the track.
3. **Tumor vs. Normal**: density plots, PCA/UMAP, and a patient similarity network based on methylation.
4. **Genome-wide Profile**: a Manhattan-style plot of tumor–normal methylation differences, flagging significant/outlier bins.
5. **Feature × Chromosome Heatmap**: methylation shift per patient and chromosome, with an optional color strip for a clinical feature.
6. **Clinical Explorer**: compare tumor methylation by mutation status (KRAS, BRAF, TP53) and browse/filter the clinical metadata table.
7. **Bin Table**: a searchable, filterable table of all bins with their stats and annotations (including gene and promoter overlaps), downloadable as a CSV.
8. **Alu Table** *(`app_alus.R` only)*: the same idea as the Bin Table, one row per individual Alu element instead of per 1 Mb bin. Each row has a Quick Links column with a one-click jump into the in-app Genome Browser plus external UCSC/Ensembl/TCGA/NCBI links, and the table is filterable/downloadable the same way as the Bin Table.

---

### How the Pipeline Works (Run These in Order!)

Neither app touches the raw data directly, each just loads pre-built `.rds` file(s). Those files are created by running the scripts below, in order.

| Step | Script | What it does |
|---|---|---|
| 1 | `Alus_and_CpGs.R` | Helper functions that read the raw per-CpG methylation tables and count CpGs/Alu elements per 1 Mb bin. You don't run this yourself, it gets `source()`-d automatically by `Data_Preprocessing.R`. |
| 2 | `Gene_Annotation.R` | Helper functions that match coordinates against the COSMIC Cancer Gene Census to add gene names/IDs. Also `source()`-d automatically, not run directly (used by both `Data_Preprocessing.R` and `Alu_Table_Preprocessing_alus.R`). |
| 3 | `Promoter_Annotation.R` | Helper functions that match coordinates against a promoter reference file to add promoter names/counts. Also `source()`-d automatically, not run directly (used by both `Data_Preprocessing.R` and `Alu_Table_Preprocessing_alus.R`). |
| 4 | `Data_Preprocessing.R` | **Run this first.** It loads the metadata, reshapes the per-sample JSON files, sources the three scripts above, computes CpG/Alu counts, gene annotation, and promoter annotation, works out prevalence stats, builds the final bin table, and saves `data_app.rds`. |
| 5 | `Alu_Table_Preprocessing_alus.R` | **Run this second, only if you want to use `app_alus.R`.** Reads the raw per-Alu methylation matrix, reshapes it, pulls Tumor/Normal + patient info in from `data_app.rds` (so step 4 must have already run), assigns each Alu to its parent 1 Mb bin, adds gene/promoter annotation (directly from COSMIC/promoter files if available, otherwise inherited from the parent bin's annotation already sitting in `data_app.rds`), and saves `data_app_alus.rds`. |
| 6 | `app.R` / `app_alus.R` | The Shiny apps themselves. `app.R` just loads `data_app.rds` (from step 4); `app_alus.R` loads both `data_app.rds` (step 4) and `data_app_alus.rds` (step 5). Launch either with `shiny::runApp()`; they're not part of the "pipeline" as such. |

**For `app.R`, you only need to run `Data_Preprocessing.R`. For `app_alus.R`, run `Data_Preprocessing.R` first, then `Alu_Table_Preprocessing_alus.R`. As long as all the `.R` scripts are sitting in the same folder, each driver script takes care of sourcing its own helpers.**

---

### Data You Need to Provide

None of the raw data is included in this repo, you need to gather the following yourself, matching these formats exactly, before running the pipeline.

1. **`Metadata_all_runs_combined.csv`**: the clinical/sample metadata table.
   - Despite the `.csv` name, it's actually **tab-delimited** (read with `read.delim(..., sep = "\t")`).
   - Needs at least a `sample_id` column and a `Sample2` column (used to build `patient_id`).

2. **`counts_bins_norm_mean/`**: a folder with one file per sample, named `counts_<sample_id>.txt`.
   - Each file is a JSON object shaped like `{chromosome: {bin_position: methylation_value}}`.
   - These get parsed into the long-format methylation table (`Methylation_long.csv`).

3. **Raw per-CpG methylation tables** (only needed for CpG/Alu counting):
   - Loose `.tsv`/`.txt` files, or `.tar.gz` / `_tar.gz` / `.tgz` archives containing such files.
   - Keep these in their own folder (called `dataset_alus_cpgs` in the script), separate from the bins folder above.

4. **COSMIC Cancer Gene Census**: used for gene annotation:
   - File: `Cosmic_CancerGeneCensus_v101_GRCh37.tsv`
   - Inside a subfolder called `Cosmic_CancerGeneCensus_Tsv_v101_GRCh37/`
   - Needs the columns `GENE_SYMBOL`, `COSMIC_GENE_ID`, `CHROMOSOME`, `GENOME_START`, `GENOME_STOP`.
   - Also used by `Alu_Table_Preprocessing_alus.R` if you want gene annotation computed directly on each Alu rather than inherited from its parent bin (see item 6 below).

5. **Promoter reference file**: used for promoter annotation:
   - File: `Promoter_reference_GRCh37.tsv`
   - Sits directly inside the `Archivo` folder (same folder as the COSMIC subfolder above).
   - Needs the columns `promoter_name`, `chr`, `start`, `end` (GRCh37 coordinates).
   - Not included in this repo, you need to build/export it yourself, e.g. from the UCSC Table Browser (`refGene` track, TSS ± N bp) or the Ensembl Regulatory Build filtered to `feature_type == "Promoter"`.
   - If this file is missing, the pipeline still runs, it just fills `promoter_count`/`promoter_names` as empty and the app shows a warning banner on the Bin Table (and Alu Table) tab.

6. **`alu_methylation_matrix_all_runs_comb_norm.tsv`**: *(only needed for `app_alus.R`)* the raw per-Alu methylation matrix, used by `Alu_Table_Preprocessing_alus.R`.
   - Wide format, tab-delimited: first column is the Alu identifier as `chr:start-end` (e.g. `1:51584-51880`), every other column is one sample's methylation value at that Alu, with the column header matching a `sample_id` in your metadata.
   - `Alu_Table_Preprocessing_alus.R` reshapes this to long format internally and joins in Tumor/Normal + patient info from `data_app.rds`, so `Data_Preprocessing.R` must be run first.
   - If `Cosmic_CancerGeneCensus_v101_GRCh37.tsv` and `Promoter_reference_GRCh37.tsv` (item 4 & 5) aren't found next to the script, gene/promoter annotation is inherited from each Alu's parent 1 Mb bin instead (already computed in `data_app.rds`) rather than being skipped entirely.

#### Suggested folder layout

Any layout works as long as you update the path variables to match (see below), but here's a layout that mirrors the original project:

```
Project/
├── R Scripts/
│   ├── Alus_and_CpGs.R
│   ├── Gene_Annotation.R
│   ├── Promoter_Annotation.R
│   ├── Data_Preprocessing.R
│   ├── Alu_Table_Preprocessing_alus.R
│   ├── app.R
│   └── app_alus.R
└── Dataset/
    ├── Metadata/
    │   └── Metadata_all_runs_combined.csv
    ├── Bins/
    │   └── counts_bins_norm_mean/
    │       ├── counts_<sample1>.txt
    │       └── ...
    ├── Alus_CpGs/
    │   ├── <raw methylation files or .tar.gz archives>
    ├── Alu_Matrix/
    │   └── alu_methylation_matrix_all_runs_comb_norm.tsv
    └── Archivo/
        ├── Cosmic_CancerGeneCensus_Tsv_v101_GRCh37/
        │   └── Cosmic_CancerGeneCensus_v101_GRCh37.tsv
        ├── Promoter_reference_GRCh37.tsv
        └── CRC_curated_genes.txt   # optional
```

Running `Data_Preprocessing.R` creates a `Data Processed/` folder (next to `Dataset/Metadata/`, i.e. `Project/Dataset/Data Processed/`) with all the intermediate CSVs plus `data_app.rds`. Running `Alu_Table_Preprocessing_alus.R` afterwards adds `data_app_alus.rds` alongside it (by default it saves to the working directory, so point `output_path` there too if you want it in the same place — see below).

---

### Paths You Need to Edit Before Running

The scripts as provided still have **hardcoded paths from the original author's laptop**. These won't work anywhere else, so you'll need to change them.

#### In `Data_Preprocessing.R`

- **Line 7**: `setwd("/Users/paulaartizduenas/Desktop/Project/R Scripts")`
  → Change this to wherever `Alus_and_CpGs.R`, `Gene_Annotation.R`, `Promoter_Annotation.R`, and `Data_Preprocessing.R` live on your machine (they all need to be in the same folder, since the `source()` calls near the top find them by relative path).

- **Lines 15–18**: the four dataset root paths (Update **all four** to point to your local copies of the data described above).
  ```r
  dataset_metadata <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Metadata"
  dataset_bins <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Bins"
  dataset_alus_cpgs <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Alus_CpGs"
  dataset_gene_annotation <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Archivo"
  ```

- **Lines 58–61**: Change this to wherever "Cosmic_CancerGeneCensus_v101_GRCh37.tsv" file is located in your machine.
  ```r
  cosmic_tsv_path <-find_first_existing(c("/Users/paulaartizduenas/Desktop/Project/Dataset/Archivo/Cosmic_CancerGeneCensus_Tsv_v101_GRCh37/Cosmic_CancerGeneCensus_v101_GRCh37.tsv"))

  ```

- **`promoter_tsv_path`**: by default this points at `Promoter_reference_GRCh37.tsv` inside `dataset_gene_annotation` (i.e. your `Archivo` folder), so it updates automatically once you fix line 15–18 above. Only change it separately if you name/place the file differently.
  ```r
  promoter_tsv_path <- file.path(dataset_gene_annotation, "Promoter_reference_GRCh37.tsv")
  ```

#### In `Alu_Table_Preprocessing_alus.R` *(only needed for `app_alus.R`)*

Unlike `Data_Preprocessing.R`, this script has no `setwd()` call and uses plain relative filenames throughout (lines 6–10), so it assumes it's being run with its working directory already set to a folder containing all of the following. Change each one to a full path if that's not the case for you:

```r
alu_matrix_path   <- "alu_methylation_matrix_all_runs_comb_norm.tsv"  # item 6 above
data_app_path     <- "data_app.rds"                                  # output of Data_Preprocessing.R (step 4)
output_path       <- "data_app_alus.rds"                              # where this script's output gets saved
cosmic_tsv_path   <- "Cosmic_CancerGeneCensus_v101_GRCh37.tsv"        # optional, same file as item 4 above
promoter_tsv_path <- "Promoter_reference_GRCh37.tsv"                  # optional, same file as item 5 above
```

Also make sure `Gene_Annotation.R` and `Promoter_Annotation.R` are in the same folder as this script (or on the `source()` search path), since it sources both at the top, same as `Data_Preprocessing.R` does.

Note `cosmic_tsv_path` and `promoter_tsv_path` are optional here: if either file isn't found, that annotation is inherited from each Alu's parent bin (already computed by `Data_Preprocessing.R`) instead of being computed fresh, so you can leave them pointing at nonexistent files and the script will still run.

#### In `app.R`

- **Line 12**: `data <- readRDS("data_app.rds")`
  → Since this is a relative path, either copy `data_app.rds` (from the pipeline's `Data Processed/` output) into the same folder as `app.R`, or change this line to point to it directly, e.g.:
  ```r
  data <- readRDS("/path/to/Data Processed/data_app.rds")
  ```

#### In `app_alus.R` *(only needed if you're running this app)*

Same idea as `app.R`, but two files, both loaded near the top of the script:

- `data <- readRDS("data_app.rds")` → point this at the `data_app.rds` produced by `Data_Preprocessing.R`.
- `data_alus <- readRDS("data_app_alus.rds")` → point this at the `data_app_alus.rds` produced by `Alu_Table_Preprocessing_alus.R`.

Either copy both `.rds` files into the same folder as `app_alus.R`, or change both lines to absolute paths, e.g.:
```r
data <- readRDS("/path/to/Data Processed/data_app.rds")
data_alus <- readRDS("/path/to/data_app_alus.rds")
```

---

### Quick Summary of What You Need

**Already in this repo:**
- `Alus_and_CpGs.R`: CpG/Alu counting helpers (sourced automatically).
- `Gene_Annotation.R`: gene annotation helpers (sourced automatically by both preprocessing scripts).
- `Promoter_Annotation.R`: promoter annotation helpers (sourced automatically by both preprocessing scripts).
- `Data_Preprocessing.R`: the bin-level pipeline driver, produces `data_app.rds`.
- `Alu_Table_Preprocessing_alus.R`: the Alu-level pipeline driver, produces `data_app_alus.rds` (needs `data_app.rds` to already exist).
- `app.R`: the bin-level Shiny app.
- `app_alus.R`: the bin- and Alu-level Shiny app (needs both `.rds` files).

**You need to supply yourself:**
- `Metadata_all_runs_combined.csv`
- `counts_bins_norm_mean/` (per-sample JSON bin files)
- Raw per-CpG methylation files/archives (for CpG/Alu counting)
- `Cosmic_CancerGeneCensus_v101_GRCh37.tsv` (download from COSMIC)
- `Promoter_reference_GRCh37.tsv` (build/export from UCSC or Ensembl, see above)
- `CRC_curated_genes.txt` (optional, not currently used)
- `alu_methylation_matrix_all_runs_comb_norm.tsv` (only if you want `app_alus.R`; wide-format per-Alu methylation matrix, see above)

**Generated automatically once you run the pipeline:**
- `data_app.rds`: the file `app.R` (and, in turn, `app_alus.R`) reads. Produced by `Data_Preprocessing.R`.
- `data_app_alus.rds`: the extra file `app_alus.R` reads. Produced by `Alu_Table_Preprocessing_alus.R`.
- Several intermediate CSVs along the way: `Metadata_clean.csv`, `Methylation_long.csv`, `CpG_Alu_bin_annotation.csv`, `Gene_bin_overlaps.csv`, `Promoter_bin_overlaps.csv`, `Bin_annotation_template.csv`, `Bin_prevalence_detection.csv`, `Bin_table.csv`.

**R packages you'll need:**
- Bin-level pipeline (`Data_Preprocessing.R` + the scripts it sources): `jsonlite`, `data.table`.
- Alu-level pipeline (`Alu_Table_Preprocessing_alus.R`): `data.table` (already listed above).
- Either app (`app.R` and `app_alus.R` use the same set): `shiny`, `bslib`, `bsicons`, `DT`, `plotly`, `ggplot2`, `patchwork`.

---

### Steps to Run Everything

1. Install the R packages listed above.
2. Put `Alus_and_CpGs.R`, `Gene_Annotation.R`, `Promoter_Annotation.R`, `Data_Preprocessing.R`, and `Alu_Table_Preprocessing_alus.R` in the same folder.
3. Edit the paths in `Data_Preprocessing.R` (line 7, and lines 15–18) so they point to your script folder and dataset folders.
4. Run the bin-level pipeline:
   ```r
   source("Data_Preprocessing.R")
   ```
   This creates `data_app.rds` inside a `Data Processed/` folder next to your metadata folder.
5. *(Only if you want `app_alus.R`)* Edit the paths in `Alu_Table_Preprocessing_alus.R` (lines 6–10) so `data_app_path` points at the `data_app.rds` from step 4, `alu_matrix_path` points at your `alu_methylation_matrix_all_runs_comb_norm.tsv`, and (optionally) `cosmic_tsv_path`/`promoter_tsv_path` point at the same files used in step 3. Then run:
   ```r
   source("Alu_Table_Preprocessing_alus.R")
   ```
   This creates `data_app_alus.rds`.
6. Copy (or symlink) `data_app.rds` (and, if you ran step 5, `data_app_alus.rds`) into the same folder as whichever app you're launching, or edit the `readRDS()` line(s) near the top of that app's script to point directly to them.
7. Launch the app you want:
   ```r
   shiny::runApp("app.R")       # bin-level only
   # or
   shiny::runApp("app_alus.R")  # bin-level + Alu Table + Alu-aware Genome Browser
   ```
