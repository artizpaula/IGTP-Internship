# Exploration of Epigenomic Data in Colorectal Cancer

##### Paula Artiz Dueñas, UPC

### What This Project Is About

This project is an interactive **R Shiny web app** for exploring and visualizing DNA methylation data in colorectal cancer. It lets users look at genomic "bins" (1 Mb chunks of the genome), compare methylation between tumor and normal samples, filter by region, and explore gene annotations, all through interactive plots and tables.

### What the App Can Do

The app has a sidebar for picking data, plus several tabs for exploring it.

**Sidebar, choosing your data:**
- Pick a chromosome and search/select specific bins (e.g. `1_1000000`).
- Add every bin from a chromosome at once, or clear your selection.
- Upload your own list of bins from a `.csv`, `.tsv`, or `.txt` file.

**Tabs, exploring and visualizing:**
1. **Overview**: summary panels comparing tumor vs. normal methylation, mutation status, clinical stage/MSI-MSS, and sex distribution for the bins you selected.
2. **Genome Browser**: scroll through methylation across the whole genome; jump to a chromosome, zoom, or search by coordinates, bin ID, or gene name.
3. **Tumor vs. Normal**: density plots, PCA/UMAP, and a patient similarity network based on methylation.
4. **Genome-wide Profile**: a Manhattan-style plot of tumor–normal methylation differences, flagging significant/outlier bins.
5. **Feature × Chromosome Heatmap**:methylation shift per patient and chromosome, with an optional color strip for a clinical feature.
6. **Clinical Explorer**: compare tumor methylation by mutation status (KRAS, BRAF, TP53) and browse/filter the clinical metadata table.
7. **Bin Table**: a searchable, filterable table of all bins with their stats and annotations, downloadable as a CSV.

---

### How the Pipeline Works (Run These in Order!)

The app never touches the raw data directly, it just loads one pre-built file, `data_app.rds`. That file is created by running three scripts, in this order:

| Step | Script | What it does |
|---|---|---|
| 1 | `Alus_and_CpGs.R` | Helper functions that read the raw per-CpG methylation tables and count CpGs/Alu elements per 1 Mb bin. You don't run this yourself, it gets `source()`-d automatically by `Week_2_With_Prevalence.R`. |
| 2 | `Gene_Annotation.R` | Helper functions that match bin coordinates against the COSMIC Cancer Gene Census to add gene names/IDs. Also `source()`-d automatically, not run directly. |
| 3 | `Week_2_With_Prevalence.R` | **This is the one you actually run.** It loads the metadata, reshapes the per-sample JSON files, sources the two scripts above, computes CpG/Alu counts and gene annotation, works out prevalence stats, builds the final bin table, and saves `data_app.rds`. |
| 4 | `app.R` | The Shiny app itself. It just loads `data_app.rds` (from step 3), launch it with `shiny::runApp()`, it's not part of the "pipeline" as such. |

**You only run `Week_2_With_Prevalence.R`. As long as all three `.R` scripts are sitting in the same folder, it takes care of the rest.**

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

#### Suggested folder layout

Any layout works as long as you update the path variables to match (see below), but here's a layout that mirrors the original project:

```
Project/
├── R Scripts/
│   ├── Alus_and_CpGs.R
│   ├── Gene_Annotation.R
│   ├── Week_2_With_Prevalence.R
│   └── app.R
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

Running the pipeline will also create a `Data Processed/` folder (next to `Dataset/Metadata/`, i.e. `Project/Dataset/Data Processed/`) with all the intermediate CSVs plus the final `data_app.rds`.

---

### Paths You Need to Edit Before Running

The scripts as provided still have **hardcoded paths from the original author's laptop**. These won't work anywhere else, so you'll need to change them.

#### In `Week_2_With_Prevalence.R`

- **Line 7**: `setwd("/Users/paulaartizduenas/Desktop/Project/R Scripts")`
  → Change this to wherever `Alus_and_CpGs.R`, `Gene_Annotation.R`, and `Week_2_With_Prevalence.R` live on your machine (they all need to be in the same folder, since lines 9 and 11 `source()` them by relative path).

- **Lines 15–18**: the four dataset root paths:
  ```r
  dataset_metadata <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Metadata"
  dataset_bins <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Bins"
  dataset_alus_cpgs <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Alus_CpGs"

  dataset_gene_annotation <- "/Users/paulaartizduenas/Desktop/Project/Dataset/Archivo"
  ```
  → Update **all four** to point to your local copies of the data described above. Everything else (`metadata_dir`, `bins_dir`, `output_dir`, `cosmic_tsv_path`, `crc_gene_list_path`) gets built automatically from these with `file.path()` so nothing past line 28 needs editing, as long as the folder/file names inside each of the four roots match exactly what's listed under "Data You Need to Provide" (e.g. the metadata file needs to be literally named `Metadata_all_runs_combined.csv`).

#### In `Alus_and_CpGs.R` and `Gene_Annotation.R`

- No hardcoded paths here, every function takes a `filepath`/`cosmic_tsv_path` argument as input. Just keep them in the same folder as `Week_2_With_Prevalence.R` so the `source()` calls above can find them.

#### In `app.R`

- **Line 12**: `data <- readRDS("data_app.rds")`
  → Since this is a relative path, either copy `data_app.rds` (from the pipeline's `Data Processed/` output) into the same folder as `app.R`, or change this line to point to it directly, e.g.:
  ```r
  data <- readRDS("/path/to/Data Processed/data_app.rds")
  ```
---

### Quick Summary of What You Need

**Already in this repo:**
- `Alus_and_CpGs.R`: CpG/Alu counting helpers (sourced automatically).
- `Gene_Annotation.R`: gene annotation helpers (sourced automatically).
- `Week_2_With_Prevalence.R`: the pipeline driver, produces `data_app.rds`.
- `app.R`: the Shiny app.

**You need to supply yourself:**
- `Metadata_all_runs_combined.csv`
- `counts_bins_norm_mean/` (per-sample JSON bin files)
- Raw per-CpG methylation files/archives (for CpG/Alu counting)
- `Cosmic_CancerGeneCensus_v101_GRCh37.tsv` (download from COSMIC)
- `CRC_curated_genes.txt` (optional, not currently used)

**Generated automatically once you run the pipeline:**
- `data_app.rds`: the file the app actually reads.
- Several intermediate CSVs along the way: `Metadata_clean.csv`, `Methylation_long.csv`, `CpG_Alu_bin_annotation.csv`, `Gene_bin_overlaps.csv`, `Bin_annotation_template.csv`, `Bin_prevalence_detection.csv`, `Bin_table.csv`.

**R packages you'll need:**
- Pipeline (`Week_2_With_Prevalence.R` + the scripts it sources): `jsonlite`, `data.table`.
- App (`app.R`): `shiny`, `bslib`, `bsicons`, `DT`, `plotly`, `ggplot2`, `patchwork`.

---

### Steps to Run Everything

1. Install the R packages listed above.
2. Put `Alus_and_CpGs.R`, `Gene_Annotation.R`, and `Week_2_With_Prevalence.R` in the same folder.
3. Edit the paths in `Week_2_With_Prevalence.R` (line 7, and lines 15–18) so they point to your script folder and dataset folders.
4. Run the pipeline:
   ```r
   source("Week_2_With_Prevalence.R")
   ```
   This creates `data_app.rds` inside a `Data Processed/` folder next to your metadata folder.
5. Copy (or symlink) `data_app.rds` into the same folder as `app.R` or edit line 12 of `app.R` to point directly to it.
6. Launch the app:
   ```r
   shiny::runApp("app.R")
   ```