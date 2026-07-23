# RScript Explanation
##### Goals:

Structure/Calculate values that will be represented in the table (e.g. average methylation, tumor-normal difference, associated genes, functional annotations, prevalence data) prior to starting the application generation.

## Overview

This R script preprocesses DNA methylation data to generate a set of precomputed tables that will later be used by a Shiny application for interactive visualization and analysis.

The pipeline integrates:

- Sample metadata
- Methylation values stored as genomic bins
- Tumor/Normal paired comparisons
- Bin prevalence statistics
- Bin annotation templates

The final output is a single summarized table (`Bin_table.csv`) together with an `.rds` object containing all precomputed data required by the Shiny application.

---

## Workflow

The script performs the following steps:

### 1. Metadata preprocessing

- Reads the metadata file.
- Verifies that the correct delimiter (tab-separated) is used.
- Replaces missing mutation information (`NA`) with `"unknown"` for:
  - KRAS
  - BRAF
  - TP53
- Normalizes sample identifiers.
- Extracts patient IDs from sample names.
- Saves cleaned metadata.

**Output**

```
Metadata_clean.csv
```

---

### 2. Loading methylation bins

Each sample is stored as a JSON file containing methylation values organized as:

```
Chromosome
    ├── Bin position
            └── Methylation value
```

The script:

- Reads every JSON file.
- Converts nested structures into a long-format data frame.
- Converts JSON `null` values into `NA`.

**Output**

```
Methylation_long.csv
```

---

### 3. Missing value classification

For every genomic bin, the script determines whether missing values correspond to:

- **complete**
  - No missing values.

- **sample_specific_missing**
  - Present in some samples but missing in others.

- **structural_gap**
  - Missing in every sample (likely biological absence).

---

### 4. Sample validation

Checks consistency between:

- Metadata sample IDs
- Methylation files

Warnings are generated whenever:

- A metadata sample has no methylation file.
- A methylation file has no metadata entry.

---

### 5. Bin annotation template

Creates a genomic annotation table containing:

- Bin ID
- Chromosome
- Genomic coordinates
- Placeholder columns for:
  - CpG count
  - Alu count
  - Associated genes
  - Functional annotation

These annotation fields are intentionally left empty and will be completed later.

**Output**

```
Bin_annotation_template.csv
```

---

### 6. Tumor–Normal pairing

Samples are matched using their patient ID.

For each patient, the script verifies that exactly one:

- Tumor sample
- Normal sample

are available.

Then, for every genomic bin:

```
Difference = Tumor methylation − Normal methylation
```

Positive values indicate increased methylation in tumor tissue.

---

### 7. Bin prevalence

The script computes, for every genomic bin:

For Tumor samples:

```
prevalence_tumor =
(number of non-missing tumor samples)
/ (total tumor samples)
```

For Normal samples:

```
prevalence_normal =
(number of non-missing normal samples)
/ (total normal samples)
```

Overall prevalence:

```
prevalence_total =
(all non-missing samples)
/ (all samples)
```

**Output**

```
Bin_prevalence_detection.csv
```

---

### 8. Final summary table

For every genomic bin, the script calculates:

- Mean tumor methylation
- Mean normal methylation
- Mean paired tumor-normal difference
- Missing-value classification
- Detection prevalence
- Annotation placeholders

All information is merged into one final table.

**Output**

```
Bin_table.csv
```

---

### 9. Shiny precomputed object

The following objects are stored together in an `.rds` file:

- metadata
- methylation_long
- bin_table

This avoids repeating preprocessing every time the Shiny application starts.

**Output**

```
app_precomputed.rds
```

---

# Input files

```
Dataset/
├── Metadata/
│   └── Metadata_all_runs_combined.csv
│
└── Bins/
    └── counts_bins_norm_mean/
        ├── counts_run01.txt
        ├── counts_run02.txt
        └── ...
```

Each methylation file is expected to be a JSON object with the structure:

```json
{
    "chr1": {
        "1000000": 0.45,
        "2000000": 0.39
    },
    "chr2": {
        "1000000": 0.61
    }
}
```

---

# Output files

```
Processed_1/
├── Metadata_clean.csv
├── Methylation_long.csv
├── Bin_annotation_template.csv
├── Bin_prevalence_detection.csv
├── Bin_table.csv
└── app_precomputed.rds
```

---

# R packages

Required package:

```r
library(jsonlite)
```

---

# Notes

Current placeholders that still require implementation include:

- Associated genes (`genes`)
- Functional annotations (`functional_annotation`)
- CpG counts (`n_cpg`)
- Alu repeat counts (`n_alu`)

These columns are generated as templates and are intended to be completed in future versions of the pipeline.

---
