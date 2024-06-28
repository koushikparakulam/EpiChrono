## Overview

**EpigChrono** R workflow designed to process, annotate, and analyze epigenetic data, particularly focusing on histone modification (e.g., H3K36me3) in the context of aging research. It integrates various steps, from raw BED file preprocessing to functional enrichment and multi-omics integratio to provide a complete view of epigenetic regulation.

## Repository Structure

```
.
├─ data/
│   └─ GSM466737_UCSD.H1.H3K36me3.LL160.bed      # Example BED file
├─ R/
│   ├─ Utilities_Setup.R                               # Utility and support functions
│   ├─ Preprocessing.R                           # BED reading and cleaning
│   ├─ Annotations.R                              # Annotation of genomic regions
│   ├─ Epigen_Stats.R                              # Statistical analyses and modeling
│   ├─ Functional_Enrichment.R                   # GO/KEGG, GSEA, hypergeometric tests
│   ├─ Omics.R                       # Expression & methylation integration
│   ├─ Visualize_Omics_Enrichment.R                           # Plotting fnctions
│   └─ Pipeline_Workflow.R                           # Main entry point
├─ analysis_results/                              # Output directory (created automatically)
```

## Getting Started

### Prerequisites

1. **Install R** (versions can work differently recc. 4.0).
2. **Clone this repo**.
3. **Place your BED file** in `./data/`.
4. **Install dependencies**:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="https://cran.r-project.org")

# Example of installing key dependencies
BiocManager::install(c("rtracklayer", "GenomicRanges", "GenomicFeatures",
                        "org.Hs.eg.db", "clusterProfiler", "BSgenome.Hsapiens.UCSC.hg19"))
install.packages(c("tidyverse", "reshape2", "lme4"))
```

5. **Run the pipeline** in an R session:

```r
source("R/main_pipeline.R")
main_pipeline(
  bed_path="data/GSM466737_UCSD.H1.H3K36me3.LL160.bed",  # path to your BED file
  blacklist_path=NULL,                                  # optional
  genome_version="hg19",
  expr_file=NULL,                                       # path to expression data .csv
  meth_file=NULL,                                       # path to methylation data .csv
  output_dir="analysis_results"
)
```

## Workflow Steps

1. **Initialization**
   - Sets random seed, configures environment.
   - Checks and installs packages if missing.

2. **Preprocessing**
   - Reads BED data, filters short/long regions, drops unusual chromosomes, and optionally removes blacklisted regions.
   - Rescales any score column for uniform 0–1 distribution (optional).

3. **Annotation**
   - Builds/loads TxDb, identifies nearest genes, calculates distance.
   - Marks promoter overlaps.
   - Measures GC content if reference genome is installed.

4. **Statistics**
   - Summarizes region counts, average lengths, etc.
   - Fits linear models and group comparisons if a grouping column is present or simulated.

5. **Functional Enrichment**
   - Uses **clusterProfiler** for GO (BP, MF, CC) and KEGG pathways.
   - Supports GSEA with a user-supplied or randomly simulated ranking.

6. **Integration with Omics**
   - Merges expression and methylation data by gene IDs.
   - Optionally runs **edgeR** / **DESeq2** for differential expression.
   - Correlates epigenetic scores with expression/methylation.
   - K-means clustering to identify multi-omics pattern clusters.

7. **Visualization**
   - Basic histograms, density plots, bar charts for enrichment, correlation matrices, and cluster plots.
   - Saves output as PNG images in `analysis_results/`.

8. **Outputs**
   - Final annotated BED data (`final_bed.csv`).
   - Summary stats (`summary_stats.csv`).
   - Differential test results (`diff_test.txt`).
   - GO/KEGG outputs, correlation matrix, GSEA results, etc.
   - Diagnostic figures and session info.

## Customization

- **Adjust thresholds** for region filtering in `preprocessing.R` (e.g., `apply_qc_filters`).
- **Toggle annotation features** (promoter ranges, GC content, etc.) in `annotation.R`.
- **Change the formulas** for linear/mixed models in `statistics.R`.
- **Customize** your background gene sets or ranking in `functional_enrichment.R`.
- **Use real data** for expression/methylation in `integration_omics.R`.
- **Add or remove plots** in `visualization.R`.

## Known Limitations

- Default references assume **hg19** for annotation and certain Bioconductor packages. Adjust if using hg38 or other species.
- Functional enrichment depends heavily on **clusterProfiler**. You may need alternate approaches if analyzing non-human data or if your gene IDs are not recognized.
- Merged multi-omics data may require consistent gene identifiers (e.g., ENTREZ, ENSEMBL, SYMBOL).
- There are still issues regarding the API calls made, and work is being done to fix them
---

### Example Output Files

- **Annotated BED File:** `analysis_results/final_bed.csv`
- **Summary Statistics:** `analysis_results/summary_stats.csv`
- **GO Enrichment Results:** `analysis_results/go_results.txt`
- **Differential Expression Results:** `analysis_results/diff_test.txt`
- **Methylation-Expression Correlation Matrix:** `analysis_results/omics_correlation.txt`
- **Visualization Files:**
  - `analysis_results/length_hist.png`
  - `analysis_results/score_density.png`
  - `analysis_results/enrichment_plot.png`
  - `analysis_results/correlation_plot.png`

---

