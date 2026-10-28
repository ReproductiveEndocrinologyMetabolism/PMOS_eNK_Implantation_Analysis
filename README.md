README
Author: Gustaw Eriksson
Date: 2026-01-21
Contact: gustaw.eriksson@ki.se

# Single-nuclei eNK-TGFβ Analysis

Single-nuclei RNA-seq pipeline for endometrial NK (eNK) cells in the PMOS endometrium project: subsetting and integration of eNK cells from the endometrial immune atlas (Control vs PMOS), followed by publication figures focused on eNK subset markers and TGF-β signalling.

## Files

### `01_Preprocessing_With_Integration.R`

Subsets and integrates eNK cells from the PMOS endometrium cell atlas.In short: 

- Loading the parent immune atlas (`.h5seurat`) and subsetting to the Control and PMOS (`PMOS_W0`, baseline) groups
- Converting the RNA assay to Seurat v5 (`Assay5`) and harmonising metadata columns (`Sample.id`, `Main_Cluster`, `Group`, `Subgroup`)
- Subsetting to the three uNK/eNK populations and excluding samples with fewer than `min_cells_per_sample` cells (default 100)
- Per-sample and per-group QC summaries (nCount, nFeature, percent mito/haemoglobin/ribosomal, cell-cycle scores) exported to Excel
- An unintegrated baseline for comparison
- Integration across samples with one of: Seurat v5 `CCAIntegration` (default), `RPCAIntegration`, `HarmonyIntegration`, or the original Seurat v3 anchor workflow (`"v3"`). Anchor/weighting parameters (`k.anchor`, `k.filter`, `k.score`, `k.weight`) are scaled automatically to the smallest sample so integration does not fail on small samples
- Clustering and UMAP on the integrated reduction, plus a side-by-side integrated vs unintegrated UMAP comparison (by sample and by cell type)
- Joining layers and saving the integrated object as both `.h5seurat` and `.rds`

### `02_Plotting.R`

Generates publication figures from the integrated object produced by script 01. Renames `uNK` → `eNK` and `PMOS_W0` → `PMOS` for display, then produces:

- UMAPs coloured by eNK subset, group, and Seurat cluster, plus UMAPs split by group
- Co-expression (blend) UMAPs for TGFB1/TGFBR1, TGFB1/TGFBR2 and TGFBR1/TGFBR2
- Feature plots, violin plots (with and without split by group) and dot plots for:
  - eNK validation markers: NCAM1 (CD56), FGFBP2 (CD16-negative), ENTPD1 (CD39, eNK 1), ITGB2 (CD18, eNK 2), CD160 (eNK 2), ITGAE (CD103, eNK 3)
  - TGF-β pathway genes: TGFB1, TGFBR1, TGFBR2, SMAD2, SMAD3, SMAD4, SMAD7 (and a short ligand/receptor set: TGFB1, TGFBR1, TGFBR2)
- A summary table of cell counts and proportions per eNK subset, broken down by group

Every figure is saved as both PDF and PNG (300 dpi).

## Usage

1. Place the PMOS endometrium immune atlas at `Data/Seurat_Immune/Endo_All_Immune_reclustered_labelled.h5seurat`, or update `seurat.path`.
2. Set the variables at the top of `01_Preprocessing_With_Integration.R`: groups and cell types to subset, `integration.method` (`"CCA"`, `"RPCA"`, `"Harmony"` or `"v3"`), `n.features`, `n.dims`, `resolution` and `min_cells_per_sample`.
3. Run `01_Preprocessing_With_Integration.R` top to bottom. Outputs:
   - `Output/0_Preprocessing/`: QC summary tables (`uNK_QC_Summary_Sample.xlsx`, `uNK_QC_Summary_Group.xlsx`) and `Integration_Comparison_<method>.pdf`
   - `Data/Subset_seurat_objects/`: `Endometrium_uNK_integrated_<method>.h5seurat` and `.rds`
4. Run `02_Plotting.R`. It reads `Endometrium_uNK_integrated_CCA.rds` by default. Update `seurat.rds.path` if a different integration method was used. Figures and summary tables are written to `Output/Publication_Figures/`.

## Input requirements

The Seurat object is expected to contain the metadata columns `orig.ident`, `Group_Stage`, `Group_Treatment`, `Labelled_Clusters_SCT.1` and `Immune_labelled`, as well as the precomputed QC columns `percent.mt`, `percent.hb`, `percent.ribo`, `S.Score` and `G2M.Score`.

## Requirements

R (≥ 4.3) with Seurat (v5), SeuratObject, SeuratDisk, ggplot2, patchwork, cowplot, dplyr, cluster, RColorBrewer, viridis, openxlsx and reticulate. `harmony` is required only when `integration.method = "Harmony"`.
