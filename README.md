# PMOS eNK Implantation Analysis
## Implantation failure linked to altered uterine NK cells in PMOS-like mice

This repository contains the analysis code and visualization scripts for human single-nuclei RNA-seq data examining uterine natural killer (uNK) cells, also referred to as endometrial NK (eNK) cells and TGF-beta signalling.

### Study Overview

**Research Question:** How do endometrial NK express TGF-beta signalling gene markers?

**Key Findings:**
- Three distinct eNK cell populations identified (eNK 1, eNK 2, eNK 3)
- TGF-β signaling pathway expression across eNK subtypes

### Summary of Methodology

- **Technology:** Single-cell RNA sequencing (scRNA-seq)
- **Cell Type:** Endometrial Natural Killer (eNK) cells
- **Samples:** Control vs PCOS groups
- **Integration Method:** Canonical Correlation Analysis (CCA) with Seurat v5
- **Analysis Framework:** R (Seurat, ggplot2, patchwork)

### Prerequisites

**Required R version:** R ≥ 4.0.0

**Required R packages:**
```r
# Core packages
install.packages(c("Seurat", "SeuratObject", "SeuratDisk"))

# Data manipulation
install.packages(c("dplyr", "openxlsx"))

# Visualization
install.packages(c("ggplot2", "patchwork", "cowplot", "RColorBrewer", "viridis"))

# Analysis
install.packages(c("cluster", "harmony"))
```

### Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/PCOS_eNK_Implantation_Analysis.git
cd PCOS_eNK_Implantation_Analysis

# Create necessary directories
mkdir -p data/Core_seurat_objects
mkdir -p data/Subset_seurat_objects
mkdir -p output/0_Preprocessing
mkdir -p output/Publication_Figures
```

### Data Setup

1. Place your original Seurat object(s) in `data/Core_seurat_objects/`
2. Expected file: `Endo_All_Immune_reclustered_labelled.h5seurat`

---

## 📋 Analysis Workflow

### Step 1: Data Preprocessing and Integration
```r
source("scripts/01_Preprocessing_With_Integration.R")
```

**This script performs:**
- Loading of parent immune cell atlas
- Subsetting for Control and PCOS groups
- Extraction of uNK/eNK cell populations (uNK 1, 2, 3)
- Quality control and filtering (min 100 cells/sample)
- Sample integration using CCA (or Harmony/RPCA)
- Dimensionality reduction (PCA, UMAP)
- Clustering and cell type annotation
- Integration quality assessment

**Key Parameters:**
- `integration.method`: "CCA" (default), "RPCA", or "Harmony"
- `n.features`: 2000 variable features
- `n.dims`: 30 principal components
- `resolution`: 0.5 for clustering
- `min_cells_per_sample`: 100

**Outputs:**
- `output/0_Preprocessing/uNK_QC_Summary_*.xlsx` - QC metrics
- `output/0_Preprocessing/Integration_Metrics.csv` - Integration quality
- `output/0_Preprocessing/Integration_Comparison_*.pdf` - Integration visualizations
- `data/Subset_seurat_objects/Endometrium_uNK_integrated_CCA.rds` - Integrated object

### Step 3: Generate Publication Figures
```r
source("scripts/02_Plotting.R")
```

**This script generates:**

1. **UMAP Visualizations**
   - Cell type composition
   - Sample distribution
   - Group comparisons (Control vs PCOS)
   - Cluster identification

2. **Marker Gene Analysis**
   - eNK validation markers (NCAM1, FGFBP2, ENTPD1, ITGB2, CD160, ITGAE)
   - TGF-β pathway genes (TGFB1, TGFBR1, TGFBR2, SMAD2-4, SMAD7)

3. **Expression Visualizations**
   - Feature plots (UMAP with gene expression)
   - Violin plots (distribution by cell type and group)
   - Dot plots (percent expressed + average expression)
   - Co-expression analysis (TGF-β ligand-receptor pairs)

4. **Summary Statistics**
   - Cell counts per type and group
   - Percentage distributions

**All figures saved to:** `output/Publication_Figures/`

---

## 🎨 Key Visualizations

### Generated Figures

| Figure | Description | File |
|--------|-------------|------|
| UMAP - Cell Types | eNK 1, 2, 3 subtypes | `UMAP_CellType.pdf` |
| UMAP - Groups | Control vs PCOS distribution | `UMAP_Group.pdf` |
| UMAP - Split | Cell types split by group | `UMAP_Split_Immune.pdf` |
| Feature Plots - eNK | Validation markers | `FeaturePlot_eNK_markers.pdf` |
| Feature Plots - TGF-β | Pathway components | `FeaturePlot_TGFb_markers.pdf` |
| Dot Plots | Expression heatmaps | `Dotplot_TGFb_markers_split.pdf` |
| Violin Plots | Expression distributions | `VlnPlot_TGFb_markers_split.pdf` |
| Co-expression | Ligand-receptor pairs | `UMAP_CoExp_TGFB1_TGFBR*.pdf` |

---

## 🧬 Marker Genes

### eNK Cell Validation Markers

| Gene | Protein | Cell Type | Function |
|------|---------|-----------|----------|
| NCAM1 | CD56 | All eNK | Pan-NK marker |
| FGFBP2 | - | eNK (CD16-) | NK secreted protein |
| ENTPD1 | CD39 | eNK 1 | Ectonucleotidase |
| ITGB2 | CD18 | eNK 2 | Integrin beta-2 |
| CD160 | CD160 | eNK 2 | NK receptor |
| ITGAE | CD103 | eNK 3 | Integrin alpha-E |

### TGF-β Signaling Pathway

| Gene | Function | Expression Pattern |
|------|----------|-------------------|
| TGFB1 | Ligand | High in eNK 1 |
| TGFBR1 | Type I receptor | Broadly expressed |
| TGFBR2 | Type II receptor | Subset-specific |
| SMAD2/3 | Signal transducers | Moderate |
| SMAD4 | Co-SMAD | Moderate |
| SMAD7 | Inhibitory SMAD | Lower expression |

---

## 📈 Results Summary

### Cell Composition
- **Total eNK cells analyzed:** ~2,000-3,000 (filtered dataset)
- **Groups:** Control (n samples), PCOS_W0 (n samples)
- **Cell populations:** 3 distinct eNK subtypes

### Key Findings

1. **eNK Heterogeneity**
   - Three transcriptionally distinct eNK populations
   - Consistent clustering across samples

2. **TGF-β Pathway Expression**
   - TGFB1 and TGFBR1 highly expressed in eNK 1
   - TGFBR2 shows subset-specific expression
   - Signaling components (SMADs) broadly present

3. **Group Differences**
   - Altered eNK composition in PCOS
   - Potential implications for implantation

4. **Integration Quality**
   - Successful batch correction across samples
   - Preservation of biological signal
   - Good sample mixing with maintained cell type separation

---

## 🔧 Customization

### Changing Integration Method

In `01_Preprocessing_With_Integration.R`, modify:

```r
# Options: "CCA", "RPCA", "Harmony", "v3"
integration.method = "CCA"
```

### Adjusting Clustering Resolution

```r
# Higher values = more clusters
resolution = 0.5  # Default
resolution = 0.3  # Fewer clusters
resolution = 0.8  # More clusters
```

### Modifying Marker Genes

In `02_Plotting.R`:

```r
# Add or remove genes from analysis
TGFb_markers <- c(
  "TGFB1", "TGFBR1", "TGFBR2",
  "YOUR_GENE_HERE"  # Add custom genes
)
```

### Color Palettes

```r
# Modify in 02_Plotting.R
CT.colors = brewer.pal(3, "Set2")      # Cell type colors
group.colors = c("#A0A0A0", "#B0D8D6") # Group colors
```

---

## 📊 Quality Control Metrics

### Outputs from Preprocessing

1. **QC Summary Tables** (`uNK_QC_Summary_*.xlsx`)
   - nCount_RNA: Total UMI per cell
   - nFeature_RNA: Unique genes per cell
   - percent.mt: Mitochondrial percentage
   - percent.ribo: Ribosomal percentage
   - percent.hb: Hemoglobin percentage
   - S.Score / G2M.Score: Cell cycle scores

2. **Integration Metrics** (`Integration_Metrics.csv`)
   - Sample mixing score (silhouette)
   - Cell type separation score (silhouette)
   - Comparison: Integrated vs Unintegrated

3. **Visualization Outputs**
   - Violin plots for QC metrics by sample/group
   - PCA plots showing batch effects
   - Integration comparison UMAPs

---

## 🤝 Contributing

If you use this pipeline and make improvements:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit your changes (`git commit -am 'Add improvement'`)
4. Push to the branch (`git push origin feature/improvement`)
5. Open a Pull Request

---

## 📝 Citation

If you use this code, please cite:

```
[Your Name] et al. (2025). Implantation Failure Linked to Altered Uterine NK 
Cells in PCOS-like Mice. [Journal Name], [Volume]([Issue]), [Pages].
```

### Related Publications

*List any related papers or preprints here*

---

## 📧 Contact

**Principal Investigator:** [PI Name]  
**Bioinformatics Analysis:** Gustaw Eriksson  
**Institution:** [Your Institution]  

For questions about the analysis:
- Open an issue on GitHub
- Contact: [your.email@institution.edu]

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- Seurat team for excellent single-cell analysis tools
- [Funding sources]
- [Collaborators]
- [Core facilities]

---

## 📚 References

### Key Methods Papers

1. Hao et al. (2021). Integrated analysis of multimodal single-cell data. *Cell*. [Seurat v4]
2. Stuart et al. (2019). Comprehensive Integration of Single-Cell Data. *Cell*. [Seurat v3]
3. Korsunsky et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*.

### Biological Context

*Add relevant references for uNK cells, PCOS, and implantation biology*

---

## 📋 Version History

- **v1.0.0** (2025-01) - Initial release
  - Complete preprocessing pipeline
  - Integration workflow with multiple methods
  - Publication-ready visualization scripts

---

## ⚠️ Notes

### Data Privacy
- Original sample identifiers have been anonymized
- Raw sequencing data is not included in this repository
- Processed objects available upon reasonable request


### Known Issues
- Very small samples (<100 cells) are filtered out
- Integration parameters are optimized for this dataset size
- Adjust k.weight, k.anchor parameters for different sample sizes

---

**Last Updated:** January 2025  
**Status:** Active Development
