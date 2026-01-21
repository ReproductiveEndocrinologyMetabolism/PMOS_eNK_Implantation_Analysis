# R Package Requirements for PCOS-uNK-implantation-analysis

# This file lists all required R packages and their recommended versions
# Install all packages by running: source("environment/requirements.R")

# ==============================================================================
# INSTALLATION INSTRUCTIONS
# ==============================================================================

# 1. Install R (version >= 4.0.0)
#    Download from: https://cran.r-project.org/

# 2. Install RStudio (optional but recommended)
#    Download from: https://posit.co/downloads/

# 3. Run this script to install all dependencies:
#    source("environment/requirements.R")

# ==============================================================================
# CORE DEPENDENCIES
# ==============================================================================

# Check R version
if (getRversion() < "4.0.0") {
  stop("R version 4.0.0 or higher is required. Please update R.")
}

# Function to install if not present
install_if_missing <- function(pkg, version = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("Installing", pkg, "...\n"))
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(paste(pkg, "is already installed.\n"))
  }
}

# ==============================================================================
# BIOCONDUCTOR SETUP
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# ==============================================================================
# SEURAT ECOSYSTEM (Single-Cell Analysis)
# ==============================================================================

cat("\n=== Installing Seurat Ecosystem ===\n")

# Seurat v5 (main package)
install_if_missing("Seurat")  # Recommended: >= 5.0.0

# Seurat utilities
install_if_missing("SeuratObject")  # Seurat object structure
install_if_missing("SeuratDisk")    # Read/write h5seurat files

# Note: If you need Seurat v4 compatibility:
# install.packages("Seurat", version = "4.4.0")

# ==============================================================================
# DATA MANIPULATION
# ==============================================================================

cat("\n=== Installing Data Manipulation Packages ===\n")

install_if_missing("dplyr")      # Data manipulation (tidyverse)
install_if_missing("tidyr")      # Data tidying (optional but useful)
install_if_missing("openxlsx")   # Excel file I/O

# ==============================================================================
# VISUALIZATION
# ==============================================================================

cat("\n=== Installing Visualization Packages ===\n")

install_if_missing("ggplot2")      # Core plotting
install_if_missing("patchwork")    # Combine plots
install_if_missing("cowplot")      # Publication-ready themes
install_if_missing("RColorBrewer") # Color palettes
install_if_missing("viridis")      # Perceptually uniform colors
install_if_missing("scales")       # Scale functions for ggplot2

# ==============================================================================
# STATISTICAL ANALYSIS
# ==============================================================================

cat("\n=== Installing Statistical Analysis Packages ===\n")

install_if_missing("cluster")  # Clustering algorithms and silhouette analysis
install_if_missing("Matrix")   # Sparse matrix operations (usually auto-installed)

# ==============================================================================
# INTEGRATION METHODS
# ==============================================================================

cat("\n=== Installing Integration Methods ===\n")

# Harmony integration (if using Harmony method)
install_if_missing("harmony")

# ==============================================================================
# VERIFY INSTALLATION
# ==============================================================================

cat("\n=== Verifying Installation ===\n")

required_packages <- c(
  "Seurat", "SeuratObject", "SeuratDisk",
  "dplyr", "openxlsx",
  "ggplot2", "patchwork", "cowplot", "RColorBrewer", "viridis",
  "cluster", "harmony"
)

missing <- c()
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing <- c(missing, pkg)
  }
}

if (length(missing) > 0) {
  cat("\n⚠️  WARNING: The following packages failed to install:\n")
  cat(paste(" -", missing, collapse = "\n"))
  cat("\n\nPlease install these manually.\n")
} else {
  cat("\n✅ All required packages installed successfully!\n")
}