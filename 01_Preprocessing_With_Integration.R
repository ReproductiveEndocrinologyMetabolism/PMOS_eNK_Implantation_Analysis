# Author: Gustaw Eriksson
# Date: 250904 (Updated: 2025)
# Description: Script to subset and integrate uNK cells from the immune atlas
# Using Seurat v5 integration methods

# Loading the packages
library(openxlsx)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cluster)
library(cowplot)
library(reticulate)

### Setting variables section ###
# Set the paths to the scRNA/snRNA-seq atlases to process
seurat.path = "Data/Seurat_Immune/Endo_All_Immune_reclustered_labelled.h5seurat"
names(seurat.path) = "Endometrium_Immune"

# Which groups to subset
subset.groups = c("Control", "PCOS_W0")
subset.immune = c("uNK 1", "uNK 2", "uNK 3")

# Set the order of the samples and groups
group.order = c("Control", "PCOS_W0")

# Add color palettes for samples and groups
sample.colors = c(rep("#A0A0A0", 5), rep("#B0D8D6", 12),
                  rep("#95BFE1", 7), rep("#65D46E", 4))
group.colors = c("#A0A0A0", "#B0D8D6")
CT.colors = brewer.pal(3, "Set2")

# Integration parameters
integration.method = "CCA"  # Options: "CCA", "RPCA", "Harmony"
n.features = 2000
n.dims = 30
resolution = 0.5
min_cells_per_sample <- 100

# Generating the data directory for quality control output
if (dir.exists(path = paste0("Output/0_Preprocessing")) == FALSE) {
  print(paste0("Generating output directory Output/0_Preprocessing"))
  dir.create(path = paste0("Output/0_Preprocessing"), recursive = TRUE)
  Output.QC.dir = paste0("Output/0_Preprocessing/")
} else {
  print("Directory exists")
  Output.QC.dir = paste0("Output/0_Preprocessing/")
}

# Generating the data directory for subset atlases
if (dir.exists(path = paste0("Data/Subset_seurat_objects")) == FALSE) {
  print(paste0("Generating output directory Data/Subset_seurat_objects"))
  dir.create(path = paste0("Data/Subset_seurat_objects"), recursive = TRUE)
  Output.seurat.dir = paste0("Data/Subset_seurat_objects/")
} else {
  print("Directory exists")
  Output.seurat.dir = paste0("Data/Subset_seurat_objects/")
}

### Seurat atlas loading and subsetting ###
# Load the Seurat object
print("Loading parent Seurat object...")
seurat.x = SeuratDisk::LoadH5Seurat(seurat.path)

# Check the group names
print("Available groups:")
print(table(seurat.x$Group_Stage))

# Subset the Seurat object based on selected groups
print(paste0("Subsetting for groups: ", paste(subset.groups, collapse = ", ")))
seurat.x = subset(x = seurat.x, subset = Group_Stage %in% subset.groups)
print(paste0("After group subsetting: ", ncol(seurat.x), " cells"))

# Update to Seurat v5
seurat.x = UpdateSeuratObject(seurat.x)
seurat.x[["RNA"]] = as(seurat.x[["RNA"]], Class = "Assay5")
DefaultAssay(seurat.x) = "RNA"

# Rename ident class for consistency
seurat.x$Sample.id = seurat.x$orig.ident
seurat.x$Main_Cluster = seurat.x$Labelled_Clusters_SCT.1
seurat.x$Group = seurat.x$Group_Stage
seurat.x$Subgroup = seurat.x$Group_Treatment

Idents(seurat.x) = seurat.x$Sample.id

# Check if "Immune_labelled" exists
if (!"Immune_labelled" %in% colnames(seurat.x@meta.data)) {
  stop("The column 'Immune_labelled' does not exist in the metadata.")
}

print("Available immune cell types:")
print(table(seurat.x$Immune_labelled))

# Subset the uNK cells
print(paste0("Subsetting for cell types: ", paste(subset.immune, collapse = ", ")))
seurat.NK <- subset(x = seurat.x, subset = Immune_labelled %in% subset.immune)
print(paste0("After cell type subsetting: ", ncol(seurat.NK), " cells"))
print(table(seurat.NK$Immune_labelled))

# Clean up
rm(seurat.x)
gc()

### QC and Sample Filtering ###
# Check cell counts per sample
print("Cell counts per sample:")
cell_counts <- table(seurat.NK$Sample.id)
print(cell_counts)

# Filter samples with too few cells
samples_to_keep <- names(cell_counts[cell_counts >= min_cells_per_sample])
samples_to_exclude <- names(cell_counts[cell_counts < min_cells_per_sample])

if (length(samples_to_exclude) > 0) {
  print(paste0("WARNING: Excluding ", length(samples_to_exclude), " samples with < ", 
               min_cells_per_sample, " cells:"))
  print(samples_to_exclude)
  print(cell_counts[samples_to_exclude])
  
  seurat.NK <- subset(seurat.NK, subset = Sample.id %in% samples_to_keep)
  
  print(paste0("Remaining samples: ", length(samples_to_keep)))
  print(table(seurat.NK$Sample.id))
} else {
  print("All samples have sufficient cells")
}

# QC summary
# Function to generate QC related plots and tables
QC_summary = function(seurat.obj, 
                      QC.features = c("nCount_RNA", "nFeature_RNA", 
                                      "percent.mt", "percent.hb", "percent.ribo",
                                      "S.Score", "G2M.Score")) {
  # Load the meta data 
  seurat.meta = seurat.obj@meta.data
  
  # Summarize per sample ID
  summary.sample = seurat.meta %>%
    group_by(Sample.id) %>%
    summarise(across(all_of(QC.features), 
                     ~ round(mean(as.numeric(.), na.rm = TRUE), 2),
                     .names = "mean_{.col}"),
              nCells = n())
  
  # Summarize per Group
  summary.group <- seurat.meta %>%
    group_by(Group) %>%
    summarise(across(all_of(QC.features), 
                     ~ round(mean(as.numeric(.), na.rm = TRUE), 2),
                     .names = "mean_{.col}"),
              nCells = n())
  
  # Violin plots per sample ID
  Vln.plots.sample <- lapply(QC.features, function(f) {
    VlnPlot(seurat.obj, features = f, group.by = "Sample.id", pt.size = 0) +
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") +
      ggtitle(paste("QC:", f))
  })
  names(Vln.plots.sample) = QC.features
  
  # Violin plots per Group
  Vln.plots.group <- lapply(QC.features, function(f) {
    VlnPlot(seurat.obj, features = f, group.by = "Group", pt.size = 0) +
      theme_cowplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") +
      ggtitle(paste("QC:", f))
  })
  names(Vln.plots.group) = QC.features
  
  test.list = list(
    summary_per_sample = summary.sample,
    summary_per_group = summary.group,
    plots_sample = wrap_plots(Vln.plots.sample) + plot_layout(ncol = 2),
    plots_group = wrap_plots(Vln.plots.group) + plot_layout(ncol = 2)
  )
  
  # Return the output as a list
  return(list(
    summary_sample = summary.sample,
    summary_group  = summary.group,
    plots_sample = Vln.plots.sample,
    plots_group  = Vln.plots.group
  ))
}

seurat.QC = QC_summary(seurat.obj = seurat.NK)

write.xlsx(x = seurat.QC$summary_sample,
           file = paste0(Output.QC.dir, "uNK_QC_Summary_Sample.xlsx"),
           asTable = FALSE, overwrite = TRUE)
write.xlsx(x = seurat.QC$summary_group,
           file = paste0(Output.QC.dir, "uNK_QC_Summary_Group.xlsx"),
           asTable = FALSE, overwrite = TRUE)

### Create Unintegrated Baseline for Comparison ###
# Create a copy for unintegrated analysis
seurat.NK.unint <- seurat.NK

# Join layers if needed
if (inherits(seurat.NK.unint[["RNA"]], "Assay5")) {
  seurat.NK.unint[["RNA"]] <- JoinLayers(seurat.NK.unint[["RNA"]])
}

# Standard workflow without integration
seurat.NK.unint <- NormalizeData(seurat.NK.unint, verbose = FALSE)
seurat.NK.unint <- FindVariableFeatures(seurat.NK.unint, nfeatures = n.features, verbose = FALSE)
seurat.NK.unint <- ScaleData(seurat.NK.unint, verbose = FALSE)
seurat.NK.unint <- RunPCA(seurat.NK.unint, npcs = n.dims, verbose = FALSE)
seurat.NK.unint <- FindNeighbors(seurat.NK.unint, dims = 1:n.dims, verbose = FALSE)
seurat.NK.unint <- FindClusters(seurat.NK.unint, resolution = resolution, verbose = FALSE)
seurat.NK.unint <- RunUMAP(seurat.NK.unint, dims = 1:n.dims, reduction = "pca", 
                           reduction.name = "umap.unintegrated", verbose = FALSE)

print("Unintegrated baseline created")

### Integration Workflow ###
print(paste0("Using integration method: ", integration.method))

# Calculate safe integration parameters
sample_counts <- table(seurat.NK$Sample.id)
min_sample_size <- min(sample_counts)

k.weight <- min(100, floor(min_sample_size * 0.66))
k.weight <- max(k.weight, 10)

k.anchor <- min(5, floor(min_sample_size * 0.1))
k.anchor <- max(k.anchor, 3)

k.filter <- min(200, floor(min_sample_size * 0.75))
k.filter <- max(k.filter, 10)

k.score <- min(30, floor(min_sample_size * 0.45))
k.score <- max(k.score, 5)

print(paste0("Integration parameters (min sample size = ", min_sample_size, "):"))
print(paste0("  k.weight = ", k.weight))
print(paste0("  k.anchor = ", k.anchor))
print(paste0("  k.filter = ", k.filter))
print(paste0("  k.score = ", k.score))

# OPTION A: Seurat v3 style 
if (integration.method == "v3") {
  
  print("Using Seurat v3 integration (original method)...")
  
  # Split into list
  seurat.NK.list <- SplitObject(seurat.NK, split.by = "Sample.id")
  
  # Normalize and find features
  seurat.NK.list <- lapply(seurat.NK.list, function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, nfeatures = n.features, verbose = FALSE)
    return(x)
  })
  
  # Find anchors
  features <- SelectIntegrationFeatures(object.list = seurat.NK.list, nfeatures = n.features)
  anchors <- FindIntegrationAnchors(object.list = seurat.NK.list, 
                                    anchor.features = features,
                                    dims = 1:n.dims,
                                    k.anchor = k.anchor,
                                    k.filter = k.filter,
                                    k.score = k.score)
  
  # Integrate
  seurat.NK.integrated <- IntegrateData(anchorset = anchors, 
                                        dims = 1:n.dims,
                                        k.weight = k.weight)
  
  DefaultAssay(seurat.NK.integrated) <- "integrated"
  integration.reduction <- "pca"
  
} else {
  
  # OPTION B: Seurat v5 style
  print("Using Seurat v5 integration...")
  
  # Split layers
  seurat.NK[["RNA"]] <- split(seurat.NK[["RNA"]], f = seurat.NK$Sample.id)
  
  # Normalize
  seurat.NK <- NormalizeData(seurat.NK, verbose = FALSE)
  seurat.NK <- FindVariableFeatures(seurat.NK, nfeatures = n.features, verbose = FALSE)
  seurat.NK <- ScaleData(seurat.NK, verbose = FALSE)
  seurat.NK <- RunPCA(seurat.NK, npcs = n.dims, verbose = FALSE)
  
  # Integrate
  if (integration.method == "CCA") {
    seurat.NK.integrated <- IntegrateLayers(
      object = seurat.NK,
      method = CCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.cca",
      dims = 1:n.dims,
      k.anchor = k.anchor,
      k.filter = k.filter,
      k.score = k.score,
      k.weight = k.weight,
      verbose = TRUE
    )
    integration.reduction <- "integrated.cca"
    
  } else if (integration.method == "RPCA") {
    seurat.NK.integrated <- IntegrateLayers(
      object = seurat.NK,
      method = RPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.rpca",
      dims = 1:n.dims,
      k.anchor = k.anchor,
      k.filter = k.filter,
      k.score = k.score,
      k.weight = k.weight,
      verbose = TRUE
    )
    integration.reduction <- "integrated.rpca"
    
  } else if (integration.method == "Harmony") {
    library(harmony)
    seurat.NK.integrated <- IntegrateLayers(
      object = seurat.NK,
      method = HarmonyIntegration,
      orig.reduction = "pca",
      new.reduction = "harmony",
      group.by.vars = "Sample.id",
      dims = 1:n.dims,
      verbose = TRUE
    )
    integration.reduction <- "harmony"
  }
}

# Scale integrated data
seurat.NK.integrated <- ScaleData(seurat.NK.integrated, verbose = FALSE)

# Run PCA (if not using v5 with new reduction)
if (integration.method == "v3") {
  seurat.NK.integrated <- RunPCA(seurat.NK.integrated, npcs = n.dims, verbose = FALSE)
}

# Clustering
seurat.NK.integrated <- FindNeighbors(seurat.NK.integrated, 
                                      reduction = integration.reduction, 
                                      dims = 1:n.dims, verbose = FALSE)
seurat.NK.integrated <- FindClusters(seurat.NK.integrated, 
                                     resolution = resolution, verbose = FALSE)

# UMAP
seurat.NK.integrated <- RunUMAP(seurat.NK.integrated, 
                                reduction = integration.reduction, 
                                dims = 1:n.dims, verbose = FALSE)

# Switch to RNA for visualization
DefaultAssay(seurat.NK.integrated) <- "RNA"

### Comparison: Integrated vs Unintegrated ###
# Visualizations
p1_unint_sample <- DimPlot(seurat.NK.unint, reduction = "umap.unintegrated", 
                           group.by = "Sample.id") +
  ggtitle("Unintegrated - by Sample")

p2_int_sample <- DimPlot(seurat.NK.integrated, reduction = "umap", 
                         group.by = "Sample.id") +
  ggtitle(paste0(integration.method, " Integrated - by Sample"))

p3_unint_celltype <- DimPlot(seurat.NK.unint, reduction = "umap.unintegrated", 
                             group.by = "Immune_labelled", cols = CT.colors) +
  ggtitle("Unintegrated - by Cell Type")

p4_int_celltype <- DimPlot(seurat.NK.integrated, reduction = "umap", 
                           group.by = "Immune_labelled", cols = CT.colors) +
  ggtitle(paste0(integration.method, " Integrated - by Cell Type"))

# Save comparison
pdf(paste0(Output.QC.dir, "Integration_Comparison_", integration.method, ".pdf"), 
    width = 16, height = 12)
print((p1_unint_sample | p2_int_sample) / (p3_unint_celltype | p4_int_celltype))
dev.off()

### Join Layers and Save ###
# Join layers
if (integration.method != "v3") {
  seurat.NK.integrated[["RNA"]] <- JoinLayers(seurat.NK.integrated[["RNA"]])
}

# Save
SaveH5Seurat(seurat.NK.integrated, 
             filename = paste0(Output.seurat.dir, "Endometrium_uNK_integrated_", 
                               integration.method, ".h5seurat"),
             overwrite = TRUE)

saveRDS(seurat.NK.integrated, 
        file = paste0(Output.seurat.dir, "Endometrium_uNK_integrated_", 
                      integration.method, ".rds"))

print("Integration complete!")
print(paste0("Method: ", integration.method))
print(paste0("Final object contains ", ncol(seurat.NK.integrated), " cells"))
print(table(seurat.NK.integrated$Immune_labelled))

# Clean up
rm(seurat.NK.unint)
gc()
