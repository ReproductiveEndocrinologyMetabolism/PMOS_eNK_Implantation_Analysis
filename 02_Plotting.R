# Author: Gustaw Eriksson
# Date: 2025-11-27
# Description: Generate publication-ready plots for eNK cell analysis

# Loading the packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(cowplot)
library(openxlsx)
display.brewer.pal(8, "Set2")

### Setting variables section ###
# Set paths
seurat.rds.path = "Data/Subset_seurat_objects/Endometrium_uNK_integrated_CCA.rds"
output.dir = "Output/Publication_Figures/"

# Create output directory if it doesn't exist
if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
  print(paste0("Created output directory: ", output.dir))
}

# Define marker gene sets
# eNK validation markers in supplementary 
eNK_markers = c("NCAM1", # eNK CD56
                "FGFBP2", # Negative CD16
                "ENTPD1", # eNK 1 CD39
                "ITGB2", # eNK 2 CD18
                "CD160", # eNK 2 C160
                "ITGAE"# eNK 3 CD103
                ) 

# Markers of interest
TGFb_markers = c(
  "TGFB1",      # Ligand production - HIGH
  "TGFBR1",     # Receptor - HIGH
  "TGFBR2",     # Receptor - subset marker?
  "SMAD2",      # Signaling - HIGH
  "SMAD3",
  "SMAD4",
  "SMAD7"      # Negative regulator - moderate
)

TGFb_markers_short = c(
  "TGFB1",      # Ligand production - HIGH
  "TGFBR1",     # Receptor - HIGH
  "TGFBR2"       # TGF-β activation - moderate
)

# Define color palettes
group.colors = c("#A0A0A0", "#B0D8D6","#95BFE1", "#65D46E")

cluster.colors = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                    "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

CT.colors = brewer.pal(3, "Set2")

Coexp.colors = c("#FFFFCC", "#CC0000", "#3399FF")

### Function definitions ###

#' Load the integrated eNK Seurat object
load_seurat_object = function(rds.path) {
  if (!file.exists(rds.path)) {
    stop(paste0("File not found: ", rds.path))
  }
  seurat.obj = readRDS(rds.path)
  print(paste0("Loaded Seurat object with ", ncol(seurat.obj), " cells"))
  print("Cell type composition:")
  print(table(seurat.obj$Immune_labelled))
  return(seurat.obj)
}

#' Generate UMAP plots colored by different metadata
plot_umap = function(seurat.obj, group.x = "Immune_labelled", 
                      colors = NULL, label = TRUE, pt.size = 0.5) {
  
  plot.x = DimPlot(seurat.obj, reduction = "umap", group.by = group.x,
               label = label, pt.size = pt.size, repel = TRUE,
              cols = colors) +
    theme_cowplot() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0.5)
    ) +
    labs(title = paste0("UMAP colored by ", group.x),
         x = "UMAP 1",
         y = "UMAP 2")
  
  return(plot.x)
}

plot_umap_split = function(seurat.obj, group.x = "Immune_labelled",
                            split.x = "Group",
                            colors = NULL, label = TRUE, pt.size = 0.5,
                            order.splits = NULL, ncol.x = NULL) {
  
  # Order splits if specified
  if (!is.null(order.splits)) {
    seurat.obj@meta.data[[split.x]] = factor(seurat.obj@meta.data[[split.x]], 
                                              levels = order.splits)
  }
  
  # If ncol not specified, auto-calculate based on number of splits
  if (is.null(ncol.x)) {
    n_splits = length(unique(seurat.obj@meta.data[[split.x]]))
    ncol.x = min(n_splits, 3)  # Max 3 columns
  }
  
  plot.x = DimPlot(seurat.obj, 
              reduction = "umap", 
              group.by = group.x,
              split.by = split.x,
              label = label, 
              pt.size = pt.size, 
              repel = TRUE,
              cols = colors,
              ncol = ncol.x) +
    theme_cowplot() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12),
      plot.title = element_text(size = 14, hjust = 0.5),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_blank()) +
    labs(title = paste0("UMAP colored by ", group.x, " (split by ", split.x, ")"),
         x = "UMAP 1",
         y = "UMAP 2")
  
  return(plot.x)
  
}

# Generate UMAP plots showing gene co-expression
plot_coexpression_umap = function(seurat.obj, gene.x, gene.y, 
                                   pt.size = 0.5, blend.threshold = 0.25,
                                   order = TRUE, 
                                   colors = c("green3", "purple", "yellow")) {
  
  # Check if genes exist in the dataset
  if (!gene.x %in% rownames(seurat.obj)) {
    stop(paste0("Gene ", gene.x, " not found in the dataset"))
  }
  if (!gene.y %in% rownames(seurat.obj)) {
    stop(paste0("Gene ", gene.y, " not found in the dataset"))
  }
  
  # Handle color inputs
  if (length(colors) == 2) {
    # If only 2 colors provided, add yellow as default co-expression color
    colors = c(colors[1], colors[2], "yellow")
    message("Using yellow as co-expression color. Provide 3 colors to customize: c(gene.x, gene.y, coexpression)")
  } else if (length(colors) != 3) {
    warning("Colors should be a vector of length 2 or 3. Using default: green, purple, yellow")
    colors = c("green3", "purple", "yellow")
  }
  
  # Create blend plot (shows co-expression)
  p_blend = FeaturePlot(seurat.obj, 
                         features = c(gene.x, gene.y),
                         blend = TRUE,
                         blend.threshold = blend.threshold,
                         pt.size = pt.size,
                         order = order,
                         cols = colors) +
    theme_cowplot() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold.italic", hjust = 0.5)
    ) +
    labs(title = paste0("Co-expression: ", gene.x, " & ", gene.y),
         x = gene.x,
         y = gene.y)
  
  return(p_blend)
}

#' Generate feature plots for marker genes
plot_features = function(seurat.obj, features, ncol = 3, 
                          pt.size = 0.5, order = TRUE) {
  plot.x = FeaturePlot(seurat.obj, features = features, 
                   ncol = ncol, pt.size = pt.size, order = order) &
    theme_cowplot() &
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 12, face = "bold.italic"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11)
    ) &
    scale_color_viridis(option = "magma", direction = -1) &
    labs(x = "UMAP 1", y = "UMAP 2")
  
  return(plot.x)
}

#' Generate violin plots for marker genes
plot_violins = function(seurat.obj, features, group.by = "Immune_labelled",
                        split.by = NULL, colors = NULL, ncol = 3) {
  plot.x = VlnPlot(seurat.obj, features = features, group.by = group.by, split.by = split.by,
               ncol = ncol, pt.size = 0, cols = colors) &
    theme_cowplot() &
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 11),
      plot.title = element_text(size = 12, face = "bold.italic")
    )
  
  if (!is.null(colors)) {
    plot.x = plot.x & scale_fill_manual(values = colors)
  }
  
  return(plot.x)
}

# Generate dot plots for marker genes
plot_dotplot = function(seurat.obj, features, group.by = "Immune_labelled",
                         split.by = NULL, scale = TRUE, dot.min = 0,
                         colors = c("lightgrey", "blue"),
                         reverse_y = FALSE) {
  
  if (!is.null(split.by)) {
    # Create a combined identity for proper visualization with color legend
    seurat.obj$temp_combined = paste(seurat.obj[[group.by]][,1], 
                                      seurat.obj[[split.by]][,1], 
                                      sep = "_")
    
    # Reverse factor levels if requested
    if (reverse_y) {
      seurat.obj$temp_combined = factor(seurat.obj$temp_combined,
                                         levels = rev(sort(unique(seurat.obj$temp_combined))))
    }
    
    plot.x = DotPlot(seurat.obj, features = features, 
                 group.by = "temp_combined",
                 cols = colors, 
                 dot.scale = 8,
                 dot.min = dot.min,
                 scale = scale) +
      theme_cowplot() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "italic"),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold"),
        legend.position = "right"
      ) +
      labs(x = "Genes", y = "Cell Type", 
           color = "Average\nExpression", 
           size = "Percent\nExpressed")
    
    # Clean up temporary column
    seurat.obj$temp_combined = NULL
    
  } else {
    # Reverse factor levels if requested
    if (reverse_y) {
      seurat.obj[[group.by]] = factor(seurat.obj[[group.by]][,1],
                                       levels = rev(sort(unique(seurat.obj[[group.by]][,1]))))
    }
    
    plot.x = DotPlot(seurat.obj, features = features, group.by = group.by,
                 scale = scale,
                 dot.min = dot.min,
                 cols = colors, col.min = -2.5, col.max = 2.5) +
      theme_cowplot() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "italic"),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11, face = "bold")
      ) +
      labs(x = "Genes", y = "Cell Type")
  }
  
  return(plot.x)
}

# Save a plot with high resolution for publication
save_publication_plot = function(plot, filename, width = 10, height = 8, dpi = 300) {
  ggsave(filename = paste0(filename, ".pdf"), plot = plot, width = width, height = height, 
         dpi = dpi, units = "in", bg = "white")
  ggsave(filename = paste0(filename, ".png"), plot = plot, width = width, height = height, 
         dpi = dpi, units = "in", bg = "white")
  print(paste0("Saved plot: ", filename))
}

### Main analysis workflow ###

# Load the Seurat object
seurat.NK = load_seurat_object(seurat.rds.path)

# Rename uNK to eNK
seurat.NK$Immune_labelled = gsub("uNK", "eNK", seurat.NK$Immune_labelled)
table(seurat.NK$Immune_labelled)

# Rename PCOS W0 to PCOS
seurat.NK$Group = gsub("PCOS_W0", "PCOS", seurat.NK$Group)
table(seurat.NK$Group)

# Set default assay to RNA for visualization
DefaultAssay(seurat.NK) = "RNA"

### Generate individual plots ###
# 1. UMAP plots
print("Generating UMAP plots...")
umap_celltype = plot_umap(seurat.NK, group.x = "Immune_labelled", label = TRUE, colors = CT.colors)
umap_group = plot_umap(seurat.NK, group.x = "Group", colors = group.colors, label = FALSE)
umap_clusters = plot_umap(seurat.NK, group.x = "seurat_clusters", 
                           colors = NULL, label = TRUE)

save_publication_plot(umap_celltype, paste0(output.dir, "UMAP_CellType"), 
                      width = 4, height = 3)
save_publication_plot(umap_group, paste0(output.dir, "UMAP_Group"), 
                      width = 4, height = 3)
save_publication_plot(umap_clusters, paste0(output.dir, "UMAP_Clusters"), 
                      width = 4, height = 3)

umap_split_group = plot_umap_split(seurat.NK, group.x = "Group", split.x = "Group", 
                                   colors = group.colors, label = FALSE)
umap_split_immune = plot_umap_split(seurat.NK, group.x = "Immune_labelled", split.x = "Group", 
                                   colors = CT.colors, label = FALSE)
save_publication_plot(umap_split_group, paste0(output.dir, "UMAP_Split_Group"), 
                      width = 7, height = 3)
save_publication_plot(umap_split_immune, paste0(output.dir, "UMAP_Split_Immune"), 
                      width = 7, height = 3)

coexp_umap = plot_coexpression_umap(seurat.NK, 
                                     gene.x = "TGFB1", 
                                     gene.y = "TGFBR1",
                                     colors = Coexp.colors,
                                     blend.threshold = 0.25)
save_publication_plot(coexp_umap, paste0(output.dir, "UMAP_CoExp_TGFB1_TGFBR1"), 
                      width = 14, height = 3)
coexp_umap = plot_coexpression_umap(seurat.NK, 
                                     gene.x = "TGFB1", 
                                     gene.y = "TGFBR2",
                                     colors = Coexp.colors,
                                     blend.threshold = 0.25)
save_publication_plot(coexp_umap, paste0(output.dir, "UMAP_CoExp_TGFB1_TGFBR2"), 
                      width = 14, height = 3)
coexp_umap = plot_coexpression_umap(seurat.NK, 
                                     gene.x = "TGFBR1", 
                                     gene.y = "TGFBR2",
                                     colors = Coexp.colors,
                                     blend.threshold = 0.25)
save_publication_plot(coexp_umap, paste0(output.dir, "UMAP_CoExp_TGFBR1_TGFBR2"), 
                      width = 14, height = 3)


# 2. Feature plots for each eNK subset markers
print("Generating feature plots...")
features_eNK_markers = plot_features(seurat.NK, features = eNK_markers, ncol = 3)
features_eNK_markers
save_publication_plot(features_eNK_markers, paste0(output.dir, "FeaturePlot_eNK_markers"),
                      width = 10, height = 6)

features_TGFb_markers = plot_features(seurat.NK, features = TGFb_markers, ncol = 4)
features_TGFb_markers
save_publication_plot(features_TGFb_markers, paste0(output.dir, "FeaturePlot_TGFb_markers"),
                      width = 12.5, height = 6)

features_TGFb_short_markers = plot_features(seurat.NK, features = TGFb_markers_short, ncol = 3)
features_TGFb_short_markers
save_publication_plot(features_TGFb_short_markers, paste0(output.dir, "FeaturePlot_TGFb_LR_markers"),
                      width = 10, height = 3)

# 3. Violin plots
print("Generating violin plots...")
violin_eNK_markers = plot_violins(seurat.NK, features = eNK_markers, colors = CT.colors,
                                  group.by = "Immune_labelled")
save_publication_plot(violin_eNK_markers, paste0(output.dir, "VlnPlot_eNK_markers"),
                      width = 10, height = 4)

violin_TGFb_markers = plot_violins(seurat.NK, features = TGFb_markers, colors = CT.colors,
                                  group.by = "Immune_labelled", ncol = 4)
save_publication_plot(violin_TGFb_markers, paste0(output.dir, "VlnPlot_TGFb_markers"),
                      width = 10, height = 4)
violin_TGFb_markers = plot_violins(seurat.NK, features = TGFb_markers, colors = group.colors,
                                   group.by = "Immune_labelled", split.by = "Group", ncol = 4)
save_publication_plot(violin_TGFb_markers, paste0(output.dir, "VlnPlot_TGFb_markers_split"),
                      width = 15, height = 4)

violin_TGFb_markers = plot_violins(seurat.NK, features = TGFb_markers_short, colors = CT.colors,
                                   group.by = "Immune_labelled", ncol = 3)
save_publication_plot(violin_TGFb_markers, paste0(output.dir, "VlnPlot_TGFb_LR_markers"),
                      width = 10, height = 3)
violin_TGFb_markers = plot_violins(seurat.NK, features = TGFb_markers_short, colors = group.colors,
                                   group.by = "Immune_labelled", split.by = "Group", ncol = 3)
save_publication_plot(violin_TGFb_markers, paste0(output.dir, "VlnPlot_TGFb_LR_markers_split"),
                      width = 16, height = 3)

# 4. Dot plots
print("Generating dot plots...")
dot_eNK_markers = plot_dotplot(seurat.NK, features = eNK_markers,
                                  group.by = "Immune_labelled",
                               colors = c("#FFEDA0", "#F03B20"),
                               scale = TRUE, reverse_y = TRUE)
save_publication_plot(dot_eNK_markers, paste0(output.dir, "Dotplot_eNK_markers"),
                      width = 4, height = 3.5)

dot_TGFb_markers = plot_dotplot(seurat.NK, features = TGFb_markers,
                                   group.by = "Immune_labelled",
                                colors = c("#FFEDA0", "#F03B20"),
                                scale = FALSE, reverse_y = TRUE)
save_publication_plot(dot_TGFb_markers, paste0(output.dir, "Dotplot_TGFb_markers"),
                      width = 4, height = 4)

dot_TGFb_markers = plot_dotplot(seurat.NK, features = TGFb_markers, scale = FALSE, 
                                colors = c("#FFEDA0", "#F03B20"), reverse_y = TRUE,
                                   group.by = "Immune_labelled", split.by = "Group")
save_publication_plot(dot_TGFb_markers, paste0(output.dir, "Dotplot_TGFb_markers_split"),
                      width = 6, height = 4)

dot_TGFb_markers = plot_dotplot(seurat.NK, features = TGFb_markers_short, reverse_y = TRUE,
                                   group.by = "Immune_labelled", scale = FALSE,
                                colors = c("#FFEDA0", "#F03B20"),)
save_publication_plot(dot_TGFb_markers, paste0(output.dir, "Dotplot_TGFb_LR_markers"),
                      width = 4, height = 4)
dot_TGFb_markers = plot_dotplot(seurat.NK, features = TGFb_markers_short, scale = FALSE,
                                colors = c("#FFEDA0", "#F03B20"), reverse_y = TRUE,
                                   group.by = "Immune_labelled", split.by = "Group")
save_publication_plot(dot_TGFb_markers, paste0(output.dir, "Dotplot_TGFb_LR_markers_split"),
                      width = 6, height = 4)

# 5. Summary statistics table
print("Generating summary statistics...")
summary_stats = data.frame(
  Cell_Type = names(table(seurat.NK$Immune_labelled)),
  Total_Cells = as.vector(table(seurat.NK$Immune_labelled)),
  Percentage = round(as.vector(prop.table(table(seurat.NK$Immune_labelled))) * 100, 2)
)

group_summary = as.data.frame.matrix(table(seurat.NK$Immune_labelled, seurat.NK$Group))
summary_stats = cbind(summary_stats, group_summary)

write.csv(summary_stats, paste0(output.dir, "eNK_cell_summary_statistics.csv"),
          row.names = FALSE)
write.xlsx(summary_stats, paste0(output.dir, "eNK_cell_summary_statistics.xlsx"))