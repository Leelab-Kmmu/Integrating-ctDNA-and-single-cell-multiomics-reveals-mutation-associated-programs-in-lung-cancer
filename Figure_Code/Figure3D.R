library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ====================================================
# 1. Order of the seven cell types
# ====================================================

cell_order <- c(
  "AT1",
  "AT2",
  "Tumor",
  "Endothelial",
  "Fibroblast",
  "Lymphoid",
  "Myeloid"
)

# ====================================================
# 2. Extract ATAC counts
# ====================================================

atac_counts <- GetAssayData(
  lung_multi,
  assay = "ATAC",
  layer = "counts"
)

# Ensure that TF peaks are present in the ATAC matrix
tf_data <- tf_peak_result[
  tf_peak_result$Peak %in% rownames(atac_counts),
  c("TF", "Peak"),
  drop = FALSE
]

tf_data <- unique(tf_data)

# ====================================================
# 3. Calculate the detection rate of each peak across the seven cell types
# ====================================================

celltype <- lung_multi$celltype
names(celltype) <- colnames(lung_multi)

celltype <- celltype[colnames(atac_counts)]

group_mat <- model.matrix(
  ~ 0 + factor(celltype, levels = cell_order)
)

colnames(group_mat) <- cell_order

binary_atac <- atac_counts
binary_atac@x[] <- 1

peak_fraction_all <- as.matrix(
  binary_atac %*% group_mat
)

peak_fraction_all <- sweep(
  peak_fraction_all,
  2,
  colSums(group_mat),
  "/"
)

# ====================================================
# 4. Extract the maximum detection rate among peaks associated with each TF
# ====================================================

tf_data$Max_accessibility <- apply(
  peak_fraction_all[
    tf_data$Peak,
    ,
    drop = FALSE
  ],
  1,
  max
)

tf_data$Mean_accessibility <- rowMeans(
  peak_fraction_all[
    tf_data$Peak,
    ,
    drop = FALSE
  ]
)

# Cell type with the highest accessibility for each peak
tf_data$Highest_celltype <- apply(
  peak_fraction_all[
    tf_data$Peak,
    ,
    drop = FALSE
  ],
  1,
  function(x) names(which.max(x))
)

# ====================================================
# 5. Select the most accessible peak for each TF
# ====================================================

tf_data <- tf_data[
  order(
    tf_data$TF,
    -tf_data$Max_accessibility,
    -tf_data$Mean_accessibility
  ),
  ,
  drop = FALSE
]

tf_best <- tf_data[
  !duplicated(tf_data$TF),
  ,
  drop = FALSE
]

# Retain all TFs, including TFs without peaks
tf_best_all <- merge(
  data.frame(TF = tf_genes),
  tf_best,
  by = "TF",
  all.x = TRUE,
  sort = FALSE
)

print(as.data.frame(tf_best_all))

write.csv(
  tf_best_all,
  "E:/Lung/TF_highest_accessibility_peaks.csv",
  row.names = FALSE
)







# ====================================================
# 6. Merge the original heatmap peaks with TF peaks
# ====================================================

original_peaks <- colnames(mat_plot)

tf_peaks <- unique(
  na.omit(tf_best_all$Peak)
)

all_peaks <- unique(
  c(original_peaks, tf_peaks)
)

# ====================================================
# 7. Recalculate the heatmap matrix
# ====================================================

mat <- peak_fraction_all[
  all_peaks,
  cell_order,
  drop = FALSE
]

mat <- log1p(mat * 100)

# Standardize each peak
mat_z <- t(scale(t(mat)))

mat_z[!is.finite(mat_z)] <- 0

mat_z <- pmax(
  pmin(mat_z, 2),
  -2
)

# Transpose
mat_plot_new <- t(mat_z)

# ====================================================
# 8. Arrange peaks by the cell type with the highest accessibility
# ====================================================

peak_group_new <- apply(
  peak_fraction_all[
    all_peaks,
    cell_order,
    drop = FALSE
  ],
  1,
  function(x) names(which.max(x))
)

peak_group_new <- factor(
  peak_group_new,
  levels = cell_order
)

ord <- order(peak_group_new)

mat_plot_new <- mat_plot_new[
  ,
  ord,
  drop = FALSE
]

peak_group_new <- peak_group_new[ord]

# Statistics
print(table(peak_group_new))



# ====================================================
# 9. Match TFs to peak positions in the heatmap
# ====================================================

tf_anno <- tf_best_all[
  !is.na(tf_best_all$Peak),
  ,
  drop = FALSE
]

tf_anno$Column <- match(
  tf_anno$Peak,
  colnames(mat_plot_new)
)

tf_anno <- tf_anno[
  !is.na(tf_anno$Column),
  ,
  drop = FALSE
]

# If multiple TFs correspond to the same peak,
# merge the TF names to avoid duplicate anno_mark positions
tf_anno2 <- aggregate(
  TF ~ Column,
  data = tf_anno,
  FUN = function(x) paste(unique(x), collapse = "/")
)

tf_anno2 <- tf_anno2[
  order(tf_anno2$Column),
  ,
  drop = FALSE
]

# ====================================================
# 10. Set the top TF annotations
# ====================================================

top_anno <- HeatmapAnnotation(
  
  TF = anno_mark(
    at = tf_anno2$Column,
    labels = tf_anno2$TF,
    side = "top",
    
    labels_gp = gpar(
      fontsize = 8,
      fontface = "italic"
    ),
    
    link_gp = gpar(
      col = "black",
      lwd = 0.7
    ),
    
    extend = unit(25, "mm")
  ),
  
  Peak_group = peak_group_new,
  
  col = list(
    Peak_group = c(
      "AT1" = "#E84D3C",
      "AT2" = "#E6A51B",
      "Tumor" = "#CC4794",
      "Endothelial" = "#4BA184",
      "Fibroblast" = "#2A9D8F",
      "Lymphoid" = "#3B78BE",
      "Myeloid" = "#6657AE"
    )
  ),
  
  show_annotation_name = FALSE
)

# ====================================================
# 11. Draw the heatmap
# ====================================================

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#267BBF", "#F5F5F5", "#C92827")
)

ht <- Heatmap(
  mat_plot_new,
  
  name = "Column Z-score",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  # Do not use column_split, ensuring correct anno_mark positions
  show_column_names = FALSE,
  
  top_annotation = top_anno,
  
  show_row_names = TRUE,
  row_names_side = "right",
  
  row_names_gp = gpar(fontsize = 12),
  
  row_title = NULL,
  column_title = NULL,
  
  use_raster = TRUE,
  
  heatmap_legend_param = list(
    title = "Column Z-score",
    at = c(-2, 0, 2),
    legend_direction = "horizontal"
  )
)

draw(
  ht,
  heatmap_legend_side = "bottom"
)



pdf(
  "E:/Lung/ATAC_Celltype_Peak_Heatmap2.pdf",
  width = 13,
  height = 5.5
)

draw(
  ht,
  heatmap_legend_side = "bottom"
)

dev.off()
