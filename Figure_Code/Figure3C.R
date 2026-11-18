library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)

# ============================================================
# 1. Define marker genes
# ============================================================

markers <- list(
  AT1 = c(
    "RTKN2", "AGER", "CAV1", "CLIC5"
  ),
  
  AT2 = c(
    "SFTPC", "SFTPB", "SFTPA1",
    "SFTPA2", "ABCA3", "NAPSA", "ETV5"
  ),
  
  Tumor = c(
    "KRT18", "MUC1", "CEMIP",
    "MUC16", "SPINK1",
    "CEACAM5", "CEACAM6",
    "LAMC2"
  ),
  
  Endothelial = c(
    "PECAM1", "FLT1", "CLDN5"
  ),
  
  Fibroblast = c(
    "COL3A1", "COL1A1", "PDGFRA"
  ),
  
  Lymphoid = c(
    "PTPRC", "CD247", "CD3D", "TRAC"
  ),
  
  Myeloid = c(
    "SPI1", "LYZ", "MARCO", "MS4A7"
  )
)

# ============================================================
# 2. Check which genes are absent
# ============================================================

all_markers <- unlist(markers)

missing_genes <- setdiff(
  all_markers,
  rownames(lung_RNA)
)

missing_genes


markers <- lapply(
  markers,
  function(x) x[x %in% rownames(lung_RNA)]
)

markers



# Retrieve again
plot_df <- p0$data

colnames(plot_df)


plot_df$gene         <- as.character(plot_df$features.plot)
plot_df$Celltype     <- as.character(plot_df$id)
plot_df$AvgExp       <- plot_df$avg.exp.scaled
plot_df$PctExp       <- plot_df$pct.exp
plot_df$marker_group <- as.character(plot_df$feature.groups)




cell_levels <- c(
  "Myeloid",
  "Lymphoid",
  "Fibroblast",
  "Endothelial",
  "Tumor",
  "AT2",
  "AT1"
)

marker_levels <- c(
  "AT1",
  "AT2",
  "Tumor",
  "Endothelial",
  "Fibroblast",
  "Lymphoid",
  "Myeloid"
)

plot_df$Celltype <- factor(
  plot_df$Celltype,
  levels = cell_levels
)

plot_df$marker_group <- factor(
  plot_df$marker_group,
  levels = marker_levels
)

plot_df$gene <- factor(
  plot_df$gene,
  levels = unlist(markers)
)



library(ggplot2)
library(grid)

p_marker <- ggplot(
  plot_df,
  aes(
    x = gene,
    y = Celltype
  )
) +
  geom_point(
    aes(
      size = PctExp,
      color = AvgExp
    )
  ) +
  facet_grid(
    cols = vars(marker_group),
    scales = "free_x",
    space = "free_x"
  ) +
  scale_size_continuous(
    name = "Percent Expressed",
    range = c(0.5, 6),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  scale_color_gradient2(
    name = "Average Expression",
    low = "#8E8CD8",
    mid = "#F7F7F7",
    high = "#C62828",
    midpoint = 0
  ) +
  labs(
    title = "Cell-type marker expression",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0
    ),
    
    strip.background = element_rect(
      fill = "#263238",
      color = "#263238"
    ),
    
    strip.text = element_text(
      color = "white",
      face = "bold",
      size = 8
    ),
    
    panel.border = element_rect(
      color = "grey30",
      fill = NA,
      linewidth = 0.8
    ),
    
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.6
    ),
    
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_text(
      angle = 50,
      hjust = 1,
      vjust = 1,
      size = 8,
      color = "grey25"
    ),
    
    axis.text.y = element_text(
      size = 9,
      face = "bold",
      color = "grey25"
    ),
    
    axis.ticks = element_blank(),
    
    legend.position = "right",
    
    panel.spacing.x = unit(0.15, "cm")
  )

p_marker


ggsave(
  filename = "E:/Lung/lung_RNA_celltype_marker_DotPlot.pdf",
  plot = p_marker,
  width = 15,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)
