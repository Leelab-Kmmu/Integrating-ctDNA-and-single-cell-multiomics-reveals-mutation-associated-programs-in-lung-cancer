Reductions(lung_multi)


lung_multi <- FindMultiModalNeighbors(
  lung_multi,
  reduction.list = list("pca", "lsi.atac"),
  dims.list = list(
    1:30,   # RNA PCA
    3:30    # ATAC LSI, based on the preceding DepthCor results
  ),
  modality.weight.name = "RNA.weight"
)


lung_multi <- RunUMAP(
  lung_multi,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)


DimPlot(
  lung_multi,
  reduction = "wnn.umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3
)


head(
  lung_multi@meta.data[
    , c("integrated.weight", "ATAC.weight")
  ]
)

summary(lung_multi$integrated.weight)
summary(lung_multi$ATAC.weight)

cor(
  lung_multi$integrated.weight,
  1 - lung_multi$ATAC.weight
)


library(dplyr)

weight_summary <- lung_multi@meta.data %>%
  group_by(celltype) %>%
  summarise(
    n = n(),
    RNA_weight_mean = mean(integrated.weight),
    RNA_weight_median = median(integrated.weight),
    ATAC_weight_mean = mean(ATAC.weight),
    ATAC_weight_median = median(ATAC.weight)
  )

weight_summary








library(dplyr)
library(ggplot2)

plot_df <- weight_summary %>%
  arrange(desc(ATAC_weight_mean)) %>%
  mutate(
    celltype = factor(
      celltype,
      levels = celltype
    )
  )

ggplot(
  plot_df,
  aes(
    x = celltype,
    y = ATAC_weight_mean
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sprintf("%.3f", ATAC_weight_mean)),
    vjust = -0.4,
    size = 3.5
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Mean ATAC modality weight"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  expand_limits(
    y = max(plot_df$ATAC_weight_mean) * 1.15
  )



library(tidyr)

plot_df2 <- weight_summary %>%
  select(
    celltype,
    RNA_weight_mean,
    ATAC_weight_mean
  ) %>%
  pivot_longer(
    cols = c(
      RNA_weight_mean,
      ATAC_weight_mean
    ),
    names_to = "Modality",
    values_to = "Weight"
  ) %>%
  mutate(
    Modality = recode(
      Modality,
      RNA_weight_mean = "RNA",
      ATAC_weight_mean = "ATAC"
    )
  )

ggplot(
  plot_df2,
  aes(
    x = celltype,
    y = Weight,
    fill = Modality
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.7
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Mean modality weight",
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )














library(Seurat)
library(ggplot2)
library(data.table)
library(scales)
library(grid)

# ============================================================
# 1. Extract the WNN UMAP
# ============================================================

umap_df <- as.data.frame(
  Embeddings(
    lung_multi,
    reduction = "wnn.umap"
  )
)

colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$Cell <- rownames(umap_df)

umap_df$Celltype <- lung_multi$celltype[
  match(
    umap_df$Cell,
    colnames(lung_multi)
  )
]

umap_df <- as.data.table(umap_df)

# ============================================================
# 2. Order
# ============================================================

cell_levels <- c(
  "AT1",
  "AT2",
  "Tumor",
  "Endothelial",
  "Fibroblast",
  "Lymphoid",
  "Myeloid"
)

umap_df[, Celltype := factor(
  Celltype,
  levels = cell_levels
)]

# ============================================================
# 3. Color palette
# ============================================================

cell_cols <- c(
  "AT1"         = "#E84D3C",
  "AT2"         = "#E6A51B",
  "Tumor"       = "#CC4794",
  "Endothelial" = "#4BA184",
  "Fibroblast"  = "#2A9D8F",
  "Lymphoid"    = "#3B78BE",
  "Myeloid"     = "#6657AE"
)

# ============================================================
# 4. Legend labels
# ============================================================

count_dt <- umap_df[, .N, by = Celltype]

count_dt[, legend_label := paste0(
  Celltype,
  "  (",
  comma(N),
  ")"
)]

legend_labels <- setNames(
  count_dt$legend_label,
  as.character(count_dt$Celltype)
)

# ============================================================
# 5. Calculate coordinate ranges
# ============================================================

xr <- range(umap_df$UMAP_1, na.rm = TRUE)
yr <- range(umap_df$UMAP_2, na.rm = TRUE)

dx <- diff(xr)
dy <- diff(yr)

# ============================================================
# 6. Calculate the center coordinates of each cell type
# ============================================================

center_dt <- umap_df[
  ,
  .(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  ),
  by = Celltype
]

print(center_dt)

# ============================================================
# 7. Label positions
# ============================================================

label_df <- copy(center_dt)

label_df[Celltype == "AT1", `:=`(
  UMAP_1 = UMAP_1 + 0.7,
  UMAP_2 = UMAP_2 + 0.6
)]

label_df[Celltype == "AT2", `:=`(
  UMAP_1 = UMAP_1 + 0.2,
  UMAP_2 = UMAP_2 + 0.5
)]

label_df[Celltype == "Tumor", `:=`(
  UMAP_1 = UMAP_1 - 0.6,
  UMAP_2 = UMAP_2 + 0.7
)]

label_df[Celltype == "Fibroblast", `:=`(
  UMAP_1 = UMAP_1 - 0.8,
  UMAP_2 = UMAP_2 + 0.7
)]

label_df[Celltype == "Endothelial", `:=`(
  UMAP_1 = UMAP_1 - 0.7,
  UMAP_2 = UMAP_2 + 0.5
)]

label_df[Celltype == "Myeloid", `:=`(
  UMAP_1 = UMAP_1 - 0.5,
  UMAP_2 = UMAP_2 - 0.5
)]

label_df[Celltype == "Lymphoid", `:=`(
  UMAP_1 = UMAP_1 + 0.6,
  UMAP_2 = UMAP_2 - 0.5
)]

# ============================================================
# 8. Plot
# ============================================================

p_wnn <- ggplot(
  umap_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = Celltype
  )
) +
  geom_point(
    size = 0.7,
    alpha = 0.85
  ) +
  
  geom_label(
    data = label_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = Celltype,
      color = Celltype
    ),
    inherit.aes = FALSE,
    fill = "white",
    size = 3.0,
    fontface = "bold",
    label.size = 0.35,
    label.padding = unit(0.08, "lines"),
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = cell_cols,
    breaks = cell_levels,
    labels = legend_labels,
    drop = FALSE
  ) +
  
  coord_equal(
    xlim = c(
      xr[1] - 0.10 * dx,
      xr[2] + 0.08 * dx
    ),
    ylim = c(
      yr[1] - 0.08 * dy,
      yr[2] + 0.08 * dy
    ),
    clip = "off"
  ) +
  
  labs(
    color = "Cell identity"
  ) +
  
  theme_void() +
  
  theme(
    legend.position = "right",
    legend.justification = "center",
    
    legend.title = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    
    legend.text = element_text(
      size = 11,
      color = "black"
    ),
    
    legend.key.height = unit(0.38, "cm"),
    
    plot.margin = margin(
      8, 10, 8, 10
    )
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4.5,
        alpha = 1
      )
    )
  )

p_wnn


#save(lung_multi,file = 'E:/Lung/lung_muit.Rdata')

ggsave(
  filename = "E:/Lung/lung_WNN_celltype_UMAP.pdf",
  plot = p_wnn,
  width = 10,
  height = 6.5,
  units = "in",
  device = cairo_pdf
)





library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# 1. Prepare the data
# ============================================================

plot_df <- weight_summary %>%
  select(
    celltype,
    RNA_weight_mean,
    ATAC_weight_mean
  ) %>%
  pivot_longer(
    cols = c(RNA_weight_mean, ATAC_weight_mean),
    names_to = "Modality",
    values_to = "Weight"
  ) %>%
  mutate(
    Modality = recode(
      Modality,
      RNA_weight_mean  = "RNA",
      ATAC_weight_mean = "ATAC"
    )
  )

# ============================================================
# 2. Sort
# ============================================================

cell_order <- weight_summary %>%
  arrange(ATAC_weight_mean) %>%
  pull(celltype)

plot_df$celltype <- factor(
  plot_df$celltype,
  levels = cell_order
)

# ============================================================
# 3. Manually set vertical offsets
# ============================================================

plot_df <- plot_df %>%
  mutate(
    y_base = as.numeric(celltype),
    y_plot = case_when(
      Modality == "RNA"  ~ y_base + 0.10,
      Modality == "ATAC" ~ y_base - 0.10
    )
  )

# ============================================================
# 4. Color palette
# ============================================================

modality_cols <- c(
  "RNA"  = "#2F67B1",
  "ATAC" = "#F05A28"
)

# ============================================================
# 5. Plot
# ============================================================

p_weight <- ggplot(plot_df) +
  
  geom_segment(
    aes(
      x = 0,
      xend = Weight,
      y = y_plot,
      yend = y_plot,
      color = Modality
    ),
    linewidth = 1.0,
    alpha = 0.9
  ) +
  
  geom_point(
    aes(
      x = Weight,
      y = y_plot,
      color = Modality
    ),
    size = 4.5
  ) +
  
  geom_text(
    aes(
      x = Weight,
      y = y_plot,
      label = sprintf("%.3f", Weight),
      color = Modality
    ),
    hjust = -0.35,
    size = 3.5,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = modality_cols,
    breaks = c("ATAC", "RNA")
  ) +
  
  scale_y_continuous(
    breaks = seq_along(cell_order),
    labels = cell_order,
    expand = expansion(add = c(0.5, 0.5))
  ) +
  
  scale_x_continuous(
    limits = c(0, 1.08),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.01))
  ) +
  
  labs(
    x = "Mean modality weight",
    y = NULL,
    color = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.y = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 10,
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 13,
      face = "bold",
      margin = margin(t = 8)
    ),
    
    axis.ticks.y = element_blank(),
    
    legend.position = "top",
    legend.text = element_text(
      size = 11,
      face = "bold"
    ),
    
    panel.grid.major.x = element_line(
      linewidth = 0.3,
      linetype = "dashed",
      color = "grey85"
    ),
    
    panel.grid.major.y = element_blank(),
    
    plot.margin = margin(
      8, 35, 8, 8
    )
  )

p_weight


ggsave(
  filename = "E:/Lung/WNN_modality_weight_lollipop.pdf",
  plot = p_weight,
  width = 7,
  height = 5,
  units = "in",
  device = cairo_pdf
)




library(Signac)
library(GenomeInfoDb)
library(GenomicRanges)

atac_assay <- lung_multi[["ATAC"]]

gr <- granges(atac_assay)

GenomeInfoDb::genome(gr) <- "hg38"

atac_assay@ranges <- gr

lung_multi[["ATAC"]] <- atac_assay
GenomeInfoDb::genome(
  granges(lung_multi[["ATAC"]])
)
Annotation(lung_multi[["ATAC"]]) <- annotations





atac_assay <- lung_multi[["ATAC"]]

atac_assay@annotation <- annotations

lung_multi[["ATAC"]] <- atac_assay
Annotation(lung_multi[["ATAC"]])
DefaultAssay(lung_multi) <- "ATAC"

gene.activities <- GeneActivity(
  lung_multi,
  assay = "ATAC"
)


#################################################################

lung_multi






