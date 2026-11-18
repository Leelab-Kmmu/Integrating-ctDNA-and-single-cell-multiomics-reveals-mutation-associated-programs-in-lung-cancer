library(data.table)
library(ggplot2)
library(ggrepel)
output_dir='E:/Lung/EGA/data'
final_results <- fread(
  file.path(output_dir, "12_final_results.tsv")
)


dt <- as.data.table(copy(final_results))

dt[, CTS_logP  := -log10(pmax(CTS_P, 1e-300))]
dt[, MGRS_logP := -log10(pmax(MGRS_P, 1e-300))]
dt[, MRS_logP  := -log10(pmax(MRS_P, 1e-300))]

p_cut <- -log10(0.05)

cell_levels <- c("AT1", "AT2", "Endothelial", "Fibroblast", "Lymphoid", "Myeloid", "Tumor")
dt[, Celltype := factor(Celltype, levels = cell_levels)]

n_per_type <- 4

label_dt <- dt[
  , {
    s1 <- .SD[StrongCandidate == TRUE][order(-MRS_null_z)]
    
    s2 <- .SD[
      Candidate == TRUE &
        !(TF %in% s1$TF)
    ][order(-MRS_null_z)]
    
    s3 <- .SD[
      MRS_P < 0.05 &
        MRS_null_z > 0 &
        !(TF %in% c(s1$TF, s2$TF))
    ][order(-MRS_null_z)]
    
    s4 <- .SD[
      !(TF %in% c(s1$TF, s2$TF, s3$TF))
    ][order(-MRS_null_z)]
    
    pool <- rbindlist(
      list(s1, s2, s3, s4),
      fill = TRUE
    )
    
    pool <- unique(pool, by = "TF")
    
    pool[
      1:min(n_per_type, .N),
      .(
        TF,
        CTS_logP,
        MGRS_logP,
        MRS_null_z
      )
    ]
  },
  by = Celltype
]


highlight_dt <- label_dt

cols <- c(
  "AT1" = "#F8766D",
  "AT2" = "#C49A00",
  "Endothelial" = "#53B400",
  "Fibroblast" = "#00C094",
  "Lymphoid" = "#00B6EB",
  "Myeloid" = "#A58AFF",
  "Tumor" = "#FB61D7"
)

p_b <- ggplot(
  dt,
  aes(x = CTS_logP, y = MGRS_logP)
) +
  # Background points
  geom_point(
    color = "grey82",
    size = 1.5,
    alpha = 0.55
  ) +
  
  # Threshold lines
  geom_vline(
    xintercept = p_cut,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey45"
  ) +
  geom_hline(
    yintercept = p_cut,
    linetype = "dashed",
    linewidth = 0.55,
    color = "grey45"
  ) +
  
  # All candidates
  geom_point(
    data = dt[Candidate == TRUE],
    aes(color = Celltype),
    size = 2.3,
    alpha = 0.85
  ) +
  
  # Include several large black-outlined points for each cell type
  geom_point(
    data = highlight_dt,
    aes(x = CTS_logP, y = MGRS_logP, fill = Celltype),
    inherit.aes = FALSE,
    shape = 21,
    color = "black",
    stroke = 0.8,
    size = 3.8
  ) +
  
  # Labels
  geom_text_repel(
    data = label_dt,
    aes(
      x = CTS_logP,
      y = MGRS_logP,
      label = TF,
      color = Celltype
    ),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "italic",
    box.padding = 0.4,
    point.padding = 0.35,
    min.segment.length = 0,
    segment.color = "grey45",
    segment.linewidth = 0.4,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = cols,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = cols,
    guide = "none",
    drop = FALSE
  ) +
  
  labs(
    x = expression(CTS~~(-log[10](P))),
    y = expression(MGRS~~(-log[10](P))),
    color = "Cell type"
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

p_b



ggsave(
  filename = "/400T/ckn/Lung_Model/PDF/Figure2A.pdf",
  plot = p_b,
  width = 8,
  height = 6,
  units = "in",
  device = cairo_pdf
)



##########################################################################################3
######Figure2B

library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ============================================================
# 1. Prepare data
# ============================================================

dt <- as.data.table(copy(final_results))

dt[, CTS_logP  := -log10(pmax(CTS_P, 1e-300))]
dt[, MGRS_logP := -log10(pmax(MGRS_P, 1e-300))]

p_cut <- -log10(0.05)

cell_levels <- c(
  "AT1",
  "AT2",
  "Endothelial",
  "Fibroblast",
  "Lymphoid",
  "Myeloid",
  "Tumor"
)

dt[, Celltype := factor(Celltype, levels = cell_levels)]

# ============================================================
# 2. Strict right-upper-corner definition
# Highlight points in red only when both CTS and MGRS are significant
# ============================================================

dt[, Strict_dual :=
     CTS_P < 0.05 &
     MGRS_P < 0.05
]

# Inspect the number of upper-right points for each cell type
strict_summary <- dt[
  ,
  .(
    N_total = .N,
    N_right_top = sum(Strict_dual, na.rm = TRUE)
  ),
  by = Celltype
]

print(strict_summary)

# To inspect specific TFs
print(
  dt[
    Strict_dual == TRUE,
    .(
      TF,
      Celltype,
      CTS_P,
      MGRS_P,
      MRS_P,
      CTS_logP,
      MGRS_logP
    )
  ][order(Celltype, CTS_P, MGRS_P)]
)

# ============================================================
# 3. Fix the coordinate ranges across all panels
# This keeps dashed-line positions and axis scales consistent
# ============================================================

x_max <- max(dt$CTS_logP, na.rm = TRUE)
y_max <- max(dt$MGRS_logP, na.rm = TRUE)

# Add a small margin
x_lim <- c(0, ceiling(x_max * 10) / 10 + 0.1)
y_lim <- c(0, ceiling(y_max * 10) / 10 + 0.1)

print(x_lim)
print(y_lim)

# ============================================================
# 4. Function for plotting a single panel
# ============================================================

plot_panel <- function(cell_name, dt, x_lim, y_lim, p_cut) {
  
  d <- dt[Celltype == cell_name]
  sig <- d[Strict_dual == TRUE]
  
  p <- ggplot(
    d,
    aes(
      x = CTS_logP,
      y = MGRS_logP
    )
  ) +
    
    # All background points: gray
    geom_point(
      color = "grey70",
      size = 1.6,
      alpha = 0.65
    ) +
    
    # Fixed threshold lines: identical across all plots
    geom_vline(
      xintercept = p_cut,
      linetype = "dashed",
      linewidth = 0.8,
      color = "#3659A7"
    ) +
    
    geom_hline(
      yintercept = p_cut,
      linetype = "dashed",
      linewidth = 0.8,
      color = "#3659A7"
    ) +
    
    # Only significant upper-right points: red
    geom_point(
      data = sig,
      aes(
        x = CTS_logP,
        y = MGRS_logP
      ),
      inherit.aes = FALSE,
      color = "#EF5350",
      size = 2.7,
      alpha = 1
    ) +
    
    # Label only significant upper-right points
    geom_text_repel(
      data = sig,
      aes(
        x = CTS_logP,
        y = MGRS_logP,
        label = TF
      ),
      inherit.aes = FALSE,
      color = "#D32F2F",
      size = 3.0,
      fontface = "italic",
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0,
      segment.color = "grey50",
      segment.linewidth = 0.35,
      max.overlaps = Inf
    ) +
    
    labs(
      title = as.character(cell_name),
      x = expression(CTS~(-log[10](P))),
      y = expression(MGRS~(-log[10](P)))
    ) +
    
    coord_cartesian(
      xlim = x_lim,
      ylim = y_lim,
      expand = FALSE
    ) +
    
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 11,
        face = "plain"
      ),
      axis.title = element_text(
        size = 10,
        color = "black"
      ),
      axis.text = element_text(
        size = 9,
        color = "black"
      ),
      panel.grid.major = element_line(
        color = "grey90",
        linewidth = 0.4
      ),
      panel.grid.minor = element_line(
        color = "grey95",
        linewidth = 0.3
      )
    )
  
  return(p)
}

# ============================================================
# 5. Generate panels for all seven cell types
# ============================================================

plot_list <- lapply(
  cell_levels,
  function(ct) plot_panel(
    cell_name = ct,
    dt = dt,
    x_lim = x_lim,
    y_lim = y_lim,
    p_cut = p_cut
  )
)

names(plot_list) <- cell_levels

# ============================================================
# 6. Combine panels
# Seven panels; a 4-column by 2-row layout is recommended
# ============================================================

p_all <- wrap_plots(
  plotlist = plot_list,
  ncol = 4
)

p_all
ggsave(
  filename = "/400T/ckn/Lung_Model/PDF/Figure2D_all_celltypes_strict_dual.pdf",
  plot = p_all,
  width = 14,
  height = 7.5,
  units = "in",
  device = cairo_pdf
)




###################################################################################3
#Figure2B
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)

dt <- as.data.table(copy(final_results))

# ------------------------------------------------------------
# 1. Select TFs to display
# ------------------------------------------------------------
# For example, display only strict-dual TFs
show_tf <- unique(
  dt[
    CTS_P < 0.05 &
      MGRS_P < 0.05,
    TF
  ]
)

# Use candidates instead to display more TFs
# show_tf <- unique(dt[Candidate == TRUE, TF])

plot_dt <- dt[
  TF %in% show_tf
]

# ------------------------------------------------------------
# 2. Construct the MRS z-score matrix
# Rows = TFs
# Columns = cell types
# ------------------------------------------------------------

mat_dt <- dcast(
  plot_dt,
  TF ~ Celltype,
  value.var = "MRS_null_z"
)

rownames_mat <- mat_dt$TF
mat <- as.matrix(
  mat_dt[
    ,
    -1,
    with = FALSE
  ]
)

rownames(mat) <- rownames_mat

# ------------------------------------------------------------
# 3. P-value matrix
# ------------------------------------------------------------

p_dt <- dcast(
  plot_dt,
  TF ~ Celltype,
  value.var = "MRS_P"
)

p_mat <- as.matrix(
  p_dt[
    ,
    -1,
    with = FALSE
  ]
)

rownames(p_mat) <- p_dt$TF

# Ensure the order matches mat
p_mat <- p_mat[
  rownames(mat),
  colnames(mat),
  drop = FALSE
]

# ------------------------------------------------------------
# 4. Asterisk matrix
# ------------------------------------------------------------

star_mat <- matrix(
  "",
  nrow = nrow(p_mat),
  ncol = ncol(p_mat),
  dimnames = dimnames(p_mat)
)

star_mat[
  !is.na(p_mat) &
    p_mat < 0.05
] <- "*"

star_mat[
  !is.na(p_mat) &
    p_mat < 0.01
] <- "**"

star_mat[
  !is.na(p_mat) &
    p_mat < 0.001
] <- "***"

# ------------------------------------------------------------
# 5. Representative cell type for each TF
#
# Defined here as the cell type with the highest MRS_null_z
# ------------------------------------------------------------

best_cell <- apply(
  mat,
  1,
  function(x) {
    
    if (all(is.na(x))) {
      return(NA_character_)
    }
    
    names(x)[
      which.max(x)
    ]
  }
)

# ------------------------------------------------------------
# 6. Colors
# ------------------------------------------------------------

cell_cols <- c(
  "AT1" = "#F8766D",
  "AT2" = "#C49A00",
  "Endothelial" = "#53B400",
  "Fibroblast" = "#00C094",
  "Lymphoid" = "#00B6EB",
  "Myeloid" = "#A58AFF",
  "Tumor" = "#FB61D7"
)

# MRS colors
col_fun <- colorRamp2(
  c(
    -2,
    0,
    2
  ),
  c(
    "#3B6AA0",
    "#F7F7F7",
    "#D73027"
  )
)

# ------------------------------------------------------------
# 7. Row annotation: dots on the left
# ------------------------------------------------------------

row_anno <- rowAnnotation(
  Celltype = anno_simple(
    best_cell,
    col = cell_cols,
    pch = 16,
    pt_size = unit(2.3, "mm"),
    simple_anno_size = unit(4, "mm")
  ),
  
  annotation_name_gp = gpar(
    fontsize = 9
  ),
  
  show_annotation_name = FALSE
)

# ------------------------------------------------------------
# 8. Heatmap
# ------------------------------------------------------------

ht <- Heatmap(
  mat,
  
  name = "Scaled MRS",
  
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_side = "right",
  
  row_names_gp = gpar(
    fontsize = 9
  ),
  
  column_names_gp = gpar(
    fontsize = 9
  ),
  
  column_names_rot = 45,
  
  rect_gp = gpar(
    col = "grey70",
    lwd = 0.5
  ),
  
  na_col = "white",
  
  # ----------------------------------------------------------
  # Asterisks
  # ----------------------------------------------------------
  cell_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
  ) {
    
    if (
      !is.na(star_mat[i, j]) &&
      star_mat[i, j] != ""
    ) {
      
      grid.text(
        star_mat[i, j],
        x,
        y,
        gp = gpar(
          fontsize = 8,
          fontface = "bold"
        )
      )
    }
  },
  
  heatmap_legend_param = list(
    title = "Scaled MRS",
    title_gp = gpar(
      fontsize = 9
    ),
    labels_gp = gpar(
      fontsize = 8
    )
  ),
  
  width = unit(
    5,
    "cm"
  )
)

# ------------------------------------------------------------
# 9. Draw
# ------------------------------------------------------------

ht_all <- row_anno + ht

draw(
  ht_all,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)



##############################################################################
###Figure2C
library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)

library(Seurat)
library(ggplot2)
library(patchwork)
library(data.table)

genes4 <- c(
  "TP63",
  "SOX18",
  "SPI1",
  "ZIC5"
)

expr <- SeuratObject::GetAssayData(
  lung_RNA,
  assay = "RNA",
  layer = "data"
)

umap_df <- as.data.frame(
  Embeddings(
    lung_RNA,
    reduction = "umap"
  )
)

colnames(umap_df)[1:2] <- c(
  "UMAP_1",
  "UMAP_2"
)

umap_df$Cell <- rownames(umap_df)

expr4 <- as.data.frame(
  t(
    as.matrix(
      expr[
        genes4,
        umap_df$Cell,
        drop = FALSE
      ]
    )
  )
)

expr4$Cell <- rownames(expr4)

plot_df <- merge(
  umap_df,
  expr4,
  by = "Cell",
  sort = FALSE
)

plot_df <- plot_df[
  match(
    umap_df$Cell,
    plot_df$Cell
  ),
]


# ============================================================
# UMAP with more prominent expression signals
# ============================================================

plot_gene_umap_strong <- function(
    gene,
    df
) {
  
  d <- data.table::as.data.table(
    data.table::copy(df)
  )
  
  d[, Expression := get(gene)]
  
  # Positive cells
  d_pos <- d[
    is.finite(Expression) &
      Expression > 0
  ]
  
  # Cap high expression values
  if (nrow(d_pos) > 0) {
    
    upper <- as.numeric(
      quantile(
        d_pos$Expression,
        0.95,
        na.rm = TRUE
      )
    )
    
    if (!is.finite(upper) || upper <= 0) {
      upper <- max(
        d_pos$Expression,
        na.rm = TRUE
      )
    }
    
    d_pos[
      ,
      Expression_plot := pmin(
        Expression,
        upper
      )
    ]
    
    data.table::setorder(
      d_pos,
      Expression_plot
    )
  }
  
  p <- ggplot() +
    
    # ========================================================
  # Use all cells to define the UMAP outline
  # Use small but slightly darker points to clarify the overall structure
  # ========================================================
  geom_point(
    data = d,
    aes(
      x = UMAP_1,
      y = UMAP_2
    ),
    color = "grey78",
    size = 0.42,
    alpha = 0.9
  )
  
  # ==========================================================
  # Overlay expressing cells
  # ==========================================================
  
  if (nrow(d_pos) > 0) {
    
    p <- p +
      
      geom_point(
        data = d_pos,
        aes(
          x = UMAP_1,
          y = UMAP_2,
          color = Expression_plot
        ),
        size = 0.55,
        alpha = 1
      ) +
      
      scale_color_gradient(
        low = "#FDB9A9",
        high = "#D7301F",
        name = "Expression"
      )
  }
  
  p +
    
    coord_fixed() +
    
    labs(
      title = gene
    ) +
    
    theme_void() +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 13,
        face = "italic",
        color = "black"
      ),
      
      legend.title = element_text(
        size = 9
      ),
      
      legend.text = element_text(
        size = 8
      ),
      
      plot.margin = margin(
        1,
        1,
        1,
        1
      )
    )
}
# ============================================================
# Plot
# ============================================================

p1 <- plot_gene_umap_strong(
  "TP63",
  plot_df
)

p2 <- plot_gene_umap_strong(
  "SOX18",
  plot_df
)

p3 <- plot_gene_umap_strong(
  "SPI1",
  plot_df
)

p4 <- plot_gene_umap_strong(
  "ZIC5",
  plot_df
)

p_gene_umap <- (
  p1 +
    p2 +
    p3 +
    p4
) +
  plot_layout(
    ncol = 4
  )

p_gene_umap


ggsave(
  filename = "/400T/ckn/Lung_Model/PDF/Figure2_TF_expression_UMAP_outline.pdf",
  plot = p_gene_umap,
  width = 12,
  height = 3.2,
  units = "in",
  device = cairo_pdf
)




