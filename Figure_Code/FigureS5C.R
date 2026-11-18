library(data.table)
library(ggplot2)

df_reg <- fread("E:/Lung/EGA/data/12_final_results.tsv")

plot_df <- df_reg[
  StrongCandidate == TRUE
]

plot_df[, Label := paste(TF, Celltype, sep = " - ")]
plot_df[, Label := factor(
  Label,
  levels = Label[order(MRS_null_z)]
)]




p_null <- ggplot(
  plot_df,
  aes(
    y = Label
  )
) +
  
  # null mean ± SD
  geom_errorbarh(
    aes(
      xmin = Null_MRS_mean - Null_MRS_sd,
      xmax = Null_MRS_mean + Null_MRS_sd
    ),
    height = 0.18,
    linewidth = 0.7,
    color = "grey50"
  ) +
  
  # null mean
  geom_point(
    aes(x = Null_MRS_mean),
    size = 2.5,
    color = "grey40"
  ) +
  
  # observed MRS
  geom_point(
    aes(x = MRS_raw),
    size = 3,
    color = "#D73027"
  ) +
  
  # observed-null mean connection
  geom_segment(
    aes(
      x = Null_MRS_mean,
      xend = MRS_raw,
      y = Label,
      yend = Label
    ),
    linewidth = 0.8
  ) +
  
  labs(
    x = "MRS",
    y = NULL
  ) +
  
  theme_classic(base_size = 12)

p_null


ggsave

ggsave(
  filename = "E:/Lung/EGA/PDF/MRS_observed_vs_null.pdf",
  plot = p_null,
  width = 7,
  height = 12,
  units = "in"
)


##################################################################################
score_df <- df_reg[
  Candidate == TRUE,
  .(
    TF,
    Celltype,
    CTS_z,
    MGRS_z,
    MRS_null_z
  )
]

score_long <- melt(
  score_df,
  id.vars = c("TF", "Celltype"),
  measure.vars = c(
    "CTS_z",
    "MGRS_z",
    "MRS_null_z"
  ),
  variable.name = "Score",
  value.name = "Z"
)


top_tf <- df_reg[
  Candidate == TRUE
][
  order(-MRS_null_z)
][1:min(30, .N), TF]

score_long <- score_long[
  TF %in% top_tf
]



ggplot(
  score_long,
  aes(
    x = Score,
    y = TF
  )
) +
  geom_point(
    aes(
      size = abs(Z),
      color = Z
    )
  ) +
  facet_grid(
    Celltype ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_color_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0
  ) +
  labs(
    x = NULL,
    y = NULL,
    color = "Z-score",
    size = "|Z|"
  ) +
  theme_classic(base_size = 11)

##################################################################
library(Seurat)
library(data.table)
library(ggplot2)

# ============================================================
# 1. Original top five for each cell type without removing duplicate TFs
# ============================================================

cell_order <- c(
  "AT1",
  "AT2",
  "Tumor",
  "Endothelial",
  "Fibroblast",
  "Lymphoid",
  "Myeloid"
)

top5_tf <- df_reg[
  is.finite(MRS_null_z),
  .SD[order(-MRS_null_z)][1:min(5, .N)],
  by = Celltype
]

top5_tf <- top5_tf[
  TF %in% rownames(lung_multi[["RNA"]])
]

top5_tf[, Celltype := factor(
  Celltype,
  levels = cell_order
)]

setorder(
  top5_tf,
  Celltype,
  -MRS_null_z
)

top5_tf[, .(
  TF,
  Celltype,
  MRS_null_z
)]


# ============================================================
# 2. Use only unique TFs in the dot plot
# ============================================================

tf_unique <- unique(top5_tf$TF)

DefaultAssay(lung_multi) <- "RNA"

dp_base <- DotPlot(
  lung_multi,
  features = tf_unique,
  assay = "RNA",
  group.by = "celltype"
)$data

dp_base <- as.data.table(dp_base)


# ============================================================
# 3. Replicate the same TF across different source cell types
# ============================================================

tf_group <- top5_tf[
  ,
  .(
    features.plot = TF,
    SourceCelltype = as.character(Celltype),
    MRS_null_z
  )
]

dp <- merge(
  dp_base,
  tf_group,
  by = "features.plot",
  allow.cartesian = TRUE
)




dp[, SourceCelltype := factor(
  SourceCelltype,
  levels = cell_order
)]

dp[, id := factor(
  id,
  levels = rev(cell_order)
)]


dp[, TF_ID := paste(
  SourceCelltype,
  features.plot,
  sep = "___"
)]


tf_order_dt <- top5_tf[
  ,
  .(
    TF_ID = paste(
      Celltype,
      TF,
      sep = "___"
    )
  )
]

tf_order <- tf_order_dt$TF_ID

dp[, TF_ID := factor(
  TF_ID,
  levels = tf_order
)]




p_tf_dot <- ggplot(
  dp,
  aes(
    x = TF_ID,
    y = id
  )
) +
  
  geom_point(
    aes(
      size = pct.exp,
      color = avg.exp.scaled
    )
  ) +
  
  facet_grid(
    . ~ SourceCelltype,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_x_discrete(
    labels = function(x) {
      sub("^.*___", "", x)
    }
  ) +
  
  scale_size_continuous(
    range = c(0.4, 5.5),
    breaks = c(0, 25, 50, 75),
    limits = c(0, 100),
    name = "Percent Expressed"
  ) +
  
  scale_color_gradient2(
    low = "#9EA4D6",
    mid = "#F7F7F7",
    high = "#C83B37",
    midpoint = 0,
    name = "Average Expression"
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    strip.background = element_rect(
      fill = "#263238",
      color = "#263238",
      linewidth = 0.6
    ),
    
    strip.text.x = element_text(
      color = "white",
      face = "bold",
      size = 9
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 9,
      color = "black",
      face = "bold"
    ),
    
    axis.ticks = element_blank(),
    
    panel.border = element_rect(
      color = "#444444",
      fill = NA,
      linewidth = 0.7
    ),
    
    panel.spacing.x = unit(
      0.08,
      "lines"
    ),
    
    legend.position = "right"
  )

p_tf_dot

ggsave(
  filename = "E:/Lung/EGA/PDF/Top5_TF_expression_dotplot.pdf",
  plot = p_tf_dot,
  width = 14,
  height = 5.5,
  units = "in"
)
