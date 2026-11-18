library(data.table)
library(ggplot2)
library(grid)

data_dir <- "E:/Lung/EGA/data"

final_results <- fread(
  file.path(
    data_dir,
    "12_final_results.tsv"
  )
)

head(final_results)
colnames(final_results)


final_results <- fread(
  file.path(data_dir, "12_final_results.tsv")
)

eRegulon_results <- fread(
  file.path(data_dir, "11_eRegulon_CTS_MGRS_MRS.tsv")
)

target_MGRS <- fread(
  file.path(data_dir, "10_TF_celltype_target_MGRS.tsv")
)

mutation_edges <- fread(
  file.path(data_dir, "09_TF_peak_gene_mutation_edges.tsv.gz")
)

peak_scores <- fread(
  file.path(data_dir, "08_peak_mutation_scores.tsv")
)

mutation_peak <- fread(
  file.path(data_dir, "07_mutation_peak.tsv.gz")
)




library(data.table)
library(ggplot2)
library(grid)

# ============================================================
# 1. Basic settings
# ============================================================

dt <- copy(final_results)

cell_order <- c(
  "AT1",
  "AT2",
  "Endothelial",
  "Fibroblast",
  "Lymphoid",
  "Myeloid",
  "Tumor"
)

top_n <- 5L

dt <- dt[
  Celltype %in% cell_order
]


# ============================================================
# 2. Top five eRegulons for each cell type
# ============================================================

top_tf <- dt[
  !is.na(MRS_raw)
][
  order(
    match(Celltype, cell_order),
    -MRS_raw
  )
][
  ,
  head(.SD, top_n),
  by = Celltype
]



# ============================================================
# 3. All significant eRegulons
# ============================================================

sig_tf <- dt[
  CTS_P < 0.05 &
    MGRS_P < 0.05 &
    MRS_P < 0.05,
  .(
    TF,
    Celltype
  )
]

sig_tf


selected_tf <- unique(
  c(
    top_tf$TF,
    sig_tf$TF
  )
)

length(selected_tf)
selected_tf


# ============================================================
# 4. Retain results for these TFs across all cell types
# ============================================================

heat_dt <- dt[
  TF %in% selected_tf &
    Celltype %in% cell_order
]

full_grid <- CJ(
  TF = selected_tf,
  Celltype = cell_order,
  unique = TRUE
)

heat_dt <- merge(
  full_grid,
  heat_dt,
  by = c("TF", "Celltype"),
  all.x = TRUE
)




# ============================================================
# 5. Standardize MRS within each cell type
# ============================================================

heat_dt[
  ,
  Scaled_MRS := {
    x <- MRS_raw
    
    if (
      sum(!is.na(x)) >= 2 &&
      is.finite(sd(x, na.rm = TRUE)) &&
      sd(x, na.rm = TRUE) > 0
    ) {
      as.numeric(scale(x))
    } else {
      rep(NA_real_, .N)
    }
  },
  by = Celltype
]

heat_dt[
  Scaled_MRS > 2,
  Scaled_MRS := 2
]

heat_dt[
  Scaled_MRS < -2,
  Scaled_MRS := -2
]


# ============================================================
# 6. Significance labels
# ============================================================

heat_dt[, star := ""]

# Significant by both criteria
heat_dt[
  CTS_P < 0.05 &
    MGRS_P < 0.05,
  star := "*"
]

# Significant by all three criteria
heat_dt[
  CTS_P < 0.05 &
    MGRS_P < 0.05 &
    MRS_P < 0.05,
  star := "**"
]



heat_dt[
  star != "",
  .(
    TF,
    Celltype,
    MRS_raw,
    Scaled_MRS,
    CTS_P,
    MGRS_P,
    MRS_P,
    star
  )
][
  order(Celltype, -MRS_raw)
]


# ============================================================
# 7. TF order
# ============================================================

tf_anchor <- rbindlist(
  list(
    top_tf[
      ,
      .(
        TF,
        Celltype,
        MRS_raw
      )
    ],
    
    dt[
      CTS_P < 0.05 &
        MGRS_P < 0.05 &
        MRS_P < 0.05,
      .(
        TF,
        Celltype,
        MRS_raw
      )
    ]
  )
)

tf_anchor <- tf_anchor[
  order(
    match(Celltype, cell_order),
    -MRS_raw
  )
]

tf_order <- unique(
  tf_anchor$TF
)

heat_dt[
  ,
  TF := factor(
    TF,
    levels = tf_order
  )
]

heat_dt[
  ,
  Celltype := factor(
    Celltype,
    levels = rev(cell_order)
  )
]


# ============================================================
# 8. Plot
# ============================================================

p_top_regulon <- ggplot(
  heat_dt,
  aes(
    x = TF,
    y = Celltype,
    fill = Scaled_MRS
  )
) +
  
  geom_tile(
    color = "#606060",
    linewidth = 0.45,
    width = 0.98,
    height = 0.98
  ) +
  
  geom_text(
    aes(label = star),
    color = "black",
    size = 3.1,
    fontface = "bold"
  ) +
  
  scale_fill_gradient2(
    low = "#5CA6C8",
    mid = "#F5F3EC",
    high = "#EB604D",
    midpoint = 0,
    limits = c(-2, 2),
    breaks = c(-2, -1, 0, 1, 2),
    na.value = "#8B8B8B",
    name = "Scaled MRS"
  ) +
  
  labs(
    title = "Top eRegulons in each cell type",
    x = NULL,
    y = NULL,
    caption = paste0(
      "** CTS P < 0.05, MGRS P < 0.05, and MRS P < 0.05\n",
      "* CTS P < 0.05 and MGRS P < 0.05"
    )
  ) +
  
  guides(
    fill = guide_colorbar(
      title.position = "left",
      title.hjust = 1,
      barwidth = unit(7, "cm"),
      barheight = unit(0.35, "cm")
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0,
      size = 15
    ),
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    
    legend.position = "bottom",
    
    plot.caption = element_text(
      hjust = 0,
      size = 9.5,
      lineheight = 1.4,
      margin = margin(t = 12)
    ),
    
    plot.margin = margin(
      15, 20, 15, 15
    )
  )

print(p_top_regulon)


ggsave(
  filename = file.path(
    data_dir,
    "Top_eRegulons_each_celltype_heatmap.pdf"
  ),
  plot = p_top_regulon,
  width = 12,
  height = 5.5
)

ggsave(
  filename = file.path(
    data_dir,
    "Top_eRegulons_each_celltype_heatmap.png"
  ),
  plot = p_top_regulon,
  width = 12,
  height = 5.5,
  dpi = 600
)
