library(data.table)
library(ggplot2)
library(ggrepel)

plot_dt <- copy(final_results)

plot_dt[, x := -log10(CTS_P)]
plot_dt[, y := -log10(MGRS_P)]

plot_dt[, sig := 
          CTS_P < 0.05 &
          MGRS_P < 0.05 &
          MRS_P < 0.05]

# Label several TFs with the strongest MRS in each cell type
plot_dt[, label := ""]
plot_dt[
  order(Celltype, -MRS_null_z),
  label := fifelse(seq_len(.N) <= 5, TF, ""),
  by = Celltype
]

ggplot(plot_dt, aes(x = x, y = y)) +
  geom_point(
    aes(shape = sig),
    size = 2.4,
    alpha = 0.75
  ) +
  geom_vline(
    xintercept = -log10(0.05),
    linetype = 2
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2
  ) +
  geom_text_repel(
    data = plot_dt[label != ""],
    aes(label = label),
    size = 3.5,
    max.overlaps = Inf
  ) +
  facet_wrap(~Celltype, scales = "free") +
  labs(
    x = expression(CTS~~-log[10](P)),
    y = expression(MGRS~~-log[10](P))
  ) +
  theme_classic(base_size = 13)









plot_dt[, sig_group := fifelse(
  CTS_P < 0.05 &
    MGRS_P < 0.05 &
    MRS_P < 0.05,
  "Significant",
  "Not significant"
)]

ggplot(plot_dt, aes(x = x, y = y)) +
  geom_point(
    aes(color = sig_group),
    size = 2.4,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = -log10(0.05),
    linetype = 2
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = 2
  ) +
  geom_text_repel(
    data = plot_dt[
      sig_group == "Significant"
    ][order(-MRS_null_z)][1:min(.N, 20)],
    aes(label = TF),
    size = 3.2
  ) +
  facet_wrap(~Celltype) +
  labs(
    x = expression(CTS~~(-log[10](P))),
    y = expression(MGRS~~(-log[10](P))),
    color = NULL
  ) +
  theme_classic(base_size = 13)













library(data.table)
library(ggplot2)

dt <- copy(final_results)

#-----------------------------
# 1. Select top TFs
#-----------------------------
top_n <- 5

show_tfs <- unique(
  dt[order(Celltype, -MRS_null_z),
     head(.SD, top_n),
     by = Celltype]$TF
)

plot_dt <- dt[TF %in% show_tfs]

#-----------------------------
# 2. Define the order
#-----------------------------
celltype_order <- c("AT1", "AT2", "Endothelial", "Fibroblast", 
                    "Lymphoid", "Myeloid", "Tumor")
celltype_order <- intersect(celltype_order, unique(plot_dt$Celltype))

# Order TFs by appearance or by maximum MRS
tf_order <- plot_dt[
  , .(max_mrs = max(MRS_null_z, na.rm = TRUE)),
  by = TF
][order(-max_mrs)]$TF

plot_dt[, Celltype := factor(Celltype, levels = rev(celltype_order))]
plot_dt[, TF := factor(TF, levels = tf_order)]

#-----------------------------
# 3. Truncate the color range
#-----------------------------
plot_dt[, MRS_plot := pmax(pmin(MRS_null_z, 2), -2)]

#-----------------------------
# 4. Significance labels
# Adjust the threshold as needed
#-----------------------------
plot_dt[, sig_label := ""]
plot_dt[
  CTS_P < 0.05 & MGRS_P < 0.05 & MRS_P < 0.05 & MRS_null_z > 3,
  sig_label := "**"
]
plot_dt[
  CTS_P < 0.05 & MGRS_P < 0.05 & MRS_P < 0.05 &
    MRS_null_z > 1.5 & MRS_null_z <= 3,
  sig_label := "*"
]

#-----------------------------
# 5. Plot
#-----------------------------
p <- ggplot(plot_dt, aes(x = TF, y = Celltype, fill = MRS_plot)) +
  geom_tile(color = "#6E6E6E", linewidth = 0.55) +
  geom_text(aes(label = sig_label), size = 3.2, fontface = "bold") +
  scale_fill_gradientn(
    colours = c("#9BB7C9", "#EDEAE3", "#E47C72"),
    values = scales::rescale(c(-2, 0, 2)),
    limits = c(-2, 2),
    name = "Scaled MRS"
  ) +
  labs(title = "Top eRegulons in each cell type", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1,
      size = 8, color = "black"
    ),
    axis.text.y = element_text(size = 9, color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(8, 8, 28, 8)
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barwidth = unit(5, "cm"),
      barheight = unit(0.35, "cm")
    )
  ) +
  coord_fixed(ratio = 0.6) +
  annotate(
    "text",
    x = 1, y = -0.95,
    label = "** CTS P < 0.05, MGRS P < 0.05, MRS P < 0.05, and MRS > 3",
    hjust = 0, size = 3
  ) +
  annotate(
    "text",
    x = 1, y = -1.35,
    label = "* CTS P < 0.05, MGRS P < 0.05, MRS P < 0.05, and MRS > 1.5",
    hjust = 0, size = 3
  ) +
  coord_cartesian(clip = "off")

p


out_dir <- "E:/Lung/EGA/PDF/gene_list"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


ggsave(
  filename = file.path(out_dir, "Top_eRegulons_heatmap.pdf"),
  plot = p,
  width = 9.5,
  height = 3.8,
  bg = "white"
)



##################################################
###############################
library(data.table)

data_dir <- "E:/Lung/EGA/data"

edge <- fread(
  file.path(
    data_dir,
    "09_TF_peak_gene_mutation_edges.tsv.gz"
  )
)

dim(edge)
colnames(edge)
head(edge)


scale_to <- function(x, to = c(1, 5)) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(mean(to), length(x)))
  ux <- unique(na.omit(x))
  if (length(ux) <= 1) return(rep(mean(to), length(x)))
  scales::rescale(x, to = to, from = range(x, na.rm = TRUE))
}

plot_tf_peak_gene_network <- function(
    edge,
    tf_name,
    celltype_name,
    top_peaks = 10,
    top_genes_per_peak = 2,
    r_peak = 1.05,
    r_gene = 1.90
) {
  
  dt <- copy(edge)[
    TF == tf_name &
      Celltype == celltype_name
  ]
  
  if (nrow(dt) == 0) {
    stop(
      paste0(
        "没有找到 ",
        tf_name,
        " - ",
        celltype_name,
        " 对应的数据。"
      )
    )
  }
  
  #-------------------------------------------
  # 1. Summarize at the peak level
  #-------------------------------------------
  peak_info <- dt[, .(
    peak_edge = max(MGRS_edge, na.rm = TRUE),
    peak_mut  = max(MutationRelevance, na.rm = TRUE),
    n_genes   = uniqueN(gene)
  ), by = peak_id]
  
  # Prioritize peaks with high MGRS_edge values
  peak_info <- peak_info[
    order(-peak_edge, -peak_mut)
  ]
  
  peak_info <- peak_info[
    1:min(.N, top_peaks)
  ]
  
  dt <- dt[
    peak_id %in% peak_info$peak_id
  ]
  
  #-------------------------------------------
  # 2. Select the strongest target genes for each peak
  #-------------------------------------------
  dt <- dt[
    order(
      peak_id,
      -PeakGeneStrength,
      -MGRS_edge
    )
  ]
  
  # First remove duplicate peak-gene pairs
  dt <- dt[
    ,
    .SD[which.max(PeakGeneStrength)],
    by = .(peak_id, gene)
  ]
  
  # Select the top genes for each peak
  dt <- dt[
    ,
    head(.SD, top_genes_per_peak),
    by = peak_id
  ]
  
  #-------------------------------------------
  # 3. Recalculate information for the retained peaks
  #-------------------------------------------
  peak_info <- dt[, .(
    peak_edge = max(MGRS_edge, na.rm = TRUE),
    peak_mut  = max(MutationRelevance, na.rm = TRUE),
    n_genes   = uniqueN(gene)
  ), by = peak_id]
  
  peak_info <- peak_info[
    order(-peak_edge, -peak_mut)
  ]
  
  n_peak <- nrow(peak_info)
  
  if (n_peak == 0) {
    stop("筛选后没有可绘制的 peak。")
  }
  
  #-------------------------------------------
  # 4. Circular layout of peaks
  #-------------------------------------------
  peak_info[
    ,
    angle :=
      pi / 2 +
      2 * pi * (0:(n_peak - 1)) / n_peak
  ]
  
  peak_info[
    ,
    x := r_peak * cos(angle)
  ]
  
  peak_info[
    ,
    y := r_peak * sin(angle)
  ]
  
  # TF-to-peak edge width
  peak_info[
    ,
    lw_tfpeak := scale_to(
      peak_edge,
      c(1.2, 7)
    )
  ]
  
  # peak label
  peak_info[
    ,
    peak_label_x :=
      (r_peak + 0.11) * cos(angle)
  ]
  
  peak_info[
    ,
    peak_label_y :=
      (r_peak + 0.11) * sin(angle)
  ]
  
  peak_info[
    ,
    peak_hjust := fifelse(
      peak_label_x >= 0,
      0,
      1
    )
  ]
  
  #-------------------------------------------
  # 5. Gene node layout
  #-------------------------------------------
  gene_dt <- merge(
    dt,
    peak_info[, .(
      peak_id,
      peak_angle = angle,
      peak_x = x,
      peak_y = y
    )],
    by = "peak_id",
    all.x = TRUE
  )
  
  # Spread genes associated with each peak across a small angular range
  gene_dt <- gene_dt[, {
    
    n <- .N
    
    sp <- min(
      0.36,
      0.13 * (n - 1)
    )
    
    off <- if (n == 1) {
      0
    } else {
      seq(
        -sp,
        sp,
        length.out = n
      )
    }
    
    cbind(.SD, offset = off)
    
  }, by = peak_id]
  
  gene_dt[
    ,
    gene_angle :=
      peak_angle + offset
  ]
  
  gene_dt[
    ,
    gene_x :=
      r_gene * cos(gene_angle)
  ]
  
  gene_dt[
    ,
    gene_y :=
      r_gene * sin(gene_angle)
  ]
  
  # Move gene labels slightly farther outward
  gene_dt[
    ,
    gene_label_x :=
      (r_gene + 0.10) * cos(gene_angle)
  ]
  
  gene_dt[
    ,
    gene_label_y :=
      (r_gene + 0.10) * sin(gene_angle)
  ]
  
  gene_dt[
    ,
    gene_hjust := fifelse(
      gene_label_x >= 0,
      0,
      1
    )
  ]
  
  # Peak-to-gene edge width
  gene_dt[
    ,
    lw_peakgene := scale_to(
      PeakGeneStrength,
      c(0.6, 3.5)
    )
  ]
  
  #-------------------------------------------
  # 6. TF-to-peak edges
  #-------------------------------------------
  edge_tf_peak <- peak_info[, .(
    x = 0,
    y = 0,
    xend = x,
    yend = y,
    linewidth = lw_tfpeak
  )]
  
  #-------------------------------------------
  # 7. Peak-to-gene edges
  #-------------------------------------------
  edge_peak_gene <- gene_dt[, .(
    x = peak_x,
    y = peak_y,
    xend = gene_x,
    yend = gene_y,
    linewidth = lw_peakgene
  )]
  
  #-------------------------------------------
  # 8. Plot
  #-------------------------------------------
  p <- ggplot() +
    
    # TF -> peak
    geom_segment(
      data = edge_tf_peak,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        linewidth = linewidth
      ),
      color = "#C2C4C6",
      alpha = 0.90,
      lineend = "round",
      show.legend = FALSE
    ) +
    
    # peak -> gene
    geom_segment(
      data = edge_peak_gene,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        linewidth = linewidth
      ),
      color = "#BBB0CF",
      alpha = 0.80,
      lineend = "round",
      show.legend = FALSE
    ) +
    
    # peak nodes
    geom_point(
      data = peak_info,
      aes(
        x = x,
        y = y,
        size = peak_edge,
        fill = peak_mut
      ),
      shape = 21,
      color = "#367A9C",
      stroke = 1.5
    ) +
    
    # gene nodes
    geom_point(
      data = gene_dt,
      aes(
        x = gene_x,
        y = gene_y
      ),
      shape = 21,
      size = 3.1,
      fill = "#F4F4F2",
      color = "#367A9C",
      stroke = 1.25
    ) +
    
    # Central TF
    geom_point(
      aes(
        x = 0,
        y = 0
      ),
      shape = 23,
      size = 8,
      fill = "#D87873",
      color = "#8F3030",
      stroke = 1.7
    ) +
    
    # TF label
    annotate(
      "text",
      x = 0,
      y = -0.14,
      label = tf_name,
      size = 5.0,
      fontface = "bold",
      color = "#111111"
    ) +
    
    # peak labels
    geom_text(
      data = peak_info,
      aes(
        x = peak_label_x,
        y = peak_label_y,
        label = peak_id,
        hjust = peak_hjust
      ),
      size = 2.5,
      color = "#326C8C"
    ) +
    
    # gene labels
    geom_text(
      data = gene_dt,
      aes(
        x = gene_label_x,
        y = gene_label_y,
        label = gene,
        hjust = gene_hjust
      ),
      size = 3.0,
      color = "#363636"
    ) +
    
    # Mutation relevance
    scale_fill_gradientn(
      colours = c(
        "#B8CDD5",
        "#EEE8E1",
        "#D98A82"
      ),
      name = "Mutation\nrelevance"
    ) +
    
    # Peak size = MGRS edge
    scale_size_continuous(
      range = c(4.5, 9),
      name = "MGRS edge"
    ) +
    
    # Use line width internally without displaying a legend
    scale_linewidth_identity() +
    
    coord_equal(
      xlim = c(-2.45, 2.45),
      ylim = c(-2.25, 2.25),
      clip = "off"
    ) +
    
    labs(
      title = paste0(
        tf_name,
        " regulon in ",
        celltype_name
      )
    ) +
    
    theme_void(base_size = 11) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14,
        margin = margin(b = 10)
      ),
      
      legend.position = "right",
      
      legend.title = element_text(
        size = 10,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 9
      ),
      
      plot.margin = margin(
        20,
        100,
        20,
        100
      )
    )
  
  return(p)
}



p_net <- plot_tf_peak_gene_network(
  edge = edge,
  tf_name = "E2F3",
  celltype_name = "Tumor",
  top_peaks = 10,
  top_genes_per_peak = 2
)

p_net





##################################

library(data.table)
library(ggplot2)
library(scales)
library(ggrepel)

scale_to <- function(x, to = c(1, 5)) {
  
  x <- as.numeric(x)
  
  if (all(is.na(x))) {
    return(rep(mean(to), length(x)))
  }
  
  r <- range(x, na.rm = TRUE)
  
  if (diff(r) == 0) {
    return(rep(mean(to), length(x)))
  }
  
  scales::rescale(x, to = to, from = r)
}

plot_tf_peak_gene_network2 <- function(
    edge,
    tf_name,
    celltype_name,
    top_peaks = 8,
    top_genes_per_peak = 2,
    r_peak = 1.00,
    r_gene = 1.85
) {
  
  dt <- copy(edge)[
    TF == tf_name &
      Celltype == celltype_name &
      !is.na(peak_id) &
      !is.na(gene)
  ]
  
  if (nrow(dt) == 0) {
    stop("没有找到对应 TF × Celltype 数据")
  }
  
  #====================================================
  # 1. Remove duplicate peak-gene pairs
  #====================================================
  
  dt <- dt[
    order(-MGRS_edge, -PeakGeneStrength)
  ]
  
  dt <- dt[
    ,
    .SD[1],
    by = .(peak_id, gene)
  ]
  
  
  #====================================================
  # 2. Peak ranking
  #====================================================
  
  peak_info <- dt[, .(
    
    # Combined mutation-linked regulatory contribution
    peak_MGRS =
      max(MGRS_edge, na.rm = TRUE),
    
    # TF/peak cell-type regulatory evidence
    peak_CellEdge =
      max(CellEdgeScore, na.rm = TRUE),
    
    # mutation relevance
    MutationRelevance =
      max(MutationRelevance, na.rm = TRUE),
    
    n_genes =
      uniqueN(gene)
    
  ), by = peak_id]
  
  
  peak_info <- peak_info[
    order(
      -peak_MGRS,
      -MutationRelevance
    )
  ]
  
  peak_info <- head(
    peak_info,
    top_peaks
  )
  
  dt <- dt[
    peak_id %in% peak_info$peak_id
  ]
  
  
  #====================================================
  # 3. Select top target genes for each peak
  #====================================================
  
  dt <- dt[
    order(
      peak_id,
      -PeakGeneStrength,
      -MGRS_edge
    )
  ]
  
  dt <- dt[
    ,
    head(.SD, top_genes_per_peak),
    by = peak_id
  ]
  
  
  #====================================================
  # 4. Update peak information
  #====================================================
  
  peak_info <- dt[, .(
    
    peak_MGRS =
      max(MGRS_edge, na.rm = TRUE),
    
    peak_CellEdge =
      max(CellEdgeScore, na.rm = TRUE),
    
    MutationRelevance =
      max(MutationRelevance, na.rm = TRUE)
    
  ), by = peak_id]
  
  peak_info <- peak_info[
    order(-peak_MGRS)
  ]
  
  n_peak <- nrow(peak_info)
  
  
  #====================================================
  # 5. Peak radial coordinates
  #====================================================
  
  # Arrange clockwise starting from the top
  peak_info[, angle :=
              pi/2 -
              2*pi*(0:(n_peak - 1))/n_peak
  ]
  
  peak_info[, peak_x :=
              r_peak * cos(angle)
  ]
  
  peak_info[, peak_y :=
              r_peak * sin(angle)
  ]
  
  
  # TF-to-peak edge width
  peak_info[, tf_peak_lwd :=
              scale_to(
                peak_CellEdge,
                c(0.8, 4.0)
              )
  ]
  
  
  # Peak size
  peak_info[, peak_size :=
              scale_to(
                peak_MGRS,
                c(4.5, 9)
              )
  ]
  
  
  # Place peak labels slightly outside the peaks
  peak_info[, peak_label_x :=
              (r_peak + 0.13) * cos(angle)
  ]
  
  peak_info[, peak_label_y :=
              (r_peak + 0.13) * sin(angle)
  ]
  
  
  #====================================================
  # 6. Gene coordinates
  #====================================================
  
  gene_dt <- merge(
    dt,
    peak_info[, .(
      peak_id,
      peak_angle = angle,
      peak_x,
      peak_y
    )],
    by = "peak_id"
  )
  
  
  # Spread genes from the same peak slightly to both sides
  gene_dt[, gene_offset := {
    
    n <- .N
    
    if (n == 1) {
      
      0
      
    } else {
      
      seq(
        -0.12,
        0.12,
        length.out = n
      )
    }
    
  }, by = peak_id]
  
  
  gene_dt[, gene_angle :=
            peak_angle + gene_offset
  ]
  
  gene_dt[, gene_x :=
            r_gene * cos(gene_angle)
  ]
  
  gene_dt[, gene_y :=
            r_gene * sin(gene_angle)
  ]
  
  
  # Move gene labels farther outward
  gene_dt[, gene_label_x :=
            (r_gene + 0.15) *
            cos(gene_angle)
  ]
  
  gene_dt[, gene_label_y :=
            (r_gene + 0.15) *
            sin(gene_angle)
  ]
  
  
  #====================================================
  # 7. Edge width
  #====================================================
  
  gene_dt[, pg_lwd :=
            scale_to(
              PeakGeneStrength,
              c(0.5, 2.3)
            )
  ]
  
  
  #====================================================
  # 8. TF → Peak edges
  #====================================================
  
  tf_peak_edges <- peak_info[, .(
    x = 0,
    y = 0,
    xend = peak_x,
    yend = peak_y,
    lwd = tf_peak_lwd
  )]
  
  
  #====================================================
  # 9. Plot
  #====================================================
  
  p <- ggplot() +
    
    #-------------------------
  # TF -> peak
  #-------------------------
  
  geom_segment(
    data = tf_peak_edges,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      linewidth = lwd
    ),
    color = "#BFC2C4",
    alpha = 0.75,
    lineend = "round",
    show.legend = FALSE
  ) +
    
    
    #-------------------------
  # peak -> gene
  #-------------------------
  
  geom_curve(
    data = gene_dt,
    aes(
      x = peak_x,
      y = peak_y,
      xend = gene_x,
      yend = gene_y,
      linewidth = pg_lwd
    ),
    curvature = 0.08,
    color = "#B8A9C8",
    alpha = 0.72,
    lineend = "round",
    show.legend = FALSE
  ) +
    
    
    #-------------------------
  # Gene nodes
  #-------------------------
  
  geom_point(
    data = gene_dt,
    aes(
      x = gene_x,
      y = gene_y
    ),
    shape = 21,
    size = 3.3,
    stroke = 1.3,
    fill = "#F8F7F4",
    color = "#467C98"
  ) +
    
    
    #-------------------------
  # Peak nodes
  #-------------------------
  
  geom_point(
    data = peak_info,
    aes(
      x = peak_x,
      y = peak_y,
      size = peak_MGRS,
      fill = MutationRelevance
    ),
    shape = 21,
    color = "#467C98",
    stroke = 1.5
  ) +
    
    
    #-------------------------
  # Center TF
  #-------------------------
  
  geom_point(
    aes(
      x = 0,
      y = 0
    ),
    shape = 23,
    size = 8.5,
    stroke = 1.7,
    fill = "#D9827C",
    color = "#93413E"
  ) +
    
    
    # Place the TF name below the node to avoid overlap
    annotate(
      "text",
      x = 0,
      y = -0.17,
      label = tf_name,
      size = 5.2,
      fontface = "bold"
    ) +
    
    
    #-------------------------
  # Peak coordinates
  #-------------------------
  
  geom_text(
    data = peak_info,
    aes(
      x = peak_label_x,
      y = peak_label_y,
      label = peak_id
    ),
    size = 2.35,
    color = "#477692"
  ) +
    
    
    #-------------------------
  # Gene labels
  #-------------------------
  
  geom_text_repel(
    data = gene_dt,
    aes(
      x = gene_label_x,
      y = gene_label_y,
      label = gene
    ),
    size = 3.2,
    color = "#353535",
    
    box.padding = 0.25,
    point.padding = 0.15,
    
    min.segment.length = 0,
    segment.color = NA,
    
    max.overlaps = Inf,
    
    seed = 123
  ) +
    
    
    #==================================================
  # Color
  #==================================================
  
  scale_fill_gradientn(
    colours = c(
      "#AFC6D1",
      "#E7E4DE",
      "#D9958D"
    ),
    name = "Mutation\nrelevance"
  ) +
    
    
    # peak size
    scale_size_continuous(
      range = c(4.5, 8.5),
      name = "MGRS\nedge"
    ) +
    
    
    scale_linewidth_identity() +
    
    
    coord_equal(
      xlim = c(-2.4, 2.4),
      ylim = c(-2.25, 2.25),
      clip = "off"
    ) +
    
    
    labs(
      title = paste0(
        tf_name,
        " regulon · ",
        celltype_name
      )
    ) +
    
    
    theme_void() +
    
    theme(
      
      plot.title = element_text(
        size = 13,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 5)
      ),
      
      legend.position = "right",
      
      legend.title = element_text(
        size = 9,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 8
      ),
      
      plot.margin = margin(
        10,
        50,
        10,
        50
      )
    )
  
  p
}


p_net <- plot_tf_peak_gene_network2(
  edge = edge,
  tf_name = "E2F3",
  celltype_name = "Tumor",
  top_peaks = 8,
  top_genes_per_peak = 2
)

p_net



library(data.table)
library(ggplot2)

#====================================================
# 1. Output directory
#====================================================

x_cut <- 0
y_cut <- 1.3

upper_right_sig <- final_results[
  CTS_null_z > x_cut &
    MRS_null_z > y_cut &
    CTS_P < 0.05 &
    MGRS_P < 0.05 &
    MRS_P < 0.05
]

# For more stringent filtering:
# upper_right_sig <- upper_right_sig[StrongCandidate == TRUE]

upper_right_sig <- upper_right_sig[
  order(Celltype, -MRS_null_z)
]

upper_right_sig[, .(
  TF,
  Celltype,
  CTS_null_z,
  MRS_null_z,
  CTS_P,
  MGRS_P,
  MRS_P,
  Candidate,
  StrongCandidate
)]

# First inspect the number of entries:

nrow(upper_right_sig)
table(upper_right_sig$Celltype)


library(ggplot2)

for (i in seq_len(nrow(upper_right_sig))) {
  
  tf_i <- upper_right_sig$TF[i]
  ct_i <- upper_right_sig$Celltype[i]
  
  message(
    "[", i, "/", nrow(upper_right_sig), "] ",
    "绘制网络图: ", tf_i, " - ", ct_i
  )
  
  # Check whether edge data are available
  n_edge <- edge[
    TF == tf_i &
      Celltype == ct_i,
    .N
  ]
  
  if (n_edge == 0) {
    message("  跳过：edge 中没有对应数据")
    next
  }
  
  # Plot
  p_net <- tryCatch(
    plot_tf_peak_gene_network2(
      edge = edge,
      tf_name = tf_i,
      celltype_name = ct_i,
      top_peaks = 8,
      top_genes_per_peak = 2
    ),
    error = function(e) {
      message("  绘图失败：", e$message)
      return(NULL)
    }
  )
  
  if (is.null(p_net)) next
  
  # File name
  pdf_name <- paste0(
    sprintf("%02d", i),
    "_",
    ct_i,
    "_",
    tf_i,
    "_upper_right_network.pdf"
  )
  
  pdf_name <- gsub("[\\\\/:*?\"<>| ]", "_", pdf_name)
  
  # Save as PDF
  ggsave(
    filename = file.path(out_dir, pdf_name),
    plot = p_net,
    width = 9,
    height = 7.5,
    units = "in",
    bg = "white"
  )
}





library(data.table)

target_dir <- "E:/Lung/EGA/PDF/gene_list/target_gene_txt"

if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE)
}

for (i in seq_len(nrow(upper_right_sig))) {
  
  tf_i <- upper_right_sig$TF[i]
  ct_i <- upper_right_sig$Celltype[i]
  
  message(
    "[", i, "/", nrow(upper_right_sig), "] ",
    tf_i, " - ", ct_i
  )
  
  # Extract all target genes for this TF-by-cell-type pair
  genes <- edge[
    TF == tf_i &
      Celltype == ct_i &
      !is.na(gene) &
      gene != "",
    unique(gene)
  ]
  
  # Place the TF in the first row and sort the remaining target genes
  genes <- c(
    tf_i,
    sort(setdiff(unique(genes), tf_i))
  )
  
  safe_tf <- gsub("[\\\\/:*?\"<>| ]", "_", tf_i)
  safe_ct <- gsub("[\\\\/:*?\"<>| ]", "_", ct_i)
  
  file_name <- paste0(
    safe_ct,
    "_",
    safe_tf,
    "_target_genes_with_TF.txt"
  )
  
  fwrite(
    data.table(gene = genes),
    file.path(target_dir, file_name),
    col.names = FALSE,
    sep = "\t"
  )
}

message("完成：", target_dir)



################################################
subgrn <- fread(
  file.path(
    data_dir,
    "06_celltype_subGRN.tsv.gz"
  )
)

colnames(subgrn)
head(subgrn)



mc <- fread(
  file.path(
    data_dir,
    "12_matched_MonteCarlo_statistics.tsv"
  )
)

colnames(mc)
head(mc)



regulon_size <- subgrn[
  ,
  .(
    N_targets = uniqueN(gene),
    N_peaks   = uniqueN(peak_id)
  ),
  by = .(TF, Celltype)
]



regulon_size <- merge(
  regulon_size,
  final_results[, .(
    TF, Celltype,
    MRS_null_z,
    MRS_P,
    Candidate,
    StrongCandidate
  )],
  by = c("TF", "Celltype"),
  all.x = TRUE
)


ggplot(
  regulon_size[Candidate == TRUE],
  aes(
    x = Celltype,
    y = N_targets
  )
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  theme_classic()


e2f3_targets <- subgrn[
  TF == "E2F3" &
    Celltype == "Tumor",
  .(
    N_peaks = uniqueN(peak_id),
    Max_CellEdgeScore = max(CellEdgeScore, na.rm = TRUE),
    Max_PeakGeneStrength = max(PeakGeneStrength, na.rm = TRUE),
    Target_CTS = max(Target_CTS, na.rm = TRUE)
  ),
  by = gene
]

setorder(e2f3_targets, -Max_CellEdgeScore)

head(e2f3_targets, 20)

plot_dt <- head(e2f3_targets, 20)

ggplot(
  plot_dt,
  aes(
    x = reorder(gene, Max_CellEdgeScore),
    y = Max_CellEdgeScore
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "CellEdgeScore",
    title = "Top targets of E2F3 regulon in Tumor"
  ) +
  theme_classic()



mc_plot <- copy(mc)

mc_plot[, label := paste(TF, Celltype, sep = " - ")]

mc_plot <- mc_plot[
  order(-MRS_null_z)
]

mc_plot[, rank := .I]



ggplot(
  mc_plot,
  aes(
    x = rank,
    y = MRS_null_z
  )
) +
  geom_point(
    aes(color = Celltype),
    alpha = 0.65
  ) +
  geom_hline(
    yintercept = 1.96,
    linetype = 2
  ) +
  theme_classic() +
  labs(
    x = "Regulon rank",
    y = "MRS null z-score"
  )



mc2 <- merge(
  mc,
  final_results[, .(
    TF,
    Celltype,
    MRS_raw
  )],
  by = c("TF", "Celltype")
)


ggplot(
  mc2,
  aes(
    x = Null_MRS_mean,
    y = MRS_raw
  )
) +
  geom_point(
    aes(color = Celltype),
    alpha = 0.7
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  theme_classic() +
  labs(
    x = "Mean MRS of matched controls",
    y = "Observed MRS"
  )


heat_dt <- mc[
  CTS_P < 0.05 &
    MGRS_P < 0.05 &
    MRS_P < 0.05
]

heat_dt[, ID := paste(TF, Celltype, sep = " | ")]


mat <- as.matrix(
  heat_dt[, .(
    CTS_null_z,
    MGRS_null_z,
    MRS_null_z
  )]
)

rownames(mat) <- heat_dt$ID


library(ComplexHeatmap)
library(circlize)

Heatmap(
  mat,
  name = "z-score",
  col = colorRamp2(
    c(-3, 0, 3),
    c("#9CB5C6", "#F2EFEA", "#D58C82")
  ),
  cluster_columns = FALSE,
  cluster_rows = TRUE
)
