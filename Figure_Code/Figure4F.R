genes_E2F3 <- scan(
  "E:/Lung/EGA/PDF/gene_list/target_gene_txt/Tumor_E2F3_target_genes_with_TF.txt",
  what = character(),
  quiet = TRUE
)

genes_E2F3



lung_RNA=lung_multi
genes_E2F3 <- intersect(
  genes_E2F3,
  rownames(lung_RNA)
)

length(genes_E2F3)
genes_E2F3


lung_RNA <- AddModuleScore(
  lung_RNA,
  features = list(genes_E2F3),
  name = "E2F3_Regulon"
)

lung_RNA$E2F3_Regulon1

FeaturePlot(
  lung_RNA,
  features = "E2F3_Regulon1",
  reduction = "wnn.umap"
)

###################################################
genes_TBX15 <- scan(
  "E:/Lung/EGA/PDF/gene_list/target_gene_txt/AT2_TBX15_target_genes_with_TF.txt",
  what = character(),
  quiet = TRUE
)

genes_TBX15
DefaultAssay(lung_RNA) <- "RNA"

lung_RNA <- AddModuleScore(
  lung_RNA,
  features = list(genes_TBX15),
  name = "TBX15_Regulon"
)


FeaturePlot(
  lung_RNA,
  features = "TBX15_Regulon1",
  reduction = "umap"
)



genes_SPI1 <- scan(
  "E:/Lung/EGA/PDF/gene_list/target_gene_txt/Myeloid_SPI1_target_genes_with_TF.txt",
  what = character(),
  quiet = TRUE
)

genes_SPI1


DefaultAssay(lung_RNA) <- "RNA"

lung_RNA <- AddModuleScore(
  lung_RNA,
  features = list(genes_SPI1),
  name = "SPI1_Regulon"
)


FeaturePlot(
  lung_RNA,
  features = "SPI1_Regulon1",
  reduction = "umap"
)




FeaturePlot(lung_RNA, features = c("SPI1", "TBX15", "E2F3"))

VlnPlot(
  lung_RNA,
  features = c("SPI1", "TBX15", "E2F3"),
  group.by = "celltype",
  pt.size = 0
)

########################################################################################

library(data.table)
library(Seurat)
library(ggplot2)
library(patchwork)

out_dir <- "E:/Lung/EGA/PDF/gene_list/regulatory_validation"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

cell_order <- c(
  "AT1",
  "AT2",
  "Endothelial",
  "Fibroblast",
  "Lymphoid",
  "Myeloid",
  "Tumor"
)

cell_cols <- c(
  "AT1"         = "#F47A70",
  "AT2"         = "#D2A000",
  "Endothelial" = "#10B99A",
  "Fibroblast"  = "#10A8D0",
  "Lymphoid"    = "#9480E8",
  "Myeloid"     = "#F04BC8",
  "Tumor"       = "#57C000"
)



plot_tf_validation <- function(
    tf_name,
    subgrn,
    lung_RNA,
    celltype_col = "celltype"
) {
  
  #========================================
  # 1. subGRN data
  #========================================
  tf_sub <- subgrn[
    TF == tf_name &
      Celltype %in% cell_order
  ]
  
  if (nrow(tf_sub) == 0) {
    stop(paste0("subgrn 中没有 ", tf_name))
  }
  
  tf_sub[, Celltype := factor(
    Celltype,
    levels = cell_order
  )]
  
  #========================================
  # 2. TF expression violin
  #========================================
  p_expr <- VlnPlot(
    lung_RNA,
    features = tf_name,
    group.by = celltype_col,
    pt.size = 0
  ) +
    scale_fill_manual(values = cell_cols) +
    labs(
      title = paste0(tf_name, " expression"),
      x = NULL,
      y = "Expression level"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  #========================================
  # 3. Mean CellEdgeScore
  #========================================
  tf_profile <- tf_sub[, .(
    Mean_CellEdgeScore = mean(
      CellEdgeScore,
      na.rm = TRUE
    )
  ), by = Celltype]
  
  tf_profile[, Celltype := factor(
    Celltype,
    levels = cell_order
  )]
  
  p_edge <- ggplot(
    tf_profile,
    aes(
      x = Celltype,
      y = Mean_CellEdgeScore,
      fill = Celltype
    )
  ) +
    geom_col(
      width = 0.68,
      color = "black",
      linewidth = 0.4
    ) +
    scale_fill_manual(values = cell_cols) +
    labs(
      title = paste0(
        tf_name,
        " regulatory strength"
      ),
      x = NULL,
      y = "Mean CellEdgeScore"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  #========================================
  # 4. Peak CTS violin
  #========================================
  p_peak <- ggplot(
    tf_sub,
    aes(
      x = Celltype,
      y = Peak_CTS,
      fill = Celltype
    )
  ) +
    geom_violin(
      trim = FALSE,
      scale = "width"
    ) +
    geom_boxplot(
      width = 0.12,
      fill = "white",
      outlier.shape = NA,
      linewidth = 0.45
    ) +
    scale_fill_manual(values = cell_cols) +
    labs(
      title = paste0(
        tf_name,
        " peak specificity"
      ),
      x = NULL,
      y = "Peak CTS"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  #========================================
  # 5. Three-panel plot
  #========================================
  p_all <- p_expr | p_edge | p_peak
  
  return(
    list(
      expression = p_expr,
      edge = p_edge,
      peak = p_peak,
      combined = p_all,
      summary = tf_profile
    )
  )
}


library(data.table)
library(ggplot2)

#====================================================
# 1. Output directory
#====================================================

out_dir <- "E:/Lung/EGA/PDF/gene_list/regulatory_validation/significant"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


#====================================================
# 2. Filter all significant TF-by-cell-type pairs
#====================================================

sig_tf <- final_results[
  CTS_P < 0.05 &
    MGRS_P < 0.05 &
    MRS_P < 0.05
]

sig_tf <- unique(
  sig_tf[, .(
    TF,
    Celltype,
    CTS_P,
    MGRS_P,
    MRS_P,
    CTS_null_z,
    MGRS_null_z,
    MRS_null_z
  )]
)

sig_tf <- sig_tf[
  order(Celltype, -MRS_null_z)
]

# Inspect the number of pairs
nrow(sig_tf)

table(sig_tf$Celltype)

sig_tf



#====================================================
# 3. Batch-plot all significant TFs
#====================================================

for (i in seq_len(nrow(sig_tf))) {
  
  tf_i <- sig_tf$TF[i]
  ct_i <- sig_tf$Celltype[i]
  
  message(
    "[", i, "/", nrow(sig_tf), "] ",
    tf_i, " - ", ct_i
  )
  
  #--------------------------------------------------
  # Check whether the TF is present in the RNA data
  #--------------------------------------------------
  
  if (!tf_i %in% rownames(lung_RNA)) {
    
    message(
      "  跳过：",
      tf_i,
      " 不存在于 lung_RNA"
    )
    
    next
  }
  
  #--------------------------------------------------
  # Check subGRN data
  #--------------------------------------------------
  
  n_sub <- subgrn[
    TF == tf_i,
    .N
  ]
  
  if (n_sub == 0) {
    
    message(
      "  跳过：subgrn 中没有 ",
      tf_i
    )
    
    next
  }
  
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  
  res <- tryCatch(
    
    plot_tf_validation(
      tf_name = tf_i,
      subgrn = subgrn,
      lung_RNA = lung_RNA
    ),
    
    error = function(e) {
      
      message(
        "  绘图失败：",
        e$message
      )
      
      return(NULL)
    }
  )
  
  if (is.null(res)) {
    next
  }
  
  #--------------------------------------------------
  # File-system-safe name
  #--------------------------------------------------
  
  safe_tf <- gsub(
    "[\\\\/:*?\"<>| ]",
    "_",
    tf_i
  )
  
  safe_ct <- gsub(
    "[\\\\/:*?\"<>| ]",
    "_",
    ct_i
  )
  
  prefix <- paste0(
    sprintf("%03d", i),
    "_",
    safe_ct,
    "_",
    safe_tf
  )
  
  #--------------------------------------------------
  # Save the three-panel plot as a PDF
  #--------------------------------------------------
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        prefix,
        "_regulatory_validation.pdf"
      )
    ),
    plot = res$combined,
    width = 14,
    height = 4.8,
    units = "in",
    bg = "white"
  )
  
  #--------------------------------------------------
  # Optional: save regulatory strength separately
  #--------------------------------------------------
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        prefix,
        "_CellEdgeScore.pdf"
      )
    ),
    plot = res$edge,
    width = 5.5,
    height = 4.5,
    units = "in",
    bg = "white"
  )
  
  #--------------------------------------------------
  # Optional: save peak CTS separately
  #--------------------------------------------------
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        prefix,
        "_Peak_CTS.pdf"
      )
    ),
    plot = res$peak,
    width = 5.5,
    height = 4.5,
    units = "in",
    bg = "white"
  )
}

message(
  "全部完成：",
  out_dir
)

ggsave(
  file.path(
    out_dir,
    "SPI1_regulatory_validation.pdf"
  ),
  res_SPI1$combined,
  width = 14,
  height = 4.5,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "TBX15_regulatory_validation.pdf"
  ),
  res_TBX15$combined,
  width = 14,
  height = 4.5,
  bg = "white"
)

ggsave(
  file.path(
    out_dir,
    "E2F3_regulatory_validation.pdf"
  ),
  res_E2F3$combined,
  width = 14,
  height = 4.5,
  bg = "white"
)
































library(Seurat)
library(ggplot2)

p <- FeaturePlot(
  lung_RNA,
  features = "E2F3_Regulon1",
  reduction = "wnn.umap",
  pt.size = 0.85,
  order = TRUE,
  min.cutoff = "q02",
  max.cutoff = "q98"
) +
  scale_color_gradientn(
    colours = c(
      "#D0D0D0",
      "#F2D2C8",
      "#F5A08C",
      "#EF6A53",
      "#E53E2F",
      "#B51F16"
    )
  ) +
  labs(
    title = expression(italic(E2F3)*"-eRegulon")
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      color = "black"
    ),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 10,
      color = "black"
    )
  ) +
  guides(
    color = guide_colorbar(
      barwidth = 0.5,
      barheight = 5.5
    )
  )

p


ggsave(
  filename = "E:/Lung/E2F3_eRegulon_WNN.pdf",
  plot = p,
  width = 8,
  height = 6,
  units = "in",
  device = cairo_pdf
)
