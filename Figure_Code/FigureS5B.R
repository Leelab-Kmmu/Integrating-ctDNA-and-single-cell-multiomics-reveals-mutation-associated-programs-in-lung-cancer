library(data.table)
library(ggplot2)

# ============================================================
# Function: draw a bubble plot of top TFs for one cell type
# ============================================================

plot_ctmrs_bubble <- function(
    df,
    celltype_use,
    top_n = 7
) {
  
  # ----------------------------------------------------------
  # 1. Select top TFs with the highest MRS_null_z in this cell type
  # ----------------------------------------------------------
  
  tmp <- copy(df[
    Celltype == celltype_use &
      is.finite(MRS_null_z) &
      is.finite(CTS_null_z) &
      is.finite(MGRS_null_z)
  ])
  
  setorder(tmp, -MRS_null_z)
  
  tmp <- head(tmp, top_n)
  
  # ----------------------------------------------------------
  # 2. Significance symbols
  # ----------------------------------------------------------
  
  tmp[, Sign := fifelse(
    CTS_P < 0.05 &
      MGRS_P < 0.05 &
      MRS_P < 0.01,
    "**",
    
    fifelse(
      CTS_P < 0.05 &
        MGRS_P < 0.05 &
        MRS_P < 0.05,
      "*",
      
      fifelse(
        CTS_P < 0.05 &
          MGRS_P < 0.10 &
          MRS_P < 0.05,
        "#",
        ""
      )
    )
  )]
  
  # ----------------------------------------------------------
  # 3. Convert to long format
  # ----------------------------------------------------------
  
  plot_df <- melt(
    tmp,
    id.vars = c(
      "TF",
      "Celltype",
      "MRS_null_z",
      "MRS_P",
      "CTS_P",
      "MGRS_P",
      "Sign"
    ),
    measure.vars = c(
      "CTS_null_z",
      "MGRS_null_z"
    ),
    variable.name = "ScoreType",
    value.name = "Score"
  )
  
  plot_df[
    ScoreType == "CTS_null_z",
    ScoreType := "CTS"
  ]
  
  plot_df[
    ScoreType == "MGRS_null_z",
    ScoreType := "MGRS"
  ]
  
  # ----------------------------------------------------------
  # 4. Order TFs by MRS
  # ----------------------------------------------------------
  
  tf_order <- tmp$TF
  
  plot_df[, TF := factor(
    TF,
    levels = tf_order
  )]
  
  # Order of the upper and lower rows
  plot_df[, ScoreType := factor(
    ScoreType,
    levels = c(
      "MGRS",
      "CTS"
    )
  )]
  
  # Prevent negative values from causing size errors
  plot_df[, BubbleSize := pmax(Score, 0)]
  
  # ----------------------------------------------------------
  # 5. Plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    plot_df,
    aes(
      x = TF,
      y = ScoreType
    )
  ) +
    
    geom_point(
      aes(
        size = BubbleSize,
        fill = ScoreType
      ),
      shape = 21,
      color = "#EF8175",
      stroke = 0.8,
      alpha = 0.90
    ) +
    
    # Significance symbols
    geom_text(
      aes(label = Sign),
      color = "#3F3F3F",
      size = 5,
      vjust = 0.35
    ) +
    
    scale_fill_manual(
      values = c(
        "CTS" = "#65B5C5",
        "MGRS" = "#C6429A"
      ),
      name = NULL
    ) +
    
    scale_size_continuous(
      range = c(2.5, 9),
      name = "Null z-score"
    ) +
    
    labs(
      title = celltype_use,
      x = NULL,
      y = NULL
    ) +
    
    theme_minimal(base_size = 13) +
    
    theme(
      panel.grid.major.x = element_line(
        color = "#EFEFEF",
        linewidth = 1
      ),
      
      panel.grid.major.y = element_line(
        color = "#EFEFEF",
        linewidth = 1
      ),
      
      panel.grid.minor = element_blank(),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 11,
        color = "black"
      ),
      
      axis.text.y = element_text(
        size = 11,
        color = "black"
      ),
      
      axis.ticks = element_blank(),
      
      legend.position = "bottom",
      
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.margin = margin(
        15, 15, 15, 15
      )
    )
  
  return(p)
}
p_lymphoid <- plot_ctmrs_bubble(
  df = df_reg,
  celltype_use = "Lymphoid",
  top_n = 7
)

p_lymphoid


library(data.table)
library(ggplot2)
library(patchwork)

plot_ctmrs_bubble <- function(
    df,
    celltype_use,
    top_n = 7
) {
  
  tmp <- copy(
    df[
      Celltype == celltype_use &
        is.finite(MRS_null_z) &
        is.finite(CTS_null_z) &
        is.finite(MGRS_null_z)
    ]
  )
  
  setorder(tmp, -MRS_null_z)
  tmp <- head(tmp, top_n)
  
  # Significance symbols
  tmp[, Sign := fifelse(
    CTS_P < 0.05 &
      MGRS_P < 0.05 &
      MRS_P < 0.01,
    "**",
    
    fifelse(
      CTS_P < 0.05 &
        MGRS_P < 0.05 &
        MRS_P < 0.05,
      "*",
      
      fifelse(
        CTS_P < 0.05 &
          MGRS_P < 0.10 &
          MRS_P < 0.05,
        "#",
        ""
      )
    )
  )]
  
  # Convert to long format
  plot_df <- melt(
    tmp,
    id.vars = c(
      "TF",
      "Celltype",
      "MRS_null_z",
      "MRS_P",
      "CTS_P",
      "MGRS_P",
      "Sign"
    ),
    measure.vars = c(
      "CTS_null_z",
      "MGRS_null_z"
    ),
    variable.name = "ScoreType",
    value.name = "Score"
  )
  
  plot_df[
    ScoreType == "CTS_null_z",
    ScoreType := "CTS"
  ]
  
  plot_df[
    ScoreType == "MGRS_null_z",
    ScoreType := "MGRS"
  ]
  
  # TF order
  tf_order <- tmp$TF
  
  plot_df[, TF := factor(
    TF,
    levels = tf_order
  )]
  
  # Upper and lower rows
  plot_df[, ScoreType := factor(
    ScoreType,
    levels = c("MGRS", "CTS")
  )]
  
  # Point size: absolute z-score
  plot_df[, BubbleSize := abs(Score)]
  
  p <- ggplot(
    plot_df,
    aes(
      x = TF,
      y = ScoreType
    )
  ) +
    
    geom_point(
      aes(
        size = BubbleSize,
        fill = ScoreType
      ),
      shape = 21,
      color = "#EF8175",
      stroke = 0.8,
      alpha = 0.9
    ) +
    
    geom_text(
      aes(label = Sign),
      color = "#404040",
      size = 4.5,
      vjust = 0.35
    ) +
    
    scale_fill_manual(
      values = c(
        "CTS" = "#65B5C5",
        "MGRS" = "#C6429A"
      ),
      name = NULL
    ) +
    
    scale_size_continuous(
      range = c(2.5, 8),
      name = "|z-score|"
    ) +
    
    labs(
      title = celltype_use,
      x = NULL,
      y = NULL
    ) +
    
    theme_minimal(base_size = 12) +
    
    theme(
      panel.grid.major.x = element_line(
        color = "#EFEFEF",
        linewidth = 1
      ),
      
      panel.grid.major.y = element_line(
        color = "#EFEFEF",
        linewidth = 1
      ),
      
      panel.grid.minor = element_blank(),
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 10,
        color = "black"
      ),
      
      axis.text.y = element_text(
        size = 10,
        color = "black"
      ),
      
      axis.ticks = element_blank(),
      
      legend.position = "bottom",
      
      plot.title = element_text(
        size = 13,
        face = "bold",
        hjust = 0.5
      ),
      
      plot.margin = margin(
        10, 10, 10, 10
      )
    )
  
  return(p)
}


celltypes <- c(
  "AT1",
  "AT2",
  "Tumor",
  "Endothelial",
  "Fibroblast",
  "Lymphoid",
  "Myeloid"
)

plot_list <- lapply(
  celltypes,
  function(ct) {
    plot_ctmrs_bubble(
      df = df_reg,
      celltype_use = ct,
      top_n = 7
    )
  }
)

names(plot_list) <- celltypes


outdir <- "E:/Lung/EGA/PDF"

for (ct in celltypes) {
  
  ggsave(
    filename = file.path(
      outdir,
      paste0(ct, "_ctMRS_bubble.svg")
    ),
    plot = plot_list[[ct]],
    width = 7.5,
    height = 4.5,
    units = "in"
  )
}
