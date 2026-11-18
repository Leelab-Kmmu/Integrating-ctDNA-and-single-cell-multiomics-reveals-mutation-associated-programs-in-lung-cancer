library(maftools)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(tidyr)
library(ggplot2)


# ==========================================
# 1. Set the output directory
# ==========================================

outdir <- "E:/Lung/Reactome_analysis"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==========================================
# 2. MAF and clinical information
# ==========================================

maf_data <- as.data.frame(maf_lung@data)

clinical <- unique(
  as.data.frame(maf_lung@clinical.data)[
    , c("Tumor_Sample_Barcode", "phenotype_group")
  ]
)

group_order <- c(
  "Lung cancer",
  "Metastasis to lung",
  "Benign tumor"
)

group_n <- clinical %>%
  dplyr::count(
    phenotype_group,
    name = "Total_samples"
  )

# ==========================================
# 3. Number of patients with gene mutations in each group
# ==========================================

gene_freq <- maf_data %>%
  dplyr::distinct(
    Tumor_Sample_Barcode,
    Hugo_Symbol
  ) %>%
  dplyr::inner_join(
    clinical,
    by = "Tumor_Sample_Barcode"
  ) %>%
  dplyr::count(
    phenotype_group,
    Hugo_Symbol,
    name = "Mutated_samples"
  ) %>%
  dplyr::left_join(
    group_n,
    by = "phenotype_group"
  ) %>%
  dplyr::mutate(
    Frequency = Mutated_samples / Total_samples
  )

# ==========================================
# 4. Filter candidate genes
# ==========================================

candidate_genes <- gene_freq %>%
  dplyr::filter(
    Frequency >= 0.10,
    Mutated_samples >= 3
  )

# Inspect the number of candidate genes in each group
candidate_genes %>%
  dplyr::count(phenotype_group)

write.csv(
  candidate_genes,
  file.path(outdir, "Recurrent_mutated_genes.csv"),
  row.names = FALSE
)


# All candidate genes
all_candidate_symbols <- unique(
  candidate_genes$Hugo_Symbol
)

# Use all genes in the current MAF as the initial background
background_symbols <- unique(
  maf_data$Hugo_Symbol
)

# Convert all genes in one step
gene_mapping <- clusterProfiler::bitr(
  background_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = TRUE
)

gene_mapping <- unique(
  gene_mapping[, c("SYMBOL", "ENTREZID")]
)

background_ids <- unique(
  gene_mapping$ENTREZID
)

# Convert the three groups to Entrez IDs separately
candidate_ids <- candidate_genes %>%
  dplyr::inner_join(
    gene_mapping,
    by = c("Hugo_Symbol" = "SYMBOL")
  ) %>%
  dplyr::distinct(
    phenotype_group,
    ENTREZID
  )

gene_list <- split(
  candidate_ids$ENTREZID,
  candidate_ids$phenotype_group
)

gene_list <- lapply(gene_list, unique)

print(lengths(gene_list))



# ==========================================
# 5. Perform enrichment analysis for each group separately
# ==========================================

reactome_results <- list()

for (grp in group_order) {
  
  genes <- gene_list[[grp]]
  
  if (is.null(genes) || length(genes) < 10) {
    message("Skipping ", grp, ": insufficient genes")
    next
  }
  
  result <- ReactomePA::enrichPathway(
    gene = genes,
    organism = "human",
    universe = background_ids,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
  
  result_df <- as.data.frame(result)
  
  if (nrow(result_df) == 0) next
  
  result_df$phenotype_group <- grp
  
  reactome_results[[grp]] <- result_df
}

# Merge
reactome_all <- dplyr::bind_rows(
  reactome_results
)

if (nrow(reactome_all) == 0) {
  stop("No Reactome pathways returned.")
}

# Calculate -log10(FDR)
reactome_all$negLog10FDR <- -log10(
  pmax(reactome_all$p.adjust, .Machine$double.xmin)
)

write.csv(
  reactome_all,
  file.path(outdir, "Reactome_enrichment_all.csv"),
  row.names = FALSE
)

# Inspect the number of significant pathways in each group
reactome_all %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::count(phenotype_group)




# All candidate genes
all_candidate_symbols <- unique(
  candidate_genes$Hugo_Symbol
)

# Use all genes in the current MAF as the initial background
background_symbols <- unique(
  maf_data$Hugo_Symbol
)

# Convert all genes in one step
gene_mapping <- clusterProfiler::bitr(
  background_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db,
  drop = TRUE
)

gene_mapping <- unique(
  gene_mapping[, c("SYMBOL", "ENTREZID")]
)

background_ids <- unique(
  gene_mapping$ENTREZID
)

# Convert the three groups to Entrez IDs separately
candidate_ids <- candidate_genes %>%
  dplyr::inner_join(
    gene_mapping,
    by = c("Hugo_Symbol" = "SYMBOL")
  ) %>%
  dplyr::distinct(
    phenotype_group,
    ENTREZID
  )

gene_list <- split(
  candidate_ids$ENTREZID,
  candidate_ids$phenotype_group
)

gene_list <- lapply(gene_list, unique)

print(lengths(gene_list))


# ==========================================
# 5. Perform enrichment analysis for each group separately
# ==========================================

reactome_results <- list()

for (grp in group_order) {
  
  genes <- gene_list[[grp]]
  
  if (is.null(genes) || length(genes) < 10) {
    message("Skipping ", grp, ": insufficient genes")
    next
  }
  
  result <- ReactomePA::enrichPathway(
    gene = genes,
    organism = "human",
    universe = background_ids,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
  )
  
  result_df <- as.data.frame(result)
  
  if (nrow(result_df) == 0) next
  
  result_df$phenotype_group <- grp
  
  reactome_results[[grp]] <- result_df
}

# Merge
reactome_all <- dplyr::bind_rows(
  reactome_results
)

if (nrow(reactome_all) == 0) {
  stop("No Reactome pathways returned.")
}

# Calculate -log10(FDR)
reactome_all$negLog10FDR <- -log10(
  pmax(reactome_all$p.adjust, .Machine$double.xmin)
)

write.csv(
  reactome_all,
  file.path(outdir, "Reactome_enrichment_all.csv"),
  row.names = FALSE
)

# Inspect the number of significant pathways in each group
reactome_all %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::count(phenotype_group)


library(ggplot2)
library(dplyr)

# ==========================================
# 6. Filter significant pathways
# ==========================================

plot_data <- reactome_all %>%
  dplyr::filter(
    !is.na(p.adjust),
    p.adjust < 0.05
  )

if (nrow(plot_data) == 0) {
  stop(
    "No pathways meet FDR < 0.05. ",
    "Inspect the results rather than relaxing thresholds solely for plotting."
  )
}

# Select the 10 pathways with the lowest FDR in each group
plot_data <- plot_data %>%
  dplyr::group_by(phenotype_group) %>%
  dplyr::slice_min(
    order_by = p.adjust,
    n = 10,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

# Remove duplicate pathway-group records
plot_data <- plot_data %>%
  dplyr::distinct(
    phenotype_group,
    ID,
    .keep_all = TRUE
  )

# Sort by significance
pathway_order <- plot_data %>%
  dplyr::group_by(Description) %>%
  dplyr::summarise(
    Max_score = max(negLog10FDR),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Max_score) %>%
  dplyr::pull(Description)

plot_data$Description <- factor(
  plot_data$Description,
  levels = pathway_order
)

plot_data$phenotype_group <- factor(
  plot_data$phenotype_group,
  levels = group_order
)

# ==========================================
# 7. Plot
# ==========================================

p_reactome <- ggplot(
  plot_data,
  aes(
    x = phenotype_group,
    y = Description
  )
) +
  geom_point(
    aes(
      size = Count,
      color = negLog10FDR
    ),
    alpha = 0.9
  ) +
  scale_color_gradient(
    low = "#85B5C3",
    high = "#B94F4F"
  ) +
  scale_size_continuous(
    range = c(3, 10)
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "Gene count",
    color = expression(-log[10](FDR))
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      size = 11,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    legend.position = "right",
    panel.grid.major.y = element_line(
      color = "grey92",
      linewidth = 0.3
    )
  )

print(p_reactome)

ggsave(
  file.path(outdir, "Reactome_enrichment_dotplot.pdf"),
  plot = p_reactome,
  width = 9,
  height = max(6, 0.28 * length(pathway_order) + 2),
  device = cairo_pdf
)
