library(maftools)

# 1. Extract all target genes for E2F3-Tumor
E2F3_targets <- unique(
  subgrn[
    TF == "E2F3" &
      Celltype == "Tumor",
    gene
  ]
)

E2F3_targets <- E2F3_targets[
  !is.na(E2F3_targets) &
    E2F3_targets != ""
]

length(E2F3_targets)
E2F3_targets


oncoplot(
  maf = maf_lung,
  genes = E2F3_targets,
  removeNonMutated = FALSE,
  showTumorSampleBarcodes = FALSE
)



genes_E2F3_oncoplot <- unique(
  c(
    "E2F3",
    E2F3_targets
  )
)

genes_E2F3_oncoplot <- c(
  "E2F3",
  "AKR1C2",
  "AKR1C3",
  "CDKN2A",
  "COL4A4",
  "CPS1",
  "DIP2C",
  "GSTCD",
  "IGF2BP3",
  "IGSF1",
  "KIAA0319",
  "NPAS2",
  "NPNT",
  "PCSK6",
  "PPARGC1A",
  "PTPRH",
  "RGS17",
  "VCAN",
  "WEE1"
)

length(genes_E2F3_oncoplot)
pdf(
  file.path(
    outdir,
    "E2F3_Tumor_target_oncoplot.pdf"
  ),
  width = 17,
  height = 8,
  useDingbats = FALSE
)

oncoplot(
  maf = maf_lung,
  genes = genes_E2F3_oncoplot,
  clinicalFeatures = clinicalFeatures,
  colors = variant_colors,
  annotationColor = ann_colors,
  sortByAnnotation = TRUE,
  showTumorSampleBarcodes = FALSE,
  drawRowBar = TRUE,
  drawColBar = TRUE,
  removeNonMutated = FALSE,
  
  gene_mar = 12,
  fontSize = 0.9,
  
  titleText = "Mutation landscape of E2F3 tumor-specific regulon targets"
)

dev.off()
