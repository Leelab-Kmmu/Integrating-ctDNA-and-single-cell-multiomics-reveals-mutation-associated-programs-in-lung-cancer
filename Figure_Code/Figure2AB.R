library(data.table)
library(maftools)

# ============================================================
# 1. Paths
# ============================================================

variant_file <- "E:/Lung/EGA/meatdata/EGA_lung_variants.txt"
clinical_file <- "E:/Lung/EGA/meatdata/mutation_sample_clinical_phenotype.csv"

outdir <- "E:/Lung/EGA/meatdata/gene_list"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Read mutation data
# ============================================================
library(data.table)
library(maftools)
mut <- fread(variant_file)

head(mut)
names(mut)
dim(mut)
# ============================================================
# 1. Prepare the MAF data directly
# ============================================================

maf_dt <- mut[
  !is.na(Tumor_Sample_Barcode) &
    Tumor_Sample_Barcode != "" &
    !is.na(Hugo_Symbol) &
    Hugo_Symbol != ""
]

# Retain fields required or commonly used by maftools
maf_dt <- maf_dt[
  ,
  .(
    Hugo_Symbol,
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele2,
    Tumor_Sample_Barcode,
    Variant_Classification,
    Variant_Type
  )
]

dim(maf_dt)
head(maf_dt)


table(
  maf_dt$Variant_Classification,
  useNA = "ifany"
)

table(
  maf_dt$Variant_Type,
  useNA = "ifany"
)



#clinical_file <- "/400T/ckn/Lung_Model/Final_GRN/EGA/mutation_sample_clinical_phenotype.csv"

clinical <- fread(clinical_file)

head(clinical)
names(clinical)
dim(clinical)




library(data.table)

clinical[, phenotype_group := fifelse(
  phenotype %in% c(
    "lung cancer",
    "lung cancer, later lung cancer",
    "lung cancer, prior lung cancer",
    "no lung cancer, later lung cancer",
    "no lung cancer, prior lung cancer",
    "cancer",
    "no lung cancer, prior lung cancer, later lung cancer"
  ),
  "Lung cancer",
  
  fifelse(
    phenotype %in% c(
      "lung cancer, other cancer",
      "no lung cancer, other cancer, later lung cancer"
    ),
    "Lung cancer + other cancer",
    
    fifelse(
      phenotype %in% c(
        "met to lung",
        "met to lung, other cancer"
      ),
      "Metastasis to lung",
      
      fifelse(
        phenotype == "No baseline cancer",
        "No baseline cancer",
        
        fifelse(
          phenotype %in% c(
            "no lung cancer, other cancer"
            
          ),
          "Other cancer",
          
          fifelse(
            phenotype == "良性肿瘤",
            "Benign tumor",
            
            fifelse(
              phenotype == "Human peripheral blood monocytes",
              "PBMC",
              "Other"
            )
          )
        )
      )
    )
  )
)]


table(
  clinical$phenotype_group,
  useNA = "ifany"
)


clinicalFeatures <- c(
  "phenotype_group",
  "biological_sex"
)


ann_colors <- list(
  phenotype_group = c(
    "Lung cancer" = "#C95B5B",
    "Lung cancer + other cancer" = "#A63603",
    "Metastasis to lung" = "#7A0177",
    "Other cancer" = "#807DBA",
    "No baseline lung cancer" = "#BDBDBD",
    "Benign tumor" = "#5AAE61",
    "PBMC" = "#4DBBD5",
    "Other" = "#D9D9D9"
  ),
  
  biological_sex = c(
    "male" = "#4DBBD5",
    "female" = "#E78AC3"
  )
)




library(data.table)
library(maftools)

# ============================================================
# 1. Standardize biological_sex
# ============================================================

clinical[
  biological_sex == "male",
  biological_sex := "Male"
]

clinical[
  biological_sex == "female",
  biological_sex := "Female"
]

# Check
table(clinical$biological_sex, useNA = "ifany")
table(clinical$phenotype_group, useNA = "ifany")


# ============================================================
# 2. Create the sample ID column required by maftools
# ============================================================

clinical[
  ,
  Tumor_Sample_Barcode := sample_id
]


# ============================================================
# 3. Check sample matching between mutation and clinical data
# ============================================================

cat(
  "Mutation samples:",
  uniqueN(mut$Tumor_Sample_Barcode),
  "\n"
)

cat(
  "Clinical samples:",
  uniqueN(clinical$Tumor_Sample_Barcode),
  "\n"
)

matched_samples <- intersect(
  unique(mut$Tumor_Sample_Barcode),
  unique(clinical$Tumor_Sample_Barcode)
)

cat(
  "Matched samples:",
  length(matched_samples),
  "\n"
)

# Inspect the phenotype distribution of samples included in the analysis
table(
  clinical[
    Tumor_Sample_Barcode %in% matched_samples,
    phenotype_group
  ],
  useNA = "ifany"
)



# ============================================================
# 4. Prepare the mutation table
# ============================================================

maf_dt <- mut[
  !is.na(Tumor_Sample_Barcode) &
    Tumor_Sample_Barcode != "" &
    !is.na(Hugo_Symbol) &
    Hugo_Symbol != ""
]

maf_dt <- maf_dt[
  ,
  .(
    Hugo_Symbol,
    Chromosome,
    Start_Position,
    End_Position,
    Reference_Allele,
    Tumor_Seq_Allele2,
    Tumor_Sample_Barcode,
    Variant_Classification,
    Variant_Type
  )
]

# Check
table(
  maf_dt$Variant_Classification,
  useNA = "ifany"
)

table(
  maf_dt$Variant_Type,
  useNA = "ifany"
)



# ============================================================
# 5. Create the MAF object
# ============================================================
library(maftools)
maf_lung <- read.maf(
  maf = maf_dt,
  clinicalData = clinical,
  verbose = TRUE
)

maf_lung

head(maf_lung@clinical.data)



# ============================================================
# 6. Clinical features
# ============================================================

clinicalFeatures <- c(
  "phenotype_group",
  "biological_sex"
)

clinicalFeatures_use <- clinicalFeatures[
  clinicalFeatures %in%
    names(maf_lung@clinical.data)
]

clinicalFeatures_use


# ============================================================
# 7. Mutation colors
# ============================================================

variant_colors <- c(
  Missense_Mutation = "#377EB8",
  Nonsense_Mutation = "#E64B35",
  Frame_Shift_Del = "#4DBBD5",
  Frame_Shift_Ins = "#00A087",
  In_Frame_Del = "#F39B7F",
  In_Frame_Ins = "#8491B4",
  Splice_Site = "#91D1C2",
  Translation_Start_Site = "#DC0000",
  Nonstop_Mutation = "#7E6148",
  Multi_Hit = "#8DD3C7",
  Silent = "#BDBDBD"
)


# ============================================================
# 8. Annotation colors
# ============================================================

ann_colors <- list(
  
  phenotype_group = c(
    "Lung cancer" = "#C95B5B",
    "Lung cancer + other cancer" = "#A63603",
    "Metastasis to lung" = "#7A0177",
    "No baseline cancer" = "#BDBDBD",
    "Other cancer" = "#807DBA",
    "Benign tumor" = "#5AAE61",
    "PBMC" = "#4DBBD5"
  ),
  
  biological_sex = c(
    "Male" = "#4DBBD5",
    "Female" = "#E78AC3"
  )
)
ann_colors_use <- ann_colors[
  clinicalFeatures_use
]



# ============================================================
# 9. Output
# ============================================================

outdir <- "E:/Lung/EGA/PDF/gene_list"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 10. Oncoplot
# ============================================================

pdf(
  file.path(
    outdir,
    "Figure2A_lung_oncoplot.pdf"
  ),
  width = 16,
  height = 9,
  useDingbats = FALSE
)

oncoplot(
  maf = maf_lung,
  
  top = 30,
  
  clinicalFeatures = clinicalFeatures_use,
  
  colors = variant_colors,
  
  annotationColor = ann_colors_use,
  
  sortByAnnotation = TRUE,
  
  showTumorSampleBarcodes = FALSE,
  
  drawRowBar = TRUE,
  
  drawColBar = TRUE,
  
  removeNonMutated = TRUE,
  gene_mar = 12,   # Important: increase the left margin for gene names
  fontSize = 0.9,
  
  titleText = "Clinical annotated ctDNA mutation landscape"
)

dev.off()






library(data.table)
library(maftools)

gene_summary <- getGeneSummary(maf_lung)

top35_genes <- gene_summary[
  !grepl("^LOC", Hugo_Symbol)
][
  order(-MutatedSamples),
  Hugo_Symbol
][1:35]

top35_genes
lung_driver_genes <- c(
  "TP53",
  "EGFR",
  "KRAS",
  "STK11",
  "KEAP1",
  "NF1",
  "RB1",
  "BRAF",
  "PIK3CA",
  "MET",
  "ERBB2",
  "ALK",
  "RET",
  "NRAS",
  "SMARCA4"
)

lung_driver_present <- lung_driver_genes[
  lung_driver_genes %in% gene_summary$Hugo_Symbol
]

genes_oncoplot <- unique(
  c(
    top35_genes,
    lung_driver_present
  )
)

genes_oncoplot
length(genes_oncoplot)



pdf(
  file.path(
    outdir,
    "Figure2A_Top35_plus_lung_driver_oncoplot.pdf"
  ),
  width = 20,
  height = 11,
  useDingbats = FALSE
)

oncoplot(
  maf = maf_lung,
  genes = genes_oncoplot,
  clinicalFeatures = clinicalFeatures,
  colors = variant_colors,
  annotationColor = ann_colors,
  sortByAnnotation = TRUE,
  showTumorSampleBarcodes = FALSE,
  drawRowBar = TRUE,
  drawColBar = TRUE,
  removeNonMutated = TRUE,
  
  gene_mar = 12,   # Important: increase the left margin for gene names
  fontSize = 0.9,
  
  titleText = "Clinical annotated ctDNA mutation landscape"
)

dev.off()























oncoplot(
  maf = maf_lung,
  
  top = 30,
  
  #clinicalFeatures = clinicalFeatures_use,
  
  colors = variant_colors,
  
  annotationColor = ann_colors_use,
  
  #sortByAnnotation = TRUE,
  
  showTumorSampleBarcodes = FALSE,
  
  drawRowBar = TRUE,
  
  drawColBar = TRUE,
  
  removeNonMutated = TRUE,
  gene_mar = 12,   # Important: increase the left margin for gene names
  fontSize = 0.9,
  
  titleText = "Clinical annotated ctDNA mutation landscape"
)
