# ============================================================================
# Final GRN complete pipeline
# RNA + ATAC + JASPAR2020 + somatic mutation -> CTS + MGRS + MRS
#
# Run in a fresh R session. This script does not use Signac::LinkPeaks.
# All intermediate and final results are written to rebuild_results/ so the
# existing process files in Final_GRN are not overwritten.
# ============================================================================

options(stringsAsFactors = FALSE)

# ----------------------------------------------------------------------------
# Step 0. Parameters
# ----------------------------------------------------------------------------
rna_rds <- "E:/Lung/GSE304587/lung_RNA.rds"
atac_rds <- "E:/Lung/lung_ATAC_with_JASPAR_motifs.rds"
gene_annotation_rds <- "E:/Lung/GRN/gene_annotation_grch38.rds"
mutation_file <- "E:/Lung/Lung_all_sample_multianno.txt"

project_dir <- "E:/Lung/GRN/Final_GRN"
output_dir <- file.path(project_dir, "rebuild_results")

rna_assay <- "RNA"
peak_assay <- "ATAC"
celltype_column <- "celltype"

variable_genes <- 5000L
peak_gene_distance <- 500000L
correlation_chunk_size <- 2000L
min_peak_gene_cor <- 0.05
max_peak_gene_fdr <- 0.05
min_targets_per_tf <- 5L
cell_edge_quantile <- 0.50

mutation_peak_buffer <- 0L
min_recurrent_samples <- 1L
theta <- 0.5
top_n <- 5L
perm_n <- 1000L       # Use 10000L for the final publication run.
seed <- 1234L
control_pool_n <- 50L
min_control_pool <- 20L
nominal_p <- 0.05
fdr_cutoff <- 0.05

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# Step 1. Packages and input checks
# ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

required_files <- c(
  rna_rds,
  atac_rds,
  gene_annotation_rds,
  mutation_file
)

if (!all(file.exists(required_files))) {
  stop(
    "Missing input files: ",
    paste(required_files[!file.exists(required_files)], collapse = ", ")
  )
}

message("Step 1 complete: inputs and packages checked.")

# ----------------------------------------------------------------------------
# Step 2. Load RNA/ATAC, retain paired cells and preserve ChromatinAssay
# ----------------------------------------------------------------------------
lung_RNA <- readRDS(rna_rds)
lung_ATAC <- readRDS(atac_rds)
gene_annotation <- readRDS(gene_annotation_rds)

if (!rna_assay %in% SeuratObject::Assays(lung_RNA)) {
  stop("RNA assay not found: ", rna_assay)
}
if (!peak_assay %in% SeuratObject::Assays(lung_ATAC)) {
  stop("ATAC assay not found: ", peak_assay)
}
if (!inherits(lung_ATAC[[peak_assay]], "ChromatinAssay")) {
  stop("ATAC assay is not a Signac ChromatinAssay.")
}
if (!celltype_column %in% colnames(lung_RNA[[]])) {
  stop("RNA metadata lacks cell-type column: ", celltype_column)
}

common_cells <- intersect(colnames(lung_RNA), colnames(lung_ATAC))
if (!length(common_cells)) stop("RNA and ATAC have no shared cell barcodes.")

lung_RNA <- subset(lung_RNA, cells = common_cells)
lung_ATAC <- subset(lung_ATAC, cells = colnames(lung_RNA))
if (!identical(colnames(lung_RNA), colnames(lung_ATAC))) {
  lung_ATAC <- lung_ATAC[, colnames(lung_RNA)]
}
stopifnot(identical(colnames(lung_RNA), colnames(lung_ATAC)))

if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  lung_RNA <- SeuratObject::JoinLayers(lung_RNA, assay = rna_assay)
}

single_cell <- lung_ATAC
single_cell[[rna_assay]] <- lung_RNA[[rna_assay]]
single_cell[[celltype_column]] <- lung_RNA[[celltype_column, drop = TRUE]]

celltype_vec <- as.character(single_cell[[celltype_column, drop = TRUE]])
if (anyNA(celltype_vec) || any(!nzchar(celltype_vec))) {
  stop("Cell-type annotations contain NA or empty values.")
}

rna_data <- SeuratObject::LayerData(
  single_cell,
  assay = rna_assay,
  layer = "data"
)
atac_counts <- SeuratObject::LayerData(
  single_cell,
  assay = peak_assay,
  layer = "counts"
)

if (!identical(colnames(rna_data), colnames(atac_counts))) {
  atac_counts <- atac_counts[, colnames(rna_data), drop = FALSE]
}
stopifnot(identical(colnames(rna_data), colnames(atac_counts)))

peak_ids <- rownames(atac_counts)
if (is.null(peak_ids) || anyDuplicated(peak_ids)) {
  stop("ATAC peak IDs are missing or duplicated.")
}

message(
  "Step 2 complete: paired cells = ", ncol(rna_data),
  "; genes = ", nrow(rna_data),
  "; peaks = ", nrow(atac_counts)
)

# ----------------------------------------------------------------------------
# Step 3. Calculate gene cell-type specificity (gene CTS)
# ----------------------------------------------------------------------------
celltype_groups <- factor(celltype_vec)
celltype_membership <- Matrix::sparseMatrix(
  i = seq_along(celltype_groups),
  j = as.integer(celltype_groups),
  x = 1,
  dims = c(length(celltype_groups), nlevels(celltype_groups)),
  dimnames = list(colnames(rna_data), levels(celltype_groups))
)

average_expression <- rna_data %*% celltype_membership
average_expression <- sweep(
  as.matrix(average_expression),
  2L,
  as.numeric(table(celltype_groups)),
  "/"
)

gene_cts_denominator <- sqrt(rowSums(average_expression^2))
gene_cts_denominator[
  !is.finite(gene_cts_denominator) | gene_cts_denominator == 0
] <- 1
gene_cts_matrix <- average_expression / gene_cts_denominator

target_scores <- data.table::data.table(
  gene = rep(rownames(gene_cts_matrix), times = ncol(gene_cts_matrix)),
  Celltype = rep(colnames(gene_cts_matrix), each = nrow(gene_cts_matrix)),
  Gene_CTS = as.vector(gene_cts_matrix)
)

data.table::fwrite(
  target_scores,
  file.path(output_dir, "00_gene_celltype_CTS.tsv"),
  sep = "\t"
)

message("Step 3 complete: gene CTS calculated for ", nlevels(celltype_groups), " cell types.")

# ----------------------------------------------------------------------------
# Step 4. Extract JASPAR motif matches: TF -> peak
# ----------------------------------------------------------------------------
motif_obj <- Signac::Motifs(single_cell[[peak_assay]])
if (is.null(motif_obj)) stop("ATAC assay has no stored Motif object.")

motif_matrix <- motif_obj@data
if (!identical(rownames(motif_matrix), peak_ids)) {
  stop("Motif matrix rows are not identical to ATAC peak rows.")
}

motif_names <- motif_obj@motif.names
if (!length(motif_names) || is.null(names(motif_names))) {
  stop("Motif object lacks motif-to-TF names.")
}

motif_map_parts <- vector("list", length(motif_names))
names(motif_map_parts) <- names(motif_names)
for (motif_id_now in names(motif_names)) {
  motif_map_parts[[motif_id_now]] <- data.table::data.table(
    motif_id = motif_id_now,
    TF_raw = as.character(motif_names[[motif_id_now]])
  )
}

motif_tf_map <- data.table::rbindlist(motif_map_parts, use.names = TRUE)
motif_tf_map <- motif_tf_map[, .(
  TF = unlist(strsplit(TF_raw, "::", fixed = TRUE))
), by = motif_id]
motif_tf_map[, TF := trimws(sub("\\(.*$", "", TF))]
motif_tf_map <- unique(motif_tf_map[nzchar(TF), .(motif_id, TF)])

rna_gene_lookup <- setNames(rownames(rna_data), toupper(rownames(rna_data)))
motif_tf_map[, TF_RNA := unname(rna_gene_lookup[toupper(TF)])]
motif_tf_map <- unique(
  motif_tf_map[!is.na(TF_RNA), .(motif_id, TF = TF_RNA)]
)

motif_sparse_summary <- Matrix::summary(motif_matrix)
tf_peak <- data.table::data.table(
  motif_id = colnames(motif_matrix)[motif_sparse_summary$j],
  peak_id = rownames(motif_matrix)[motif_sparse_summary$i]
)
tf_peak <- merge(
  unique(tf_peak),
  motif_tf_map,
  by = "motif_id",
  allow.cartesian = TRUE
)
tf_peak <- unique(tf_peak[, .(motif_id, peak_id, TF)])

data.table::fwrite(
  tf_peak,
  file.path(output_dir, "01_TF_peak.tsv.gz"),
  sep = "\t"
)

message(
  "Step 4 complete: TF-peak motif rows = ", nrow(tf_peak),
  "; TFs = ", uniqueN(tf_peak$TF)
)

# ----------------------------------------------------------------------------
# Step 5. Calculate full-genome peak-gene links without LinkPeaks
# ----------------------------------------------------------------------------
SeuratObject::DefaultAssay(single_cell) <- rna_assay
single_cell <- Seurat::FindVariableFeatures(
  single_cell,
  assay = rna_assay,
  selection.method = "vst",
  nfeatures = variable_genes,
  verbose = FALSE
)
genes_use <- intersect(
  SeuratObject::VariableFeatures(single_cell[[rna_assay]]),
  rownames(rna_data)
)
if (!length(genes_use)) stop("No variable gene is present in the RNA data matrix.")

atac_depth <- Matrix::colSums(atac_counts)
if (any(atac_depth <= 0)) stop("ATAC matrix contains cells with zero depth.")

atac_norm <- t(t(atac_counts) / atac_depth)
atac_norm <- atac_norm * 1e4
atac_norm@x <- log1p(atac_norm@x)

peak_gr <- Signac::granges(single_cell[[peak_assay]])
if (length(peak_gr) != nrow(atac_norm)) {
  stop("Peak ranges and ATAC matrix have different lengths.")
}
S4Vectors::mcols(peak_gr)$peak <- peak_ids

gene_anno_dt <- data.table::data.table(
  chr = as.character(GenomicRanges::seqnames(gene_annotation)),
  start = GenomicRanges::start(gene_annotation),
  end = GenomicRanges::end(gene_annotation),
  strand = as.character(GenomicRanges::strand(gene_annotation)),
  gene = as.character(S4Vectors::mcols(gene_annotation)$gene_name),
  type = as.character(S4Vectors::mcols(gene_annotation)$type)
)

gene_anno_dt <- gene_anno_dt[
  type == "exon" & !is.na(gene) & nzchar(gene) &
    gene %in% genes_use & strand %in% c("+", "-")
]
if (!nrow(gene_anno_dt)) stop("No variable gene has a valid exon annotation.")

gene_dt <- gene_anno_dt[, .(
  gene_start = min(start),
  gene_end = max(end)
), by = .(gene, chr, strand)]
gene_dt[, gene_width := gene_end - gene_start + 1L]
data.table::setorder(gene_dt, gene, -gene_width)
gene_dt <- gene_dt[!duplicated(gene)]
gene_dt[, gene_width := NULL]
gene_dt[, tss := ifelse(strand == "+", gene_start, gene_end)]

gene_window_dt <- gene_dt[, .(
  chr,
  start = pmax(1L, tss - peak_gene_distance),
  end = tss + peak_gene_distance,
  gene,
  tss
)]

peak_dt <- data.table::data.table(
  peak_idx = seq_along(peak_gr),
  peak = peak_ids,
  chr = as.character(GenomicRanges::seqnames(peak_gr)),
  start = GenomicRanges::start(peak_gr),
  end = GenomicRanges::end(peak_gr)
)
peak_dt[, peak_center := as.integer(round((start + end) / 2))]

data.table::setkey(gene_window_dt, chr, start, end)
peak_gene_pairs <- data.table::foverlaps(
  x = peak_dt,
  y = gene_window_dt,
  by.x = c("chr", "start", "end"),
  by.y = c("chr", "start", "end"),
  type = "any",
  nomatch = 0L
)
peak_gene_pairs[, distance := peak_center - tss]
peak_gene_pairs <- unique(peak_gene_pairs[, .(
  peak_idx,
  peak,
  gene,
  distance
)])
peak_gene_pairs <- peak_gene_pairs[
  peak %in% rownames(atac_norm) & gene %in% rownames(rna_data)
]

message(
  "Candidate peak-gene pairs = ", nrow(peak_gene_pairs),
  "; genes = ", uniqueN(peak_gene_pairs$gene)
)

pair_chunk <- ceiling(seq_len(nrow(peak_gene_pairs)) / correlation_chunk_size)
pair_chunks <- split(seq_len(nrow(peak_gene_pairs)), pair_chunk)
peak_gene_result <- vector("list", length(pair_chunks))
n_cells <- ncol(rna_data)

for (chunk_i in seq_along(pair_chunks)) {
  idx <- pair_chunks[[chunk_i]]
  pairs_now <- peak_gene_pairs[idx]

  peak_row_idx <- match(pairs_now$peak, rownames(atac_norm))
  gene_row_idx <- match(pairs_now$gene, rownames(rna_data))

  peak_block <- atac_norm[peak_row_idx, , drop = FALSE]
  gene_block <- rna_data[gene_row_idx, , drop = FALSE]

  peak_block_sq <- peak_block
  peak_block_sq@x <- peak_block_sq@x^2
  gene_block_sq <- gene_block
  gene_block_sq@x <- gene_block_sq@x^2

  sum_peak <- Matrix::rowSums(peak_block)
  sum_gene <- Matrix::rowSums(gene_block)
  sum_peak_sq <- Matrix::rowSums(peak_block_sq)
  sum_gene_sq <- Matrix::rowSums(gene_block_sq)
  sum_cross <- Matrix::rowSums(peak_block * gene_block)

  covariance_numerator <- sum_cross - sum_peak * sum_gene / n_cells
  peak_ss <- pmax(sum_peak_sq - sum_peak^2 / n_cells, 0)
  gene_ss <- pmax(sum_gene_sq - sum_gene^2 / n_cells, 0)
  correlation <- covariance_numerator / sqrt(peak_ss * gene_ss)
  correlation[!is.finite(correlation)] <- NA_real_
  correlation <- pmax(pmin(correlation, 0.999999), -0.999999)

  t_stat <- correlation * sqrt(
    (n_cells - 2) / pmax(1 - correlation^2, .Machine$double.eps)
  )
  p_value <- 2 * stats::pt(-abs(t_stat), df = n_cells - 2)

  peak_gene_result[[chunk_i]] <- data.table::data.table(
    peak = pairs_now$peak,
    gene = pairs_now$gene,
    distance = pairs_now$distance,
    cor = correlation,
    pvalue = p_value
  )[is.finite(cor) & is.finite(pvalue)]

  if (chunk_i %% 10L == 0L || chunk_i == length(pair_chunks)) {
    message("Peak-gene correlation chunk ", chunk_i, "/", length(pair_chunks))
  }
}

peak_gene_links <- data.table::rbindlist(peak_gene_result, use.names = TRUE)
peak_gene_links[, FDR := stats::p.adjust(pvalue, method = "BH")]
peak_gene_links_sig <- peak_gene_links[
  abs(cor) >= min_peak_gene_cor & FDR <= max_peak_gene_fdr
]
peak_gene_links_positive <- peak_gene_links[
  cor >= min_peak_gene_cor & FDR <= max_peak_gene_fdr
]

if (!nrow(peak_gene_links_positive)) {
  stop("No positive peak-gene link passed correlation and FDR thresholds.")
}

data.table::fwrite(
  peak_gene_links,
  file.path(output_dir, "02_peak_gene_all.tsv.gz"),
  sep = "\t"
)
data.table::fwrite(
  peak_gene_links_sig,
  file.path(output_dir, "02_peak_gene_significant.tsv"),
  sep = "\t"
)
data.table::fwrite(
  peak_gene_links_positive,
  file.path(output_dir, "02_peak_gene_positive.tsv"),
  sep = "\t"
)

message(
  "Step 5 complete: all links = ", nrow(peak_gene_links),
  "; positive significant links = ", nrow(peak_gene_links_positive)
)

# ----------------------------------------------------------------------------
# Step 6. Construct the full-genome TF-peak-gene global GRN
# ----------------------------------------------------------------------------
global_grn_raw <- merge(
  peak_gene_links_positive,
  tf_peak,
  by.x = "peak",
  by.y = "peak_id",
  allow.cartesian = TRUE
)
global_grn_raw <- global_grn_raw[TF != gene]

# Collapse multiple JASPAR motifs for the same TF at the same regulatory edge.
grn <- global_grn_raw[, .(
  motif_id = paste(sort(unique(motif_id)), collapse = ";"),
  distance = distance[1L],
  PeakGeneCor = cor[1L],
  PeakGeneP = pvalue[1L],
  PeakGeneFDR = FDR[1L]
), by = .(TF, gene, peak_id = peak)]

grn[, PeakGeneStrength :=
  PeakGeneCor *
    (pmin(-log10(pmax(PeakGeneFDR, .Machine$double.xmin)), 50) / 50)]

tf_target_counts <- unique(grn[, .(TF, gene)])[, .(N_targets = .N), by = TF]
valid_tf <- tf_target_counts[N_targets >= min_targets_per_tf, TF]
grn <- grn[TF %in% valid_tf]
if (!nrow(grn)) stop("No TF retained at least the configured number of targets.")

data.table::setcolorder(grn, c(
  "TF", "gene", "peak_id", "motif_id", "distance",
  "PeakGeneCor", "PeakGeneP", "PeakGeneFDR", "PeakGeneStrength"
))
data.table::fwrite(
  grn,
  file.path(output_dir, "03_TF_peak_gene_global_GRN.tsv"),
  sep = "\t"
)

message(
  "Step 6 complete: global GRN edges = ", nrow(grn),
  "; TFs = ", uniqueN(grn$TF),
  "; targets = ", uniqueN(grn$gene),
  "; chromosomes = ", uniqueN(sub("-.*$", "", grn$peak_id))
)

# ----------------------------------------------------------------------------
# Step 7. Peak cell-type accessibility and cell-type-specific sub-GRN
# ----------------------------------------------------------------------------
celltypes <- levels(celltype_groups)
peak_accessibility_list <- vector("list", length(celltypes))
names(peak_accessibility_list) <- celltypes

for (celltype_now in celltypes) {
  cells_now <- colnames(atac_counts)[celltype_vec == celltype_now]
  peak_fraction <- Matrix::rowMeans(
    atac_counts[, cells_now, drop = FALSE] > 0
  )
  peak_accessibility_list[[celltype_now]] <- data.table::data.table(
    peak_id = rownames(atac_counts),
    Celltype = celltype_now,
    PeakAccessibility = as.numeric(peak_fraction)
  )
}

peak_accessibility <- data.table::rbindlist(peak_accessibility_list)
peak_accessibility[, PeakNorm := sqrt(sum(PeakAccessibility^2)), by = peak_id]
peak_accessibility[
  !is.finite(PeakNorm) | PeakNorm == 0,
  PeakNorm := 1
]
peak_accessibility[, Peak_CTS := PeakAccessibility / PeakNorm]
peak_accessibility[, PeakNorm := NULL]

data.table::fwrite(
  peak_accessibility,
  file.path(output_dir, "04_peak_celltype_accessibility.tsv.gz"),
  sep = "\t"
)

celltype_grn <- merge(
  grn,
  target_scores,
  by = "gene",
  allow.cartesian = TRUE
)
data.table::setnames(celltype_grn, "Gene_CTS", "Target_CTS")

tf_cts <- data.table::copy(target_scores)
data.table::setnames(tf_cts, c("gene", "Gene_CTS"), c("TF", "TF_CTS"))
celltype_grn <- merge(
  celltype_grn,
  tf_cts,
  by = c("TF", "Celltype"),
  all.x = TRUE
)
celltype_grn[is.na(TF_CTS), TF_CTS := 0]

celltype_grn <- merge(
  celltype_grn,
  peak_accessibility,
  by = c("peak_id", "Celltype"),
  all.x = TRUE
)
celltype_grn[is.na(PeakAccessibility), PeakAccessibility := 0]
celltype_grn[is.na(Peak_CTS), Peak_CTS := 0]

celltype_grn[, CellEdgeScore :=
  PeakGeneStrength * Peak_CTS *
    sqrt(pmax(TF_CTS, 0) * pmax(Target_CTS, 0))]
celltype_grn <- celltype_grn[
  is.finite(CellEdgeScore) & CellEdgeScore > 0
]
celltype_grn[, EdgeQuantile := stats::quantile(
  CellEdgeScore,
  probs = cell_edge_quantile,
  na.rm = TRUE
), by = .(TF, Celltype)]

celltype_subgrn <- celltype_grn[CellEdgeScore >= EdgeQuantile]
subgrn_target_count <- unique(
  celltype_subgrn[, .(TF, Celltype, gene)]
)[, .(N_targets = .N), by = .(TF, Celltype)]
valid_eregulons <- subgrn_target_count[N_targets >= min_targets_per_tf]
celltype_subgrn <- merge(
  celltype_subgrn,
  valid_eregulons[, .(TF, Celltype)],
  by = c("TF", "Celltype")
)

data.table::fwrite(
  celltype_grn,
  file.path(output_dir, "05_celltype_GRN.tsv.gz"),
  sep = "\t"
)
data.table::fwrite(
  celltype_subgrn,
  file.path(output_dir, "06_celltype_subGRN.tsv.gz"),
  sep = "\t"
)

message(
  "Step 7 complete: cell-type sub-GRN edges = ", nrow(celltype_subgrn),
  "; eRegulons = ", uniqueN(celltype_subgrn[, .(TF, Celltype)])
)

# ----------------------------------------------------------------------------
# Step 8. Read somatic mutations and calculate peak MutationRelevance
#
# Mutation relevance is determined only by cohort recurrence:
#   VariantMutationRelevance = log(1 + n_samples)
# Functional annotations are retained for description, but no manually
# assigned functional-impact weight is used in MGRS.
# ----------------------------------------------------------------------------
mutation_required <- c(
  "Chr", "Start", "End", "Ref", "Alt", "Func.refGene",
  "Gene.refGene", "ExonicFunc.refGene", "Tumor_Sample_Barcode"
)
mutation_header <- names(data.table::fread(mutation_file, nrows = 0L, sep = "\t"))
mutation_missing <- setdiff(mutation_required, mutation_header)
if (length(mutation_missing)) {
  stop("Mutation table lacks columns: ", paste(mutation_missing, collapse = ", "))
}

mut <- data.table::fread(
  mutation_file,
  sep = "\t",
  select = mutation_required,
  na.strings = c("", ".", "NA"),
  showProgress = TRUE
)
data.table::setnames(mut, mutation_required, c(
  "chr", "start", "end", "ref", "alt", "func",
  "gene_anno", "exonic_func", "sample_id"
))

mut[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
mut[chr == "chr23", chr := "chrX"]
mut[chr == "chr24", chr := "chrY"]
mut[chr %in% c("chrM", "chrMT"), chr := "chrM"]
mut <- mut[
  grepl("^chr(?:[0-9]+|X|Y|M)$", chr) &
    is.finite(start) & is.finite(end) & start >= 1L & end >= start &
    !is.na(sample_id)
]
mut[, variant_id := paste(chr, start, ref, alt, sep = ":")]
mut <- unique(mut, by = c("variant_id", "sample_id"))

mut[, impact_class := "other"]
mut[grepl("intergenic", tolower(fifelse(is.na(func), "", func))), impact_class := "intergenic"]
mut[grepl("intronic", tolower(fifelse(is.na(func), "", func))), impact_class := "intronic"]
mut[grepl("utr|upstream|downstream", tolower(fifelse(is.na(func), "", func))), impact_class := "regulatory"]
mut[
  grepl("splic", tolower(fifelse(is.na(func), "", func))) |
    grepl("splic", tolower(fifelse(is.na(exonic_func), "", exonic_func))),
  impact_class := "splice"
]
mut[grepl("synonymous", tolower(fifelse(is.na(exonic_func), "", exonic_func))), impact_class := "synonymous"]
mut[grepl("nonsynonymous|nonframeshift", tolower(fifelse(is.na(exonic_func), "", exonic_func))), impact_class := "missense"]
mut[grepl("stopgain|stoploss|frameshift", tolower(fifelse(is.na(exonic_func), "", exonic_func))), impact_class := "high_impact"]

variants <- mut[, .(
  chr = chr[1L],
  start = start[1L],
  end = end[1L],
  ref = ref[1L],
  alt = alt[1L],
  n_samples = uniqueN(sample_id),
  impact_class = impact_class[1L],
  func = func[1L],
  exonic_func = exonic_func[1L],
  gene_anno = gene_anno[1L]
), by = variant_id]
variants <- variants[n_samples >= min_recurrent_samples]
variants[, VariantMutationRelevance := log1p(n_samples)]

peak_match <- stringr::str_match(
  peak_ids,
  "^(chr(?:[0-9]+|X|Y|M|MT))[:-]([0-9]+)[:-]([0-9]+)$"
)
valid_peak <- complete.cases(peak_match[, 2:4])
peak_gr_standard <- GenomicRanges::GRanges(
  seqnames = peak_match[valid_peak, 2],
  ranges = IRanges::IRanges(
    start = as.integer(peak_match[valid_peak, 3]),
    end = as.integer(peak_match[valid_peak, 4])
  ),
  peak_id = peak_ids[valid_peak]
)
peak_gr_query <- peak_gr_standard
GenomicRanges::start(peak_gr_query) <- pmax(
  1L,
  GenomicRanges::start(peak_gr_query) - mutation_peak_buffer
)
GenomicRanges::end(peak_gr_query) <-
  GenomicRanges::end(peak_gr_query) + mutation_peak_buffer

mutation_gr <- GenomicRanges::GRanges(
  seqnames = variants$chr,
  ranges = IRanges::IRanges(variants$start, variants$end),
  variant_id = variants$variant_id
)
mutation_peak_hits <- GenomicRanges::findOverlaps(
  mutation_gr,
  peak_gr_query,
  ignore.strand = TRUE
)
if (!length(mutation_peak_hits)) {
  stop("No mutation overlaps an ATAC peak. Check genome build and chromosome names.")
}

mutation_to_peak <- data.table::data.table(
  variant_id = S4Vectors::mcols(mutation_gr)$variant_id[
    S4Vectors::queryHits(mutation_peak_hits)
  ],
  peak_id = S4Vectors::mcols(peak_gr_query)$peak_id[
    S4Vectors::subjectHits(mutation_peak_hits)
  ]
)
mutation_to_peak <- merge(
  mutation_to_peak,
  variants,
  by = "variant_id",
  all.x = TRUE
)

peak_scores <- mutation_to_peak[, .(
  n_variants = uniqueN(variant_id),
  n_mutated_samples = max(n_samples),
  mutation_peak_raw = sum(VariantMutationRelevance, na.rm = TRUE)
), by = peak_id]

peak_cap <- as.numeric(stats::quantile(
  peak_scores$mutation_peak_raw,
  0.99,
  na.rm = TRUE
))
peak_scores[, mutation_peak_raw_capped := pmin(mutation_peak_raw, peak_cap)]

# Normalize relative to zero, not relative to the smallest mutated peak.
# Therefore every genuinely mutated peak retains a score greater than zero.
peak_score_denominator <- log1p(max(
  peak_scores$mutation_peak_raw_capped,
  na.rm = TRUE
))
if (is.finite(peak_score_denominator) && peak_score_denominator > 0) {
  peak_scores[, MutationPeakScore :=
    log1p(mutation_peak_raw_capped) / peak_score_denominator]
} else {
  peak_scores[, MutationPeakScore := as.numeric(mutation_peak_raw_capped > 0)]
}

# Formula-oriented aliases; MutationPeakScore is kept so Step 9 is unchanged.
peak_scores[, `:=`(
  MutationRelevanceRaw = mutation_peak_raw,
  MutationRelevanceCapped = mutation_peak_raw_capped,
  MutationRelevance = MutationPeakScore
)]

data.table::fwrite(
  mutation_to_peak,
  file.path(output_dir, "07_mutation_peak.tsv.gz"),
  sep = "\t"
)
data.table::fwrite(
  peak_scores,
  file.path(output_dir, "08_peak_mutation_scores.tsv"),
  sep = "\t"
)

message(
  "Step 8 complete: variants overlapping peaks = ",
  uniqueN(mutation_to_peak$variant_id),
  "; mutated peaks = ", nrow(peak_scores)
)

# ----------------------------------------------------------------------------
# Step 9. Mutation-weighted regulatory edges and target MGRS
# ----------------------------------------------------------------------------
mutation_edges <- merge(
  celltype_subgrn,
  peak_scores,
  by = "peak_id",
  all.x = TRUE
)
mutation_edges[is.na(MutationPeakScore), MutationPeakScore := 0]
mutation_edges[is.na(n_variants), n_variants := 0L]
mutation_edges[is.na(n_mutated_samples), n_mutated_samples := 0L]
mutation_edges[is.na(mutation_peak_raw), mutation_peak_raw := 0]
mutation_edges[is.na(mutation_peak_raw_capped), mutation_peak_raw_capped := 0]

mutation_edges[, MGRS_edge := MutationPeakScore * CellEdgeScore]

target_mutation_status <- mutation_edges[, .(
  Target_MGRS = max(MGRS_edge, na.rm = TRUE)
), by = .(TF, Celltype, gene)]
target_mutation_status[!is.finite(Target_MGRS), Target_MGRS := 0]

data.table::fwrite(
  mutation_edges,
  file.path(output_dir, "09_TF_peak_gene_mutation_edges.tsv.gz"),
  sep = "\t"
)
data.table::fwrite(
  target_mutation_status,
  file.path(output_dir, "10_TF_celltype_target_MGRS.tsv"),
  sep = "\t"
)

message("Step 9 complete: mutation regulatory edges = ", nrow(mutation_edges))

# ----------------------------------------------------------------------------
# Step 10. Aggregate eRegulon CTS, MGRS and MRS
# ----------------------------------------------------------------------------
eregulon_target_info <- celltype_subgrn[, .(
  TF_CTS = TF_CTS[1L],
  Target_CTS = max(Target_CTS, na.rm = TRUE),
  MeanPeakCTS = mean(Peak_CTS, na.rm = TRUE),
  MaxPeakCTS = max(Peak_CTS, na.rm = TRUE),
  MaxCellEdgeScore = max(CellEdgeScore, na.rm = TRUE)
), by = .(TF, Celltype, gene)]

eregulon_target_info <- merge(
  eregulon_target_info,
  target_mutation_status,
  by = c("TF", "Celltype", "gene"),
  all.x = TRUE
)
eregulon_target_info[is.na(Target_MGRS), Target_MGRS := 0]

score_table <- eregulon_target_info[, {
  mgrs_values <- Target_MGRS
  mgrs_values[!is.finite(mgrs_values)] <- 0
  mgrs_sorted <- sort(mgrs_values, decreasing = TRUE)
  mgrs_top <- mean(head(mgrs_sorted, min(top_n, length(mgrs_sorted))))
  if (!is.finite(mgrs_top)) mgrs_top <- 0

  .(
    N_targets = uniqueN(gene),
    CTS_raw = TF_CTS[1L] + theta * sum(Target_CTS, na.rm = TRUE) /
      sqrt(uniqueN(gene)),
    N_mutation_linked_targets = sum(mgrs_values > 0),
    MutationLinkedFraction = mean(mgrs_values > 0),
    MGRS_raw = mgrs_top + theta * sum(mgrs_values, na.rm = TRUE) /
      sqrt(uniqueN(gene))
  )
}, by = .(TF, Celltype)]

score_table[, CTS_mean_celltype := mean(CTS_raw, na.rm = TRUE), by = Celltype]
score_table[, CTS_sd_celltype := stats::sd(CTS_raw, na.rm = TRUE), by = Celltype]
score_table[, MGRS_mean_celltype := mean(MGRS_raw, na.rm = TRUE), by = Celltype]
score_table[, MGRS_sd_celltype := stats::sd(MGRS_raw, na.rm = TRUE), by = Celltype]
score_table[!is.finite(CTS_sd_celltype) | CTS_sd_celltype == 0, CTS_sd_celltype := 1]
score_table[!is.finite(MGRS_sd_celltype) | MGRS_sd_celltype == 0, MGRS_sd_celltype := 1]
score_table[, CTS_z := (CTS_raw - CTS_mean_celltype) / CTS_sd_celltype]
score_table[, MGRS_z := (MGRS_raw - MGRS_mean_celltype) / MGRS_sd_celltype]
score_table[, Penalty := abs(CTS_z - MGRS_z) / sqrt(2)]
score_table[, MRS_raw := CTS_z + MGRS_z - Penalty]
score_table[, c(
  "CTS_mean_celltype", "CTS_sd_celltype",
  "MGRS_mean_celltype", "MGRS_sd_celltype"
) := NULL]

data.table::fwrite(
  score_table,
  file.path(output_dir, "11_eRegulon_CTS_MGRS_MRS.tsv"),
  sep = "\t"
)

message("Step 10 complete: eRegulons = ", nrow(score_table))

# ----------------------------------------------------------------------------
# Step 11. Matched-control Monte Carlo and final statistics
# ----------------------------------------------------------------------------
eregulon_metadata <- eregulon_target_info[, .(
  N_targets = uniqueN(gene),
  TF_CTS = TF_CTS[1L],
  MeanTargetCTS = mean(Target_CTS, na.rm = TRUE),
  SDTargetCTS = stats::sd(Target_CTS, na.rm = TRUE),
  MeanPeakCTS = mean(MeanPeakCTS, na.rm = TRUE)
), by = .(TF, Celltype)]
eregulon_metadata[!is.finite(SDTargetCTS), SDTargetCTS := 0]

celltype_scale <- score_table[, .(
  CTS_mean = mean(CTS_raw, na.rm = TRUE),
  CTS_sd = stats::sd(CTS_raw, na.rm = TRUE),
  MGRS_mean = mean(MGRS_raw, na.rm = TRUE),
  MGRS_sd = stats::sd(MGRS_raw, na.rm = TRUE)
), by = Celltype]
celltype_scale[!is.finite(CTS_sd) | CTS_sd == 0, CTS_sd := 1]
celltype_scale[!is.finite(MGRS_sd) | MGRS_sd == 0, MGRS_sd := 1]

eregulon_target_info[, RegulonID := paste(TF, Celltype, sep = "___")]
target_list <- split(eregulon_target_info, eregulon_target_info$RegulonID)

set.seed(seed)
mc_results <- vector("list", nrow(score_table))

for (obs_i in seq_len(nrow(score_table))) {
  tf_now <- score_table$TF[obs_i]
  celltype_now <- score_table$Celltype[obs_i]
  n_targets_now <- score_table$N_targets[obs_i]
  observed_cts_raw <- score_table$CTS_raw[obs_i]
  observed_mgrs_raw <- score_table$MGRS_raw[obs_i]
  observed_mrs <- score_table$MRS_raw[obs_i]

  obs_meta <- eregulon_metadata[TF == tf_now & Celltype == celltype_now]
  if (nrow(obs_meta) != 1L) next

  candidates <- data.table::copy(
    eregulon_metadata[Celltype == celltype_now & TF != tf_now]
  )
  if (!nrow(candidates)) next

  size_lower <- max(5L, floor(n_targets_now * 0.70))
  size_upper <- ceiling(n_targets_now * 1.30)
  candidates_size <- candidates[
    N_targets >= size_lower & N_targets <= size_upper
  ]
  if (nrow(candidates_size) >= min_control_pool) {
    candidates <- candidates_size
  }

  size_sd <- stats::sd(log1p(candidates$N_targets), na.rm = TRUE)
  tfcts_sd <- stats::sd(candidates$TF_CTS, na.rm = TRUE)
  targetcts_sd <- stats::sd(candidates$MeanTargetCTS, na.rm = TRUE)
  peakcts_sd <- stats::sd(candidates$MeanPeakCTS, na.rm = TRUE)
  if (!is.finite(size_sd) || size_sd == 0) size_sd <- 1
  if (!is.finite(tfcts_sd) || tfcts_sd == 0) tfcts_sd <- 1
  if (!is.finite(targetcts_sd) || targetcts_sd == 0) targetcts_sd <- 1
  if (!is.finite(peakcts_sd) || peakcts_sd == 0) peakcts_sd <- 1

  candidates[, MatchDistance := sqrt(
    ((log1p(N_targets) - log1p(obs_meta$N_targets)) / size_sd)^2 +
      ((TF_CTS - obs_meta$TF_CTS) / tfcts_sd)^2 +
      ((MeanTargetCTS - obs_meta$MeanTargetCTS) / targetcts_sd)^2 +
      ((MeanPeakCTS - obs_meta$MeanPeakCTS) / peakcts_sd)^2
  )]
  candidates <- candidates[is.finite(MatchDistance)]
  data.table::setorder(candidates, MatchDistance)
  candidates <- head(candidates, min(control_pool_n, nrow(candidates)))
  control_tfs <- candidates$TF
  if (!length(control_tfs)) next

  scale_now <- celltype_scale[Celltype == celltype_now]
  null_cts_raw <- rep(NA_real_, perm_n)
  null_mgrs_raw <- rep(NA_real_, perm_n)
  null_mrs <- rep(NA_real_, perm_n)

  for (perm_i in seq_len(perm_n)) {
    random_tf <- sample(control_tfs, 1L)
    random_id <- paste(random_tf, celltype_now, sep = "___")
    random_targets <- target_list[[random_id]]
    if (is.null(random_targets) || !nrow(random_targets)) next

    sampled_rows <- sample(
      seq_len(nrow(random_targets)),
      size = n_targets_now,
      replace = nrow(random_targets) < n_targets_now
    )
    random_targets_now <- random_targets[sampled_rows]

    random_cts_raw <- random_targets_now$TF_CTS[1L] +
      theta * sum(random_targets_now$Target_CTS, na.rm = TRUE) /
      sqrt(n_targets_now)

    random_mgrs_values <- random_targets_now$Target_MGRS
    random_mgrs_values[!is.finite(random_mgrs_values)] <- 0
    random_mgrs_sorted <- sort(random_mgrs_values, decreasing = TRUE)
    random_mgrs_top <- mean(head(
      random_mgrs_sorted,
      min(top_n, length(random_mgrs_sorted))
    ))
    if (!is.finite(random_mgrs_top)) random_mgrs_top <- 0
    random_mgrs_raw <- random_mgrs_top +
      theta * sum(random_mgrs_values, na.rm = TRUE) / sqrt(n_targets_now)

    random_cts_z <- (random_cts_raw - scale_now$CTS_mean) / scale_now$CTS_sd
    random_mgrs_z <- (random_mgrs_raw - scale_now$MGRS_mean) / scale_now$MGRS_sd
    random_penalty <- abs(random_cts_z - random_mgrs_z) / sqrt(2)

    null_cts_raw[perm_i] <- random_cts_raw
    null_mgrs_raw[perm_i] <- random_mgrs_raw
    null_mrs[perm_i] <- random_cts_z + random_mgrs_z - random_penalty
  }

  valid_null <- is.finite(null_cts_raw) &
    is.finite(null_mgrs_raw) & is.finite(null_mrs)
  null_cts_raw <- null_cts_raw[valid_null]
  null_mgrs_raw <- null_mgrs_raw[valid_null]
  null_mrs <- null_mrs[valid_null]
  n_valid <- length(null_mrs)
  if (n_valid < 10L) next

  cts_p <- (1 + sum(null_cts_raw >= observed_cts_raw)) / (1 + n_valid)
  mgrs_p <- (1 + sum(null_mgrs_raw >= observed_mgrs_raw)) / (1 + n_valid)
  mrs_p <- (1 + sum(null_mrs >= observed_mrs)) / (1 + n_valid)

  cts_null_z <- if (stats::sd(null_cts_raw) > 0) {
    (observed_cts_raw - mean(null_cts_raw)) / stats::sd(null_cts_raw)
  } else NA_real_
  mgrs_null_z <- if (stats::sd(null_mgrs_raw) > 0) {
    (observed_mgrs_raw - mean(null_mgrs_raw)) / stats::sd(null_mgrs_raw)
  } else NA_real_
  mrs_null_z <- if (stats::sd(null_mrs) > 0) {
    (observed_mrs - mean(null_mrs)) / stats::sd(null_mrs)
  } else NA_real_

  mc_results[[obs_i]] <- data.table::data.table(
    TF = tf_now,
    Celltype = celltype_now,
    N_controls = length(control_tfs),
    N_valid_permutations = n_valid,
    MeanMatchDistance = mean(candidates$MatchDistance),
    MinMatchDistance = min(candidates$MatchDistance),
    CTS_P = cts_p,
    MGRS_P = mgrs_p,
    MRS_P = mrs_p,
    CTS_null_z = cts_null_z,
    MGRS_null_z = mgrs_null_z,
    MRS_null_z = mrs_null_z,
    Null_CTS_mean = mean(null_cts_raw),
    Null_CTS_sd = stats::sd(null_cts_raw),
    Null_MGRS_mean = mean(null_mgrs_raw),
    Null_MGRS_sd = stats::sd(null_mgrs_raw),
    Null_MRS_mean = mean(null_mrs),
    Null_MRS_sd = stats::sd(null_mrs)
  )

  if (obs_i %% 100L == 0L || obs_i == nrow(score_table)) {
    message("Monte Carlo eRegulon ", obs_i, "/", nrow(score_table))
  }
}

mc_results_dt <- data.table::rbindlist(mc_results, fill = TRUE)
final_results <- merge(
  score_table,
  mc_results_dt,
  by = c("TF", "Celltype"),
  all.x = TRUE
)

final_results[, CTS_FDR := stats::p.adjust(CTS_P, method = "BH")]
final_results[, MGRS_FDR := stats::p.adjust(MGRS_P, method = "BH")]
final_results[, MRS_FDR := stats::p.adjust(MRS_P, method = "BH")]
final_results[, MRS_FDR_Celltype :=
  stats::p.adjust(MRS_P, method = "BH"), by = Celltype]
final_results[, Significant :=
  !is.na(MRS_P) & !is.na(MRS_FDR) &
    MRS_P < nominal_p & MRS_FDR < fdr_cutoff]
final_results[, ConcordantPositive := CTS_z > 0 & MGRS_z > 0]
final_results[, HighConfidence := Significant & ConcordantPositive]

# Exploratory labels are never substitutes for FDR-significant results.
final_results[, Candidate := !is.na(MRS_P) & MRS_P < nominal_p]
final_results[, StrongCandidate := Candidate & ConcordantPositive]
data.table::setorder(final_results, MRS_FDR, MRS_P, -MRS_raw)

data.table::fwrite(
  mc_results_dt,
  file.path(output_dir, "12_matched_MonteCarlo_statistics.tsv"),
  sep = "\t"
)
data.table::fwrite(
  final_results,
  file.path(output_dir, "12_final_results.tsv"),
  sep = "\t"
)
data.table::fwrite(
  final_results[HighConfidence == TRUE],
  file.path(output_dir, "12_final_high_confidence.tsv"),
  sep = "\t"
)

saveRDS(
  list(
    settings = list(
      variable_genes = variable_genes,
      peak_gene_distance = peak_gene_distance,
      min_peak_gene_cor = min_peak_gene_cor,
      max_peak_gene_fdr = max_peak_gene_fdr,
      min_targets_per_tf = min_targets_per_tf,
      cell_edge_quantile = cell_edge_quantile,
      theta = theta,
      top_n = top_n,
      perm_n = perm_n,
      seed = seed,
      control_pool_n = control_pool_n
    ),
    grn = grn,
    peak_accessibility = peak_accessibility,
    celltype_subgrn = celltype_subgrn,
    mutation_to_peak = mutation_to_peak,
    peak_scores = peak_scores,
    mutation_edges = mutation_edges,
    target_mutation_status = target_mutation_status,
    score_table = score_table,
    monte_carlo = mc_results_dt,
    final_results = final_results
  ),
  file.path(output_dir, "Final_GRN_complete_objects.rds")
)

run_summary <- data.table::data.table(
  Metric = c(
    "paired_cells", "ATAC_peaks", "motif_TF_peak_rows",
    "peak_gene_tests", "positive_significant_peak_gene_links",
    "global_GRN_edges", "global_GRN_TFs", "global_GRN_targets",
    "celltype_subGRN_edges", "eRegulons", "mutated_peaks",
    "final_tests", "significant_MRS", "high_confidence_MRS"
  ),
  Value = c(
    ncol(rna_data), nrow(atac_counts), nrow(tf_peak),
    nrow(peak_gene_links), nrow(peak_gene_links_positive),
    nrow(grn), uniqueN(grn$TF), uniqueN(grn$gene),
    nrow(celltype_subgrn), nrow(score_table), nrow(peak_scores),
    nrow(final_results), sum(final_results$Significant, na.rm = TRUE),
    sum(final_results$HighConfidence, na.rm = TRUE)
  )
)
data.table::fwrite(
  run_summary,
  file.path(output_dir, "00_run_summary.tsv"),
  sep = "\t"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "00_sessionInfo.txt")
)

message("Step 11 complete: final eRegulon tests = ", nrow(final_results))
message("Significant MRS = ", sum(final_results$Significant, na.rm = TRUE))
message("High-confidence MRS = ", sum(final_results$HighConfidence, na.rm = TRUE))
message("All outputs written to: ", output_dir)
print(head(final_results, 20L))
