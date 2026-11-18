# ================================
# 0. Load R packages
# ================================
library(Seurat)
library(Signac)

library(TFBSTools)
library(JASPAR2020)

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

library(GenomeInfoDb)
library(GenomicRanges)
library(IRanges)
library(Biostrings)

library(motifmatchr)


# ================================
# 1. Set the ATAC assay
# ================================
DefaultAssay(lung_ATAC) <- "ATAC"

lung_ATAC


# ================================
# 2. Load the hg38 genome
# ================================
hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

hg38


# ================================
# 3. Check whether BSgenome methods work correctly
# ================================
methods::selectMethod("seqnames", "BSgenome")
methods::selectMethod("seqinfo", "BSgenome")

# Inspect hg38 chromosome information
GenomeInfoDb::seqinfo(hg38)


# ================================
# 4. Extract peaks from lung_ATAC
# ================================
# Because granges(ChromatinAssay) previously had a method conflict in this environment,
# extracting directly from the ranges slot is the most reliable approach
peak.granges <- lung_ATAC[["ATAC"]]@ranges

class(peak.granges)
length(peak.granges)

head(peak.granges)


# ================================
# 5. Check peak chromosome names
# ================================
head(as.character(seqnames(peak.granges)))

unique(as.character(seqnames(peak.granges)))


# ================================
# 6. Check whether all peaks are present in hg38
# ================================
peak_chr <- unique(as.character(seqnames(peak.granges)))

hg38_chr <- seqnames(GenomeInfoDb::seqinfo(hg38))

# Chromosomes absent from hg38
setdiff(
  peak_chr,
  hg38_chr
)


# ================================
# 7. Test extraction of peak sequences from hg38
# ================================
test.seq <- BSgenome::getSeq(
  hg38,
  peak.granges[1:10]
)

test.seq


# ================================
# 8. Obtain JASPAR human CORE motifs
# ================================
pfm <- TFBSTools::getMatrixSet(
  JASPAR2020::JASPAR2020,
  opts = list(
    species = 9606,
    collection = "CORE",
    all_versions = FALSE
  )
)

length(pfm)

# Inspect the first few motifs
pfm[1:3]


# ================================
# 9. Check whether motifmatchr works correctly
# ================================
packageVersion("motifmatchr")

# Test motif matching with a small set of peaks
test.motif <- motifmatchr::matchMotifs(
  pwms = pfm,
  subject = peak.granges[1:100],
  genome = hg38
)

test.motif


# ================================
# 10. Add motifs to lung_ATAC
# ================================
lung_ATAC <- AddMotifs(
  object = lung_ATAC,
  genome = hg38,
  pfm = pfm
)


# ================================
# 11. Check whether motifs were added successfully
# ================================
lung_ATAC

Motifs(lung_ATAC)


motif_obj <- Motifs(lung_ATAC)

slotNames(motif_obj)

motif_matrix <- motif_obj@data

dim(motif_matrix)

motif_matrix[1:5, 1:5]


# ================================
# 12. Save the results
# ================================
saveRDS(
  lung_ATAC,
  file = "E:/Lung/lung_ATAC_with_JASPAR_motifs.rds"
)
