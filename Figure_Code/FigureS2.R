# ============================================================
# 96-channel trinucleotide mutation spectrum
# hg38 / GRCh38
# ============================================================

library(maftools)
library(dplyr)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(GenomicRanges)
library(scales)

# ============================================================
# 1. Extract SNVs
# ============================================================

snv <- maf_lung@data %>%
  as.data.frame() %>%
  filter(
    Variant_Type == "SNP",
    nchar(Reference_Allele) == 1,
    nchar(Tumor_Seq_Allele2) == 1,
    Reference_Allele %in% c("A", "C", "G", "T"),
    Tumor_Seq_Allele2 %in% c("A", "C", "G", "T"),
    Reference_Allele != Tumor_Seq_Allele2
  ) %>%
  mutate(
    Chromosome = as.character(Chromosome),
    Start_Position = as.integer(Start_Position)
  )

nrow(snv)


# ============================================================
# 2. Standardize chromosome names
# ============================================================

snv$Chromosome <- ifelse(
  grepl("^chr", snv$Chromosome),
  snv$Chromosome,
  paste0("chr", snv$Chromosome)
)

# Common mitochondrial chromosome names
snv$Chromosome[
  snv$Chromosome %in% c("chrMT", "chrM")
] <- "chrM"

# Retain only chromosomes present in the reference genome
snv <- snv %>%
  filter(
    Chromosome %in%
      seqnames(BSgenome.Hsapiens.UCSC.hg38)
  )


# ============================================================
# 3. Construct +/-1-bp genomic ranges
# ============================================================

gr <- GRanges(
  seqnames = snv$Chromosome,
  ranges = IRanges(
    start = snv$Start_Position - 1,
    end   = snv$Start_Position + 1
  )
)

# Extract trinucleotide sequences
tri_seq <- as.character(
  getSeq(
    BSgenome.Hsapiens.UCSC.hg38,
    gr
  )
)

snv$trinucleotide <- toupper(tri_seq)

head(
  snv[, c(
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "trinucleotide"
  )]
)



# ============================================================
# 4. reverse complement function
# ============================================================

rev_comp <- function(x) {
  as.character(
    reverseComplement(
      DNAStringSet(x)
    )
  )
}



# ============================================================
# 5. Standardize to C/T as the reference base
# ============================================================

snv96 <- snv %>%
  mutate(
    ref = Reference_Allele,
    alt = Tumor_Seq_Allele2,
    context = trinucleotide
  )

# Use the reverse complement when the reference base is A/G
idx <- snv96$ref %in% c("A", "G")

snv96$context[idx] <-
  rev_comp(
    snv96$context[idx]
  )

# nucleotide complement
comp_base <- c(
  A = "T",
  T = "A",
  C = "G",
  G = "C"
)

snv96$ref[idx] <-
  comp_base[
    snv96$ref[idx]
  ]

snv96$alt[idx] <-
  comp_base[
    snv96$alt[idx]
  ]


# ============================================================
# 6. Substitution class
# ============================================================

snv96 <- snv96 %>%
  mutate(
    substitution = paste0(
      ref,
      ">",
      alt
    )
  ) %>%
  filter(
    substitution %in%
      c(
        "C>A",
        "C>G",
        "C>T",
        "T>A",
        "T>C",
        "T>G"
      )
  )

table(
  snv96$substitution
)


# ============================================================
# 7. Flanking context
# ============================================================

snv96 <- snv96 %>%
  mutate(
    left_base =
      substr(
        context,
        1,
        1
      ),
    
    right_base =
      substr(
        context,
        3,
        3
      ),
    
    flank = paste0(
      left_base,
      ".",
      right_base
    )
  )




# ============================================================
# 8. Define 16-context order
# ============================================================

bases <- c(
  "A",
  "C",
  "G",
  "T"
)

context_order <- c(
  "A.A", "A.C", "A.G", "A.T",
  "C.A", "C.C", "C.G", "C.T",
  "G.A", "G.C", "G.G", "G.T",
  "T.A", "T.C", "T.G", "T.T"
)

substitution_order <- c(
  "C>A",
  "C>G",
  "C>T",
  "T>A",
  "T>C",
  "T>G"
)



# ============================================================
# 9. Count 96 channels
# ============================================================

sbs96 <- snv96 %>%
  count(
    substitution,
    flank,
    name = "Count"
  ) %>%
  tidyr::complete(
    substitution =
      substitution_order,
    
    flank =
      context_order,
    
    fill =
      list(
        Count = 0
      )
  )

sbs96 <- sbs96 %>%
  mutate(
    Relative_contribution =
      Count /
      sum(Count)
  )



# ============================================================
# 10. Plot
# ============================================================

p96 <- ggplot(
  sbs96,
  aes(
    x = flank,
    y = Relative_contribution,
    fill = substitution
  )
) +
  
  geom_col(
    width = 0.88,
    color = "black",
    linewidth = 0.35
  ) +
  
  facet_grid(
    . ~ substitution,
    scales = "free_x",
    space = "free_x",
    switch = "x"
  ) +
  
  scale_fill_manual(
    values = c(
      "C>A" = "#43B5D8",
      "C>G" = "#222222",
      "C>T" = "#D73027",
      "T>A" = "#C9C9C9",
      "T>C" = "#9AC44D",
      "T>G" = "#E6C9C9"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.06
      )
    )
  ) +
  
  labs(
    x = NULL,
    y = "Relative contribution"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    
    legend.position = "none",
    
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 7,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 11
    ),
    
    axis.line = element_line(
      color = "black",
      linewidth = 0.5
    ),
    
    axis.ticks = element_line(
      color = "black"
    ),
    
    strip.background = element_rect(
      fill = "#D0D0D0",
      color = "black",
      linewidth = 0.6
    ),
    
    strip.text = element_text(
      size = 10,
      color = "black",
      margin = margin(
        3,
        3,
        3,
        3
      )
    ),
    
    panel.spacing.x = unit(
      0.08,
      "cm"
    ),
    
    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )

p96
ggsave(
  "E:/Lung/EGA/PDF/SBS96_trinucleotide_spectrum.pdf",
  p96,
  width = 11,
  height = 4.2
)






# ============================================================
# Co-occurrence and mutual exclusivity
# Top15 mutated genes + 5 well-known lung cancer genes
# ============================================================

library(maftools)
library(dplyr)
library(ggplot2)
library(grid)


# ============================================================
# 1. Obtain mutation frequencies for all genes
# ============================================================

gene_tbl_all <- getGeneSummary(maf_lung) %>%
  as.data.frame() %>%
  dplyr::select(
    Hugo_Symbol,
    MutatedSamples
  ) %>%
  mutate(
    Hugo_Symbol = trimws(as.character(Hugo_Symbol))
  ) %>%
  arrange(
    desc(MutatedSamples)
  )


# ============================================================
# 2. Top 15 frequently mutated genes
# ============================================================

top15_genes <- gene_tbl_all %>%
  slice_head(n = 15) %>%
  pull(Hugo_Symbol)

top15_genes


# ============================================================
# 3. Five canonical lung cancer-related genes
# ============================================================

famous_genes <- c(
  "TP53",
  "EGFR",
  "KRAS",
  "STK11",
  "KEAP1"
)

# Retain only genes present in the data
famous_genes <- famous_genes[
  famous_genes %in% gene_tbl_all$Hugo_Symbol
]

famous_genes


# ============================================================
# 4. Combine genes
# Use unique to prevent duplication when a canonical gene is already in the top 15
# ============================================================

genes <- unique(
  c(
    top15_genes,
    famous_genes
  )
)

genes
length(genes)


# ============================================================
# 5. Inspect the number of mutated samples for the final gene set
# ============================================================

gene_tbl <- gene_tbl_all %>%
  filter(
    Hugo_Symbol %in% genes
  ) %>%
  mutate(
    gene_order = match(
      Hugo_Symbol,
      genes
    )
  ) %>%
  arrange(
    gene_order
  )

gene_tbl


# ============================================================
# 6. Obtain all samples
# ============================================================

samples <- as.character(
  getSampleSummary(maf_lung)$Tumor_Sample_Barcode
)

length(samples)


# ============================================================
# 7. Extract mutation data for the final gene set
# ============================================================

mut_df <- maf_lung@data %>%
  as.data.frame() %>%
  transmute(
    Hugo_Symbol =
      trimws(
        as.character(Hugo_Symbol)
      ),
    
    Tumor_Sample_Barcode =
      as.character(
        Tumor_Sample_Barcode
      )
  ) %>%
  distinct() %>%
  filter(
    Hugo_Symbol %in% genes,
    Tumor_Sample_Barcode %in% samples
  )


# ============================================================
# 8. Create a complete gene-by-sample binary matrix
# ============================================================

mat <- matrix(
  0,
  nrow = length(genes),
  ncol = length(samples),
  dimnames = list(
    genes,
    samples
  )
)


gene_idx <- match(
  mut_df$Hugo_Symbol,
  rownames(mat)
)

sample_idx <- match(
  mut_df$Tumor_Sample_Barcode,
  colnames(mat)
)

valid_idx <-
  !is.na(gene_idx) &
  !is.na(sample_idx)


mat[
  cbind(
    gene_idx[valid_idx],
    sample_idx[valid_idx]
  )
] <- 1


# ============================================================
# 9. Check
# ============================================================

dim(mat)

rowSums(mat)


# ============================================================
# 10. Pairwise Fisher exact test
# ============================================================

pair_list <- list()


for (i in seq_along(genes)) {
  
  for (j in seq_along(genes)) {
    
    if (i <= j) next
    
    g1 <- genes[i]
    g2 <- genes[j]
    
    
    # --------------------------------------------------------
    # 2 x 2 table
    # --------------------------------------------------------
    
    a <- sum(
      mat[g1, ] == 1 &
        mat[g2, ] == 1
    )
    
    b <- sum(
      mat[g1, ] == 1 &
        mat[g2, ] == 0
    )
    
    c <- sum(
      mat[g1, ] == 0 &
        mat[g2, ] == 1
    )
    
    d <- sum(
      mat[g1, ] == 0 &
        mat[g2, ] == 0
    )
    
    
    tab <- matrix(
      c(
        a, b,
        c, d
      ),
      nrow = 2,
      byrow = TRUE
    )
    
    
    ft <- fisher.test(
      tab
    )
    
    
    OR <- unname(
      ft$estimate
    )
    
    
    # --------------------------------------------------------
    # Define the direction
    #
    # OR >= 1  => co-occurrence
    # OR < 1   => mutual exclusivity
    # --------------------------------------------------------
    
    if (
      is.na(OR) ||
      is.nan(OR) ||
      is.infinite(OR)
    ) {
      
      observed <-
        a / ncol(mat)
      
      expected <-
        mean(
          mat[g1, ]
        ) *
        mean(
          mat[g2, ]
        )
      
      direction <- ifelse(
        observed >= expected,
        1,
        -1
      )
      
    } else {
      
      direction <- ifelse(
        OR >= 1,
        1,
        -1
      )
    }
    
    
    pair_list[[length(pair_list) + 1]] <- data.frame(
      gene1 = g1,
      gene2 = g2,
      both_mutated = a,
      OR = OR,
      p = ft$p.value,
      direction = direction,
      signed_logp = direction * (-log10(ft$p.value + 1e-300))
    )
  }
}


pair_df <- bind_rows(
  pair_list
)


# ============================================================
# 11. Construct labels
# Gene [mutation sample count]
# ============================================================

label_map <- setNames(
  
  paste0(
    gene_tbl$Hugo_Symbol,
    " [",
    gene_tbl$MutatedSamples,
    "]"
  ),
  
  gene_tbl$Hugo_Symbol
)


label_map


# ============================================================
# 12. Add labels and significance indicators
# ============================================================

pair_df <- pair_df %>%
  mutate(
    
    gene1_lab =
      unname(
        label_map[gene1]
      ),
    
    gene2_lab =
      unname(
        label_map[gene2]
      ),
    
    
    # --------------------------------------------------------
    # significance
    # --------------------------------------------------------
    
    sig = case_when(
      
      p < 0.01 ~ "**",
      
      p < 0.05 ~ "*",
      
      TRUE ~ ""
    ),
    
    
    # --------------------------------------------------------
    # Limit the color scale to -3 to +3
    # --------------------------------------------------------
    
    plot_val =
      pmax(
        pmin(
          signed_logp,
          3
        ),
        -3
      )
  )


# ============================================================
# 13. Check missing values
# ============================================================

sum(
  is.na(
    pair_df$gene1_lab
  )
)

sum(
  is.na(
    pair_df$gene2_lab
  )
)

# All values should be zero


# ============================================================
# 14. Set the gene order
# ============================================================

label_order <-
  unname(
    label_map[
      genes
    ]
  )


pair_df <- pair_df %>%
  mutate(
    
    gene1_lab =
      factor(
        gene1_lab,
        levels =
          label_order
      ),
    
    gene2_lab =
      factor(
        gene2_lab,
        levels =
          rev(
            label_order
          )
      )
  )


# ============================================================
# 15. Plot
# ============================================================

p <- ggplot(
  
  pair_df,
  
  aes(
    x = gene1_lab,
    y = gene2_lab,
    fill = plot_val
  )
  
) +
  
  # ----------------------------------------------------------
# Tile
# ----------------------------------------------------------

geom_tile(
  
  color = "white",
  
  linewidth = 0.8,
  
  width = 0.96,
  
  height = 0.96
  
) +
  
  
  # ----------------------------------------------------------
# Significance
# ----------------------------------------------------------

geom_text(
  
  aes(
    label = sig
  ),
  
  fontface = "bold",
  
  size = 4,
  
  color = "black"
  
) +
  
  
  # ----------------------------------------------------------
# Color scale
# ----------------------------------------------------------

scale_fill_gradient2(
  
  # mutual exclusivity
  low = "#D6A03A",
  
  # neutral
  mid = "#F3F1EC",
  
  # co-occurrence
  high = "#5A9E8A",
  
  midpoint = 0,
  
  limits = c(
    -3,
    3
  ),
  
  breaks = c(
    -3,
    -2,
    -1,
    0,
    1,
    2,
    3
  ),
  
  labels = c(
    ">3",
    "2",
    "1",
    "0",
    "1",
    "2",
    ">3"
  ),
  
  name =
    expression(
      paste(
        "signed ",
        -log[10],
        "(P)"
      )
    )
) +
  
  
  coord_fixed() +
  
  
  # ----------------------------------------------------------
# Labels
# ----------------------------------------------------------

labs(
  
  x = NULL,
  
  y = NULL,
  
  title =
    "Co-occurrence and mutual exclusivity",
  
  caption =
    "Green: co-occurrence; orange: mutual exclusivity.  * P < 0.05;  ** P < 0.01"
  
) +
  
  
  # ----------------------------------------------------------
# Theme
# ----------------------------------------------------------

theme_classic(
  base_size = 12
) +
  
  
  theme(
    
    # axis
    axis.line =
      element_blank(),
    
    axis.ticks =
      element_blank(),
    
    
    axis.text.x =
      element_text(
        
        size = 10,
        
        angle = 45,
        
        hjust = 1,
        
        vjust = 1,
        
        color = "black",
        
        margin =
          margin(
            t = 3
          )
      ),
    
    
    axis.text.y =
      element_text(
        
        size = 10.5,
        
        color = "black",
        
        margin =
          margin(
            r = 5
          )
      ),
    
    
    # title
    plot.title =
      element_text(
        
        face = "bold",
        
        size = 15,
        
        hjust = 0,
        
        color = "black",
        
        margin =
          margin(
            b = 10
          )
      ),
    
    
    # caption
    plot.caption =
      element_text(
        
        size = 9,
        
        hjust = 0,
        
        color = "grey30",
        
        margin =
          margin(
            t = 8
          )
      ),
    
    
    # legend
    legend.position =
      "right",
    
    
    legend.title =
      element_text(
        size = 10
      ),
    
    
    legend.text =
      element_text(
        size = 9
      ),
    
    
    legend.key.height =
      unit(
        0.5,
        "cm"
      ),
    
    
    # background
    panel.background =
      element_blank(),
    
    
    plot.background =
      element_rect(
        
        fill = "white",
        
        color = NA
        
      ),
    
    
    plot.margin =
      margin(
        12,
        18,
        12,
        12
      )
  )


# ============================================================
# 16. Display
# ============================================================

p


# ============================================================
# 17. Export significant gene pairs
# ============================================================

sig_pairs <- pair_df %>%
  
  filter(
    p < 0.05
  ) %>%
  
  mutate(
    
    interaction =
      ifelse(
        direction > 0,
        "Co-occurrence",
        "Mutual exclusivity"
      )
  ) %>%
  
  arrange(
    p
  ) %>%
  
  dplyr::select(
    
    gene1,
    
    gene2,
    
    both_mutated,
    
    OR,
    
    p,
    
    interaction
    
  )


sig_pairs


# ============================================================
# 18. Save as PDF
# ============================================================

ggsave(
  
  filename =
    "E:/Lung/EGA/PDF/lung_top15_plus5_cooccurrence.pdf",
  
  plot = p,
  
  width = 10,
  
  height = 6.2
  
)






#####################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)

# ============================================================
# 1. Clinical groups
# ============================================================

clinical_lung <- maf_lung@clinical.data %>%
  as.data.frame() %>%
  dplyr::select(
    Tumor_Sample_Barcode,
    phenotype_group
  ) %>%
  distinct()


# ============================================================
# 2. Merge phenotype data into the previously generated snv96 table
# snv96 should contain:
# Tumor_Sample_Barcode
# substitution
# flank
# ============================================================

snv96_group <- snv96 %>%
  left_join(
    clinical_lung,
    by = "Tumor_Sample_Barcode"
  ) %>%
  filter(
    !is.na(phenotype_group)
  )


# ============================================================
# 3. Set the order
# ============================================================

substitution_order <- c(
  "C>A",
  "C>G",
  "C>T",
  "T>A",
  "T>C",
  "T>G"
)

context_order <- c(
  "A.A", "A.C", "A.G", "A.T",
  "C.A", "C.C", "C.G", "C.T",
  "G.A", "G.C", "G.G", "G.T",
  "T.A", "T.C", "T.G", "T.T"
)

phenotype_order <- c(
  "Lung cancer",
  "Metastasis to lung",
  "Benign tumor"
)


# ============================================================
# 4. Summarize the 96 channels for each phenotype
# ============================================================

sbs96_group <- snv96_group %>%
  count(
    phenotype_group,
    substitution,
    flank,
    name = "Count"
  ) %>%
  complete(
    phenotype_group = phenotype_order,
    substitution = substitution_order,
    flank = context_order,
    fill = list(
      Count = 0
    )
  ) %>%
  group_by(
    phenotype_group
  ) %>%
  mutate(
    Relative_contribution =
      Count / sum(Count)
  ) %>%
  ungroup()


# ============================================================
# 5. Factor order
# ============================================================

sbs96_group$phenotype_group <- factor(
  sbs96_group$phenotype_group,
  levels = phenotype_order
)

sbs96_group$substitution <- factor(
  sbs96_group$substitution,
  levels = substitution_order
)

sbs96_group$flank <- factor(
  sbs96_group$flank,
  levels = context_order
)


# ============================================================
# 6. Plot
# Each row = phenotype
# Each column = substitution class
# ============================================================

p_sbs96_group <- ggplot(
  sbs96_group,
  aes(
    x = flank,
    y = Relative_contribution,
    fill = substitution
  )
) +
  
  geom_col(
    width = 0.88,
    color = "black",
    linewidth = 0.25
  ) +
  
  facet_grid(
    phenotype_group ~ substitution,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_fill_manual(
    values = c(
      "C>A" = "#43B5D8",
      "C>G" = "#222222",
      "C>T" = "#D73027",
      "T>A" = "#C9C9C9",
      "T>C" = "#9AC44D",
      "T>G" = "#E6C9C9"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Relative contribution"
  ) +
  
  theme_classic(
    base_size = 10
  ) +
  
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 5.5,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 10
    ),
    
    strip.background = element_rect(
      fill = "#D0D0D0",
      color = "black",
      linewidth = 0.4
    ),
    
    strip.text = element_text(
      size = 8,
      color = "black"
    ),
    
    panel.spacing.x = unit(
      0.05,
      "cm"
    ),
    
    panel.spacing.y = unit(
      0.10,
      "cm"
    )
  )

p_sbs96_group

ggsave(
  "E:/Lung/EGA/PDF/SBS96_by_phenotype.pdf",
  p_sbs96_group,
  width = 12,
  height = 7
)



#############################################


# If not installed
# install.packages("RIdeogram")

library(RIdeogram)
library(dplyr)
library(tidyr)


mut_pos <- maf_lung@data %>%
  as.data.frame() %>%
  dplyr::select(
    Chromosome,
    Start_Position
  ) %>%
  mutate(
    Chromosome = as.character(Chromosome),
    Chromosome = gsub("^chr", "", Chromosome),
    Start_Position = as.numeric(Start_Position)
  ) %>%
  filter(
    Chromosome %in% c(as.character(1:22), "X"),
    !is.na(Start_Position)
  )

head(mut_pos)


chr_length <- data.frame(
  chr = c(as.character(1:22), "X"),
  length = c(
    248956422,
    242193529,
    198295559,
    190214555,
    181538259,
    170805979,
    159345973,
    145138636,
    138394717,
    133797422,
    135086622,
    133275309,
    114364328,
    107043718,
    101991189,
    90338345,
    83257441,
    80373285,
    58617616,
    64444167,
    46709983,
    50818468,
    156040895
  )
)


window_size <- 5e6


windows <- lapply(
  seq_len(nrow(chr_length)),
  function(i) {
    
    chr_i <- chr_length$chr[i]
    len_i <- chr_length$length[i]
    
    starts <- seq(
      1,
      len_i,
      by = window_size
    )
    
    ends <- pmin(
      starts + window_size - 1,
      len_i
    )
    
    data.frame(
      chr = chr_i,
      start = starts,
      end = ends
    )
  }
) %>%
  bind_rows()


mut_pos <- mut_pos %>%
  mutate(
    start = floor(
      (Start_Position - 1) / window_size
    ) * window_size + 1
  )


density_count <- mut_pos %>%
  dplyr::count(
    Chromosome,
    start,
    name = "mutation_count"
  ) %>%
  dplyr::rename(
    chr = Chromosome
  )

head(density_count)


density_data <- windows %>%
  left_join(
    density_count,
    by = c(
      "chr",
      "start"
    )
  ) %>%
  mutate(
    mutation_count =
      replace_na(
        mutation_count,
        0
      )
  )



density_data <- density_data %>%
  mutate(
    window_Mb =
      (end - start + 1) / 1e6,
    
    density =
      mutation_count /
      window_Mb
  )

summary(density_data$density)


density_data <- density_data %>%
  mutate(
    value = log10(
      density + 1
    )
  )


ideo_density <- density_data %>%
  dplyr::select(
    chr,
    start,
    end,
    value
  )

head(ideo_density)


data(
  human_karyotype,
  package = "RIdeogram"
)

head(human_karyotype)


human_karyotype$Chr <- gsub(
  "^chr",
  "",
  human_karyotype$Chr
)


library(dplyr)
library(RIdeogram)

# ============================================================
# 1. Convert to the format required by RIdeogram
# ============================================================

ideo_input <- ideo_density %>%
  transmute(
    Chr   = as.character(chr),
    Start = as.numeric(start),
    End   = as.numeric(end),
    Value = as.numeric(value)
  )

head(ideo_input)
str(ideo_input)
ideo_input <- ideo_input %>%
  filter(
    !is.na(Chr),
    is.finite(Start),
    is.finite(End),
    is.finite(Value)
  )

ideogram(
  karyotype = human_karyotype,
  overlaid = ideo_input,
  label = NULL,
  output = "E:/Lung/EGA/PDF/lung_mutation_density_ideogram.svg"
)










##################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ============================================================
# 1. Clinical information
# ============================================================

clinical_lung <- maf_lung@clinical.data %>%
  as.data.frame() %>%
  dplyr::select(
    Tumor_Sample_Barcode,
    phenotype_group
  ) %>%
  distinct()


# ============================================================
# 2. Merge mutation and phenotype data
# ============================================================

var_group <- maf_lung@data %>%
  as.data.frame() %>%
  dplyr::select(
    Tumor_Sample_Barcode,
    Variant_Classification
  ) %>%
  left_join(
    clinical_lung,
    by = "Tumor_Sample_Barcode"
  ) %>%
  filter(
    !is.na(phenotype_group)
  )


# ============================================================
# 3. Summarize
# ============================================================

var_stat <- var_group %>%
  count(
    phenotype_group,
    Variant_Classification,
    name = "Count"
  ) %>%
  group_by(
    phenotype_group
  ) %>%
  mutate(
    Proportion = Count / sum(Count)
  ) %>%
  ungroup()


# ============================================================
# 4. Set the order
# ============================================================

var_stat$phenotype_group <- factor(
  var_stat$phenotype_group,
  levels = c(
    "Lung cancer",
    "Metastasis to lung",
    "Benign tumor"
  )
)


# ============================================================
# 5. Plot
# ============================================================

p_var <- ggplot(
  var_stat,
  aes(
    x = phenotype_group,
    y = Proportion,
    fill = Variant_Classification
  )
) +
  
  geom_col(
    width = 0.68,
    color = "white",
    linewidth = 0.35
  ) +
  
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  
  labs(
    x = NULL,
    y = "Proportion of variants",
    fill = "Variant classification"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1,
      color = "black"
    ),
    
    axis.text.y = element_text(
      color = "black"
    ),
    
    legend.position = "right"
  )

p_var

