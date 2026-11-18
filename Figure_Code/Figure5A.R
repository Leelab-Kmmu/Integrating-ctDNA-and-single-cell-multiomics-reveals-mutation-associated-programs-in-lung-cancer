# ============================================================
# E2F3-centred TF-target network (complete standalone script)
#
# Required objects already present in the R environment:
#   final_results : TF significance table
#   subgrn        : TF-target regulatory network table
#
# Main result:
#   E2F3 is fixed at the exact centre. The remaining network uses
#   graphlayouts::layout_with_focus() and a compact Metascape-like style.
#   Network content is not changed by the layout or colours.
# ============================================================


# ------------------------------------------------------------
# 0. User settings
# ------------------------------------------------------------

out_dir <- "E:/Lung/EGA/PDF/gene_list/network_graphopt/code_network"

top_n <- 7L
shared_target_top_k <- 2L

central_tf <- "E2F3"
preferred_central_celltype <- "Tumor"

# "tf"  : label TFs only (recommended for the compact figure)
# "all" : label TFs and every target gene
# "none": label only E2F3
label_mode <- "tf"

# FALSE resembles the compact Metascape network; TRUE shows regulation arrows.
show_arrows <- FALSE

# TRUE colours target genes by their strongest displayed regulator's cell type.
colour_targets_by_regulator <- TRUE

seed <- 1234L

width_mm <- 183
height_mm <- 150
dpi <- 600

# Central E2F3 styling. It should remain visually important without becoming
# a large bullseye or covering its first-neighbour nodes.
central_outer_size <- 7.2
central_inner_size <- 5.5
central_border_width <- 0.55
central_label_size <- 2.65
central_label_nudge_y <- -0.30

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

required_packages <- c(
  "data.table",
  "igraph",
  "ggraph",
  "graphlayouts",
  "ggplot2",
  "ggrepel",
  "scales",
  "svglite",
  "ragg"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them with:\n",
    "install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggraph)
  library(graphlayouts)
  library(ggplot2)
  library(ggrepel)
})


# ------------------------------------------------------------
# 2. Input checks
# ------------------------------------------------------------

if (!exists("final_results", inherits = TRUE)) {
  stop("Object 'final_results' was not found. Load the TF result table first.")
}

if (!exists("subgrn", inherits = TRUE)) {
  stop("Object 'subgrn' was not found. Load the GRN table first.")
}

final_results_dt <- data.table::as.data.table(final_results)
subgrn_dt <- data.table::as.data.table(subgrn)

required_final_cols <- c(
  "TF", "Celltype", "CTS_P", "MGRS_P", "MRS_P", "MRS_null_z"
)

required_grn_cols <- c(
  "TF", "Celltype", "gene", "CellEdgeScore",
  "Target_CTS", "Peak_CTS", "peak_id"
)

missing_final_cols <- setdiff(required_final_cols, names(final_results_dt))
missing_grn_cols <- setdiff(required_grn_cols, names(subgrn_dt))

if (length(missing_final_cols) > 0L) {
  stop(
    "Missing columns in final_results: ",
    paste(missing_final_cols, collapse = ", ")
  )
}

if (length(missing_grn_cols) > 0L) {
  stop(
    "Missing columns in subgrn: ",
    paste(missing_grn_cols, collapse = ", ")
  )
}

if (!label_mode %in% c("tf", "all", "none")) {
  stop("label_mode must be one of: 'tf', 'all', 'none'.")
}


# ------------------------------------------------------------
# 3. Safe summary helpers
# ------------------------------------------------------------

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}


# ------------------------------------------------------------
# 4. Significant TFs and Top-N targets
# ------------------------------------------------------------

sig_tf <- unique(
  final_results_dt[
    !is.na(CTS_P) & !is.na(MGRS_P) & !is.na(MRS_P) &
      CTS_P < 0.05 & MGRS_P < 0.05 & MRS_P < 0.05,
    .(TF, Celltype, MRS_null_z)
  ]
)

if (nrow(sig_tf) == 0L) {
  stop("No TF passed CTS_P < 0.05, MGRS_P < 0.05 and MRS_P < 0.05.")
}

net_raw <- merge(
  subgrn_dt,
  sig_tf,
  by = c("TF", "Celltype"),
  all = FALSE
)

if (nrow(net_raw) == 0L) {
  stop("No GRN rows remained after merging significant TFs with subgrn.")
}

net_gene <- net_raw[
  !is.na(gene) & nzchar(gene),
  .(
    EdgeScore = safe_max(CellEdgeScore),
    MeanEdgeScore = safe_mean(CellEdgeScore),
    Target_CTS = safe_max(Target_CTS),
    Peak_CTS = safe_max(Peak_CTS),
    N_peaks = data.table::uniqueN(peak_id)
  ),
  by = .(TF, Celltype, gene)
]

net_gene <- net_gene[is.finite(EdgeScore)]

data.table::setorder(net_gene, TF, Celltype, -EdgeScore, gene)

net_top <- net_gene[
  , head(.SD, top_n),
  by = .(TF, Celltype)
]

net_top[, TF_node := paste(TF, Celltype, sep = " | ")]


# ------------------------------------------------------------
# 5. Complete edge and node tables
# ------------------------------------------------------------

edge_complete <- net_top[
  , .(
    from = TF_node,
    to = gene,
    EdgeScore,
    MeanEdgeScore,
    Target_CTS,
    Peak_CTS,
    N_peaks,
    TF,
    Celltype
  )
]

tf_nodes <- unique(
  net_top[
    , .(
      name = TF_node,
      label = TF,
      Type = Celltype,
      node_type = "TF"
    )
  ]
)

gene_nodes <- unique(
  net_top[
    , .(
      name = gene,
      label = gene,
      Type = "Target",
      node_type = "Gene"
    )
  ]
)

node_df <- data.table::rbindlist(
  list(tf_nodes, gene_nodes),
  use.names = TRUE,
  fill = TRUE
)

node_df <- unique(node_df, by = "name")

gene_degree <- edge_complete[
  , .(degree = data.table::uniqueN(from)),
  by = .(name = to)
]

node_df <- merge(node_df, gene_degree, by = "name", all.x = TRUE)
node_df[is.na(degree), degree := 0L]


# ------------------------------------------------------------
# 6. Locate the E2F3 TF node
# ------------------------------------------------------------

preferred_name <- paste(
  central_tf,
  preferred_central_celltype,
  sep = " | "
)

central_candidates <- node_df[
  node_type == "TF" & label == central_tf,
  name
]

if (length(central_candidates) == 0L) {
  stop("The central TF was not found: ", central_tf)
}

if (preferred_name %in% central_candidates) {
  e2f3_node <- preferred_name
} else {
  e2f3_node <- central_candidates[1L]
  warning(
    "Preferred central node was absent; using ", e2f3_node, " instead."
  )
}


# ------------------------------------------------------------
# 7. Displayed edges
#
# All single-parent edges are retained. E2F3 edges are retained. For targets
# with multiple parents, only the strongest K edges are displayed. The complete
# edge table is always exported unchanged.
# ------------------------------------------------------------

parent_count <- edge_complete[
  , .(N_parent = data.table::uniqueN(from)),
  by = to
]

edge_ranked <- merge(edge_complete, parent_count, by = "to", all.x = TRUE)
data.table::setorder(edge_ranked, to, -EdgeScore, from)

edge_single <- edge_ranked[N_parent == 1L]
edge_central <- edge_ranked[from == e2f3_node]
edge_shared <- edge_ranked[N_parent >= 2L]

edge_shared_top <- edge_shared[
  , head(.SD, shared_target_top_k),
  by = to
]

edge_plot <- unique(
  data.table::rbindlist(
    list(edge_single, edge_central, edge_shared_top),
    use.names = TRUE,
    fill = TRUE
  ),
  by = c("from", "to")
)

used_names <- unique(c(edge_plot$from, edge_plot$to))
nodes_plot <- node_df[name %in% used_names]

graph_plot <- igraph::graph_from_data_frame(
  d = as.data.frame(edge_plot),
  vertices = as.data.frame(nodes_plot),
  directed = TRUE
)


# ------------------------------------------------------------
# 8. Focus layout: E2F3 at the exact centre
# ------------------------------------------------------------

graph_undirected <- igraph::as_undirected(
  graph_plot,
  mode = "collapse",
  edge.attr.comb = "first"
)

focus_index <- match(e2f3_node, igraph::V(graph_undirected)$name)

if (is.na(focus_index)) {
  stop("E2F3 is absent from the displayed graph: ", e2f3_node)
}

# layout_with_focus requires one connected graph. If the displayed network has
# multiple components, add layout-only helper links from E2F3 to one vertex in
# each remote component. These links are never added to graph_plot or exported.
graph_for_layout <- graph_undirected
component_info <- igraph::components(graph_for_layout)

if (component_info$no > 1L) {
  focus_component <- component_info$membership[focus_index]
  remote_components <- setdiff(
    seq_len(component_info$no),
    focus_component
  )

  representative_indices <- vapply(
    remote_components,
    function(component_id) {
      which(component_info$membership == component_id)[1L]
    },
    integer(1)
  )

  helper_edge_vector <- as.vector(
    rbind(
      rep(focus_index, length(representative_indices)),
      representative_indices
    )
  )

  graph_for_layout <- igraph::add_edges(
    graph_for_layout,
    helper_edge_vector,
    attr = list(layout_helper = TRUE)
  )
}

focus_layout_result <- tryCatch(
  graphlayouts::layout_with_focus(
    graph_for_layout,
    v = focus_index
  ),
  error = function(e) {
    warning(
      "layout_with_focus failed; using stress layout. Reason: ",
      conditionMessage(e)
    )
    graphlayouts::layout_with_stress(graph_for_layout)
  }
)

# graphlayouts >= 1.2 returns list(xy = ..., distance = ...) for the focus
# layout, while stress layout returns a matrix.
if (is.list(focus_layout_result) && !is.null(focus_layout_result$xy)) {
  layout_xy <- focus_layout_result$xy
} else {
  layout_xy <- focus_layout_result
}

layout_xy <- as.matrix(layout_xy)

# Put E2F3 at (0, 0), then apply equal scaling to both axes.
layout_xy[, 1] <- layout_xy[, 1] - layout_xy[focus_index, 1]
layout_xy[, 2] <- layout_xy[, 2] - layout_xy[focus_index, 2]

layout_scale <- max(abs(layout_xy), na.rm = TRUE)
if (is.finite(layout_scale) && layout_scale > 0) {
  layout_xy <- layout_xy / layout_scale * 6
}

coordinates <- data.table::data.table(
  name = igraph::V(graph_plot)$name,
  x = layout_xy[, 1],
  y = layout_xy[, 2]
)

# Absolute final guarantee.
coordinates[name == e2f3_node, `:=`(x = 0, y = 0)]

lay <- ggraph::create_layout(
  graph_plot,
  layout = "manual",
  x = coordinates$x,
  y = coordinates$y
)

node_meta <- data.table::as.data.table(nodes_plot)

lay$label <- node_meta$label[match(lay$name, node_meta$name)]
lay$node_type <- node_meta$node_type[match(lay$name, node_meta$name)]
lay$Type <- node_meta$Type[match(lay$name, node_meta$name)]
lay$degree <- node_meta$degree[match(lay$name, node_meta$name)]
lay$is_central <- lay$name == e2f3_node


# ------------------------------------------------------------
# 9. Target colours and communities
# ------------------------------------------------------------

strongest_regulator <- data.table::copy(edge_plot)
data.table::setorder(strongest_regulator, to, -EdgeScore, from)

strongest_regulator <- strongest_regulator[
  , .SD[1L],
  by = to
][
  , .(name = to, regulator_type = Celltype)
]

target_type_map <- stats::setNames(
  strongest_regulator$regulator_type,
  strongest_regulator$name
)

lay$PlotType <- lay$Type

if (isTRUE(colour_targets_by_regulator)) {
  gene_rows <- lay$node_type == "Gene"
  lay$PlotType[gene_rows] <- unname(
    target_type_map[lay$name[gene_rows]]
  )
}

lay$PlotType[is.na(lay$PlotType) | !nzchar(lay$PlotType)] <- "Target"

if (igraph::ecount(graph_undirected) > 0L && igraph::vcount(graph_undirected) > 1L) {
  community_fit <- igraph::cluster_louvain(graph_undirected)
  membership_vector <- igraph::membership(community_fit)
  lay$community <- as.integer(membership_vector[lay$name])
} else {
  lay$community <- 1L
}

lay$plot_size <- ifelse(
  lay$node_type == "TF",
  4.2,
  scales::rescale(
    pmax(lay$degree, 1),
    to = c(1.8, 3.2),
    from = range(pmax(lay$degree, 1), na.rm = TRUE)
  )
)

lay$plot_size[lay$is_central] <- central_inner_size


# ------------------------------------------------------------
# 10. Colour palette
# ------------------------------------------------------------

type_cols <- c(
  "AT1" = "#E64B35",
  "AT2" = "#D99A00",
  "Endothelial" = "#3C9FB7",
  "Fibroblast" = "#405C91",
  "Lymphoid" = "#6968B3",
  "Myeloid" = "#C94A87",
  "Tumor" = "#70A650",
  "Target" = "#C9CED6"
)

present_types <- sort(unique(lay$PlotType))
unknown_types <- setdiff(present_types, names(type_cols))

if (length(unknown_types) > 0L) {
  extra_cols <- scales::hue_pal(l = 60, c = 70)(length(unknown_types))
  names(extra_cols) <- unknown_types
  type_cols <- c(type_cols, extra_cols)
}

legend_breaks <- intersect(
  c(
    "AT1", "AT2", "Endothelial", "Fibroblast",
    "Lymphoid", "Myeloid", "Tumor", "Target",
    unknown_types
  ),
  present_types
)


# ------------------------------------------------------------
# 11. Draw the network
# ------------------------------------------------------------

arrow_spec <- if (isTRUE(show_arrows)) {
  grid::arrow(length = grid::unit(0.8, "mm"), type = "closed")
} else {
  NULL
}

plot_limit <- max(abs(c(lay$x, lay$y)), na.rm = TRUE) * 1.12

p_network <- ggraph::ggraph(lay) +
  ggraph::geom_edge_link(
    colour = "#7566A8",
    linewidth = 0.28,
    alpha = 0.38,
    lineend = "round",
    arrow = arrow_spec,
    end_cap = ggraph::circle(1.2, "mm"),
    show.legend = FALSE
  ) +
  ggraph::geom_node_point(
    data = function(x) x[x$node_type == "Gene", ],
    ggplot2::aes(fill = PlotType, size = plot_size),
    shape = 21,
    colour = "white",
    stroke = 0.28,
    alpha = 0.90,
    show.legend = FALSE
  ) +
  ggraph::geom_node_point(
    data = function(x) x[x$node_type == "TF" & !x$is_central, ],
    ggplot2::aes(fill = PlotType, size = plot_size),
    shape = 21,
    colour = "white",
    stroke = 0.65,
    alpha = 1,
    show.legend = TRUE
  ) +
  ggraph::geom_node_point(
    data = function(x) x[x$is_central, ],
    shape = 21,
    size = central_outer_size,
    fill = "white",
    colour = "#353535",
    stroke = central_border_width,
    show.legend = FALSE
  ) +
  ggraph::geom_node_point(
    data = function(x) x[x$is_central, ],
    ggplot2::aes(fill = PlotType),
    shape = 21,
    size = central_inner_size,
    colour = "#353535",
    stroke = central_border_width,
    show.legend = FALSE
  ) +
  # White text halo keeps the label readable over nearby links.
  ggraph::geom_node_text(
    data = function(x) x[x$is_central, ],
    ggplot2::aes(label = label),
    colour = "white",
    family = "Arial",
    fontface = "bold",
    size = central_label_size + 0.55,
    nudge_y = central_label_nudge_y,
    show.legend = FALSE
  ) +
  ggraph::geom_node_text(
    data = function(x) x[x$is_central, ],
    ggplot2::aes(label = label),
    colour = "#C93C32",
    family = "Arial",
    fontface = "bold",
    size = central_label_size,
    nudge_y = central_label_nudge_y,
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = type_cols,
    breaks = legend_breaks,
    name = "Cell type",
    drop = FALSE
  ) +
  ggplot2::scale_size_identity() +
  ggplot2::coord_fixed(
    xlim = c(-plot_limit, plot_limit),
    ylim = c(-plot_limit, plot_limit),
    clip = "off"
  ) +
  ggplot2::theme_void(base_family = "Arial", base_size = 7) +
  ggplot2::theme(
    legend.position = "right",
    legend.title = ggplot2::element_text(
      family = "Arial", face = "bold", size = 7
    ),
    legend.text = ggplot2::element_text(
      family = "Arial", size = 6.2
    ),
    legend.key.height = grid::unit(0.40, "cm"),
    legend.spacing.y = grid::unit(0.02, "cm"),
    plot.margin = ggplot2::margin(7, 12, 7, 7),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA)
  )


# ------------------------------------------------------------
# 12. Optional labels
# ------------------------------------------------------------

if (label_mode %in% c("tf", "all")) {
  p_network <- p_network +
    ggrepel::geom_text_repel(
      data = as.data.frame(lay[lay$node_type == "TF" & !lay$is_central, ]),
      ggplot2::aes(x = x, y = y, label = label),
      family = "Arial",
      fontface = "bold",
      size = 2.5,
      colour = "#252525",
      box.padding = 0.20,
      point.padding = 0.15,
      force = 0.65,
      force_pull = 0.08,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.colour = "#AFAFAF",
      segment.linewidth = 0.16,
      seed = seed,
      inherit.aes = FALSE
    )
}

if (label_mode == "all") {
  p_network <- p_network +
    ggrepel::geom_text_repel(
      data = as.data.frame(lay[lay$node_type == "Gene", ]),
      ggplot2::aes(x = x, y = y, label = label),
      family = "Arial",
      size = 1.65,
      colour = "#555555",
      box.padding = 0.08,
      point.padding = 0.05,
      force = 0.30,
      force_pull = 0.03,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.colour = "#C4C4C4",
      segment.linewidth = 0.10,
      seed = seed,
      inherit.aes = FALSE
    )
}

print(p_network)


# ------------------------------------------------------------
# 13. Export editable, vector and raster files
# ------------------------------------------------------------

export_base <- file.path(
  out_dir,
  paste0("E2F3_Top", top_n, "_focus_network")
)

width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

svglite::svglite(
  paste0(export_base, ".svg"),
  width = width_in,
  height = height_in,
  bg = "white"
)
print(p_network)
grDevices::dev.off()

grDevices::cairo_pdf(
  paste0(export_base, ".pdf"),
  width = width_in,
  height = height_in,
  family = "Arial",
  bg = "white"
)
print(p_network)
grDevices::dev.off()

ragg::agg_png(
  paste0(export_base, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = dpi,
  background = "white"
)
print(p_network)
grDevices::dev.off()

ragg::agg_tiff(
  paste0(export_base, ".tiff"),
  width = width_in,
  height = height_in,
  units = "in",
  res = dpi,
  compression = "lzw",
  background = "white"
)
print(p_network)
grDevices::dev.off()


# ------------------------------------------------------------
# 14. Export source tables and reproducibility objects
# ------------------------------------------------------------

coordinate_export <- data.table::as.data.table(lay)[
  , .(
    name,
    label,
    node_type,
    Type,
    PlotType,
    degree,
    community,
    x,
    y,
    is_central
  )
]

data.table::fwrite(
  edge_complete,
  paste0(export_base, "_edges_complete.tsv"),
  sep = "\t",
  na = ""
)

data.table::fwrite(
  edge_plot,
  paste0(export_base, "_edges_displayed.tsv"),
  sep = "\t",
  na = ""
)

data.table::fwrite(
  coordinate_export,
  paste0(export_base, "_nodes_coordinates.tsv"),
  sep = "\t",
  na = ""
)

saveRDS(
  list(
    plot = p_network,
    graph = graph_plot,
    layout = coordinate_export,
    edges_complete = edge_complete,
    edges_displayed = edge_plot,
    nodes = nodes_plot,
    settings = list(
      top_n = top_n,
      shared_target_top_k = shared_target_top_k,
      central_node = e2f3_node,
      label_mode = label_mode,
      show_arrows = show_arrows,
      seed = seed
    )
  ),
  paste0(export_base, ".rds")
)


# ------------------------------------------------------------
# 15. Summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("E2F3 focus network completed\n")
cat("Central node:", e2f3_node, "\n")
cat("Central coordinates: x=0, y=0\n")
cat("TF nodes:", sum(lay$node_type == "TF"), "\n")
cat("Target genes:", sum(lay$node_type == "Gene"), "\n")
cat("Complete edges:", nrow(edge_complete), "\n")
cat("Displayed edges:", nrow(edge_plot), "\n")
cat("Detected communities:", data.table::uniqueN(lay$community), "\n")
cat("Label mode:", label_mode, "\n")
cat("Output base:", export_base, "\n")
cat("============================================\n")









ggsave(
  filename = file.path(
    out_dir,
    "E2F3_Top7_focus_network_large.pdf"
  ),
  plot = p_network,
  width = 14,
  height = 11,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
