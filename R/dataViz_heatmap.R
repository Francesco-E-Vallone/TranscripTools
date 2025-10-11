#' Create a Customized Heatmap for Count Data
#'
#' This function generates a heatmap from a normalized count matrix (e.g., TMM or TPM data), with the option to restrict
#' the visualization to a set of differentially expressed genes (DEGs). In addition, it allows full customization of the heatmap
#' colors, the top annotation, and the display of row (gene names) and column names.
#'
#' @param tmm Numeric matrix or data.frame (rows = genes, cols = samples).
#' @param title Character. Heatmap title.
#' @param deg Optional character vector of genes to include.
#' @param heatmap_colors Color vector for the heatmap gradient.
#' @param metadata Data frame with sample metadata (rows = samples).
#' @param metadata_col Character. Column in \code{metadata} to use for top annotation.
#' @param annotation_colors Optional named vector for annotation levels.
#' @param show_row_names,show_column_names Logical. Show gene/sample names. Default FALSE.
#'
#' @return A \code{ComplexHeatmap::Heatmap} object (drawn and returned invisibly).
#' @export
hmap <- function(tmm, title, deg = NULL,
                 heatmap_colors = c("#5a78b3", "#93b2d3", "#e7f4f2",
                                    "#faf0ab", "#db7d54", "#b9332c"),
                 metadata = NULL,
                 metadata_col = NULL,
                 annotation_colors = NULL,
                 show_row_names = FALSE,
                 show_column_names = FALSE) {
  
  if (is.data.frame(tmm)) tmm <- as.matrix(tmm)
  if (!is.null(deg)) tmm <- tmm[rownames(tmm) %in% deg, , drop = FALSE]
  
  #z-scores per row
  zfun <- function(x) if (sd(x) == 0) rep(0, length(x)) else (x - mean(x)) / sd(x)
  tmm_z <- t(apply(tmm, 1, zfun))
  
  #top annotation (if metadata + column provided)
  topAnn <- NULL
  if (!is.null(metadata) && !is.null(metadata_col) && metadata_col %in% colnames(metadata)) {
    ann_df <- data.frame(group = metadata[[metadata_col]])
    if (!is.null(annotation_colors))
      topAnn <- ComplexHeatmap::HeatmapAnnotation(group = ann_df$group,
                                                  col = list(group = annotation_colors))
    else
      topAnn <- ComplexHeatmap::HeatmapAnnotation(group = ann_df$group)
  }
  
  ht <- ComplexHeatmap::Heatmap(
    tmm_z,
    name = "Z.score",
    col = heatmap_colors,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    show_row_dend = F,
    border = T,
    top_annotation = topAnn,
    column_title = title
  )
  
  ComplexHeatmap::draw(ht)
  return(ht)
}