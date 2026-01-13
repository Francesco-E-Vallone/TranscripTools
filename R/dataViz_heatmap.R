#' Create a Customized Heatmap for Count Data
#'
#' This function generates a heatmap from a normalized count matrix (e.g., TMM or TPM data),
#' with the option to restrict the visualization to a set of differentially expressed genes (DEGs).
#' It supports a top annotation built from sample metadata and lets you control whether row/column
#' names are shown, plus their font size and font face (via ComplexHeatmap's \code{*_names_gp}).
#'
#' If \code{annotation_colors} is provided but does not cover all annotation levels,
#' missing levels are filled with \code{"grey80"} and a warning is emitted.
#'
#' @param tmm Numeric matrix or data.frame (rows = genes, cols = samples).
#' @param title Character. Heatmap title.
#' @param deg Optional character vector of genes to include (must match \code{rownames(tmm)}).
#' @param heatmap_colors Color vector for the heatmap gradient.
#' @param metadata Optional data.frame with sample metadata (rows = samples).
#' @param metadata_col Optional character. Column in \code{metadata} to use for top annotation.
#' @param annotation_colors Optional named vector mapping annotation levels to colors.
#' @param show_row_names,show_column_names Logical. Show gene/sample names. Default FALSE.
#' @param row_names_fontsize Numeric or NULL. Font size (points) for row names. NULL = default.
#' @param row_names_fontface Character/integer or NULL. Font face for row names
#'   (e.g., \code{"plain"}, \code{"bold"}, \code{"italic"}, \code{"bold.italic"}). NULL = default.
#' @param column_names_fontsize Numeric or NULL. Font size (points) for column names. NULL = default.
#' @param column_names_fontface Character/integer or NULL. Font face for column names
#'   (e.g., \code{"plain"}, \code{"bold"}, \code{"italic"}, \code{"bold.italic"}). NULL = default.
#'
#' @return A \code{ComplexHeatmap::Heatmap} object (drawn and returned).
#' @export
hmap <- function(tmm,
                 title,
                 deg = NULL,
                 heatmap_colors = c("#5a78b3", "#93b2d3", "#e7f4f2",
                                    "#faf0ab", "#db7d54", "#b9332c"),
                 metadata = NULL,
                 metadata_col = NULL,
                 annotation_colors = NULL,
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 row_names_fontsize = NULL,
                 row_names_fontface = NULL,
                 column_names_fontsize = NULL,
                 column_names_fontface = NULL) {
  
  # ---- dependency check ----
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop(
      "Package 'ComplexHeatmap' is required for hmap(). ",
      "Install it with: BiocManager::install('ComplexHeatmap')"
    )
  }
  
  # ---- input checks / coercion ----
  if (is.data.frame(tmm)) tmm <- as.matrix(tmm)
  if (!is.matrix(tmm)) stop("`tmm` must be a matrix or data.frame coercible to a matrix.")
  
  # subset to DEGs if provided
  if (!is.null(deg)) {
    tmm <- tmm[rownames(tmm) %in% deg, , drop = FALSE]
  }
  
  if (nrow(tmm) == 0) stop("`tmm` has 0 rows after filtering (check `deg`).")
  if (ncol(tmm) == 0) stop("`tmm` has 0 columns (no samples).")
  
  # ---- z-score per gene (row) ----
  zfun <- function(x) {
    sx <- stats::sd(x, na.rm = TRUE)
    if (is.na(sx) || sx == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / sx
  }
  tmm_z <- t(apply(tmm, 1, zfun))
  
  # ---- top annotation (optional) ----
  topAnn <- NULL
  if (!is.null(metadata) && !is.null(metadata_col)) {
    
    if (!is.data.frame(metadata)) stop("`metadata` must be a data.frame.")
    if (is.null(rownames(metadata))) stop("`metadata` must have rownames = sample IDs.")
    if (!(metadata_col %in% colnames(metadata))) {
      stop("`metadata_col` not found in metadata: ", metadata_col)
    }
    
    # align metadata to heatmap columns
    missing_meta <- setdiff(colnames(tmm_z), rownames(metadata))
    if (length(missing_meta) > 0) {
      stop("`metadata` is missing these samples: ", paste(missing_meta, collapse = ", "))
    }
    metadata <- metadata[colnames(tmm_z), , drop = FALSE]
    
    ann_df <- data.frame(group = metadata[[metadata_col]])
    rownames(ann_df) <- colnames(tmm_z)
    
    if (!is.null(annotation_colors)) {
      
      if (is.null(names(annotation_colors)) || anyNA(names(annotation_colors)) || any(names(annotation_colors) == "")) {
        stop("`annotation_colors` must be a *named* vector: names = levels, values = colors.")
      }
      
      levs <- unique(as.character(ann_df$group))
      missing_cols <- setdiff(levs, names(annotation_colors))
      
      if (length(missing_cols) > 0) {
        warning(
          "`annotation_colors` missing colors for: ",
          paste(missing_cols, collapse = ", "),
          ". Using 'grey80' for missing levels."
        )
        annotation_colors <- c(
          annotation_colors,
          stats::setNames(rep("grey80", length(missing_cols)), missing_cols)
        )
      }
      
      # ensure order matches levels present
      annotation_colors <- annotation_colors[levs]
      
      topAnn <- ComplexHeatmap::HeatmapAnnotation(
        df  = ann_df,
        col = list(group = annotation_colors)
      )
      
    } else {
      topAnn <- ComplexHeatmap::HeatmapAnnotation(df = ann_df)
    }
  }
  
  # ---- ALWAYS pass valid gpar() objects (ComplexHeatmap doesn't like NULL) ----
  row_gp <- grid::gpar()
  if (!is.null(row_names_fontsize)) row_gp$fontsize <- row_names_fontsize
  if (!is.null(row_names_fontface)) row_gp$fontface <- row_names_fontface
  
  col_gp <- grid::gpar()
  if (!is.null(column_names_fontsize)) col_gp$fontsize <- column_names_fontsize
  if (!is.null(column_names_fontface)) col_gp$fontface <- column_names_fontface
  
  # ---- heatmap ----
  ht <- ComplexHeatmap::Heatmap(
    tmm_z,
    name = "Z.score",
    col = heatmap_colors,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    row_names_gp = row_gp,
    column_names_gp = col_gp,
    show_row_dend = FALSE,
    border = TRUE,
    top_annotation = topAnn,
    column_title = title
  )
  
  ComplexHeatmap::draw(ht)
  return(ht)
}