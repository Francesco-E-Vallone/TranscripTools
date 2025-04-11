#' Create a Customized Heatmap for Count Data
#'
#' This function generates a heatmap from a normalized count matrix (e.g., TMM or TPM data), with the option to restrict
#' the visualization to a set of differentially expressed genes (DEGs). In addition, it allows full customization of the heatmap
#' colors, the top annotation, and the display of row (gene names) and column names.
#'
#' @param tmm A numeric matrix or data frame with row names representing genes.
#' @param title A character string to be used as the heatmap's column title.
#' @param deg An optional character vector of gene names. If provided, only the rows (genes) present in \code{deg} will be used.
#'   If \code{NULL} (default), the heatmap will be created using all rows in \code{tmm}.
#' @param heatmap_colors A character vector specifying the color palette for the heatmap. Default is a blue-to-red gradient.
#' @param top_annotation_df An optional data frame for top annotation (for example, sample grouping information). If provided,
#'   it will be used to annotate columns; if \code{NULL} (default) and a global \code{metadata$samples} exists, that will be used.
#' @param top_annotation_colors An optional named list of colors for the annotation. If \code{NULL} (default), the annotation
#'   will use the default (random) colors.
#' @param annotation_width A \code{unit} object specifying the annotation width. Default is \code{unit(c(1, 4), 'cm')}.
#' @param show_row_names Logical. Whether to display row (gene) names on the heatmap. Default is \code{FALSE}.
#' @param show_column_names Logical. Whether to display column names on the heatmap. Default is \code{FALSE}.
#'
#' @return A \code{Heatmap} object (from the ComplexHeatmap package) that is drawn and returned.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Optionally subsets the count matrix \code{tmm} to include only rows with gene names in \code{deg} (if provided).
#'   \item Computes z-scores for each gene across samples.
#'   \item Creates a top annotation if \code{top_annotation_df} is provided. If \code{top_annotation_colors} is also provided,
#'         it is used; otherwise, default colors are applied.
#'   \item Draws a heatmap using the \code{Heatmap} function from the ComplexHeatmap package with user-specified options for
#'         displaying row and column names.
#' }
#'
#' @examples
#' \dontrun{
#' # Example: Create a heatmap using all genes
#' ht_all <- hmap(tmm_avg, title = "Heatmap of All Genes")
#'
#' # Example: Create a heatmap using only a subset of DEGs, with custom annotation colors and showing column names:
#' deg_genes <- rownames(res_PB_vs_all)
#' custom_colors <- list(samples = c("NM_IgM" = "darkred", "NM_NT" = "darkorange", 
#'                                   "NWT_IgM" = "darkblue", "NWT_NT" = "lightblue"))
#' ht_deg <- hmap(tmm_avg, title = "PB vs ALL - DEGs", deg = deg_genes,
#'                heatmap_colors = c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"),
#'                top_annotation_df = data.frame(samples = metadata$samples),
#'                top_annotation_colors = custom_colors,
#'                show_column_names = TRUE)
#' }
#'
#' @import ComplexHeatmap
#' @import grid
#' @export
hmap <- function(tmm, title, deg = NULL,
                 heatmap_colors = c("#5a78b3", "#6381b8", "#93b2d3", "#c6dbea", "#e7f4f2",
                                    "#faf0ab", "#f6e298", "#db7d54", "#c14536", "#b9332c"),
                 top_annotation_df = NULL,
                 top_annotation_colors = NULL,
                 annotation_width = unit(c(1, 4), 'cm'),
                 show_row_names = FALSE,
                 show_column_names = FALSE) {
  #if a deg vector is provided, subset the count matrix
  if (!is.null(deg)) {
    tmm <- tmm[rownames(tmm) %in% deg, ]
  }
  
  #compute z-scores for each gene
  cal_z_score <- function(x) { (x - mean(x)) / sd(x) }
  tmm_n <- apply(tmm, c(1, 2), as.numeric)
  tmm_norm <- t(apply(tmm_n, 1, cal_z_score))
  tmm_norm <- na.omit(tmm_norm)
  
  # Use global metadata if top annotation not provided and metadata exists
  if (is.null(top_annotation_df) && exists("metadata") && !is.null(metadata$samples)) {
    top_annotation_df <- data.frame(samples = metadata$samples)
  }
  
  #create top annotation if top_annotation_df is provided
  topAnn <- NULL
  if (!is.null(top_annotation_df)) {
    if (!is.null(top_annotation_colors)) {
      topAnn <- HeatmapAnnotation(df = top_annotation_df,
                                  col = top_annotation_colors,
                                  annotation_width = annotation_width)
    } else {
      topAnn <- HeatmapAnnotation(df = top_annotation_df,
                                  annotation_width = annotation_width)
    }
  }
  
  #draw heatmap
  ht <- Heatmap(tmm_norm,
                show_row_names = show_row_names,
                show_column_names = show_column_names,
                show_row_dend = FALSE,
                border = TRUE,
                clustering_distance_rows = "euclidean",
                clustering_distance_columns = "euclidean",
                row_names_gp = gpar(fontsize = 6, fontface = 3),
                top_annotation = topAnn,
                col = heatmap_colors,
                name = "Z.score",
                column_title = title)
  ht <- draw(ht)
  return(ht)
}
