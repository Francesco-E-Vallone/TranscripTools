#' Average Specified Column Groups in a Count Matrix
#'
#' This function computes the row-wise average of specified groups of columns from a count matrix.
#' It returns a new data frame with averaged values and retains the original row names (typically gene names).
#'
#' @param tmm A numeric matrix or data frame with row names representing genes and columns representing samples.
#' @param pairs A list where each element is a character vector of column names in \code{tmm} that should be averaged together.
#'
#' @return A data frame with the same row names as \code{tmm} and one column per element in \code{pairs}. Each column contains
#'   the average of the corresponding columns from \code{tmm}.
#'
#' @details
#' The function loops over each group specified in \code{pairs}, computes the row-wise mean (ignoring \code{NA} values),
#' and assigns a new column name by concatenating the first two elements of the pair using the format \code{"col1_n_col2"}.
#' This is useful when you want to simplify and summarize data from replicated samples prior to visualization, such as in a heatmap.
#'
#' @examples
#' \dontrun{
#' # Assume tmm is your count matrix and you have defined column groups to average:
#' pairs <- list(c("counts_RS1050_1_KIDNEY", "counts_RS1316_3_SPLEEN"),
#'               c("counts_RS1050_1_PB", "counts_RS1316_2_PB"),
#'               c("counts_RS1316_1_BM", "counts_RS9737_2_BM"))
#' averaged_tmm <- average_cols(tmm, pairs)
#' }
#'
#' @export
average_cols <- function(tmm, pairs) {
  #create an empty data frame with the same row names as tmm
  averaged_df <- data.frame(row.names = rownames(tmm))

  #loop over each group of columns defined in 'pairs'
  for (pair in pairs) {
    #calculate the row-wise average for the specified columns
    avg_values <- rowMeans(tmm[, pair, drop = FALSE], na.rm = TRUE)

    #create a new column name using the first two column names
    new_col <- paste(pair[1], pair[2], sep = "_n_")

    #append the averaged values to the data frame
    averaged_df[[new_col]] <- avg_values
  }

  return(averaged_df)
}


#' Create a Customized Heatmap for Count Data
#'
#' This function generates a heatmap from a count matrix (e.g., TMM or TPM data), with the option to restrict to a set
#' of differentially expressed genes (DEGs). In addition, it allows full customization of the heatmap colors and the top annotation.
#'
#' @param tmm A numeric matrix or data frame with row names representing genes.
#' @param title A character string to be used as the heatmap's column title.
#' @param deg An optional character vector of gene names. If provided, only the rows (genes) present in \code{deg} will be used.
#'   If \code{NULL} (default), the heatmap will be created using all rows in \code{tmm}.
#' @param heatmap_colors A character vector specifying the color palette for the heatmap. Default is a blue-to-red gradient.
#' @param top_annotation_df An optional data frame for top annotation. For example, this could contain sample groupings.
#'   If \code{NULL} (default), no top annotation is added unless a global \code{metadata$samples} exists.
#' @param top_annotation_colors A named list of colors for the annotation. Default is \code{NULL}, meaning no annotation colors are applied.
#' @param annotation_width A \code{unit} object specifying the annotation width. Default is \code{unit(c(1, 4), 'cm')}.
#'
#' @return A \code{Heatmap} object (from the ComplexHeatmap package) that is drawn and returned.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Optionally subsets the count matrix \code{tmm} to include only rows in \code{deg} (if provided).
#'   \item Computes z-scores for each gene across samples.
#'   \item Optionally creates a top annotation if both \code{top_annotation_df} and \code{top_annotation_colors} are provided.
#'   \item Draws a heatmap using the \code{Heatmap} function from the ComplexHeatmap package.
#' }
#' This flexibility allows users to visualize either the entire dataset or just a subset of DEGs, with full control over aesthetics.
#'
#' @examples
#' \dontrun{
#' # Example: Create a heatmap using all genes
#' ht_all <- hmap(tmm_avg, title = "Heatmap of All Genes")
#'
#' # Example: Create a heatmap using only a subset of DEGs
#' deg_genes <- rownames(res_PB_vs_all)
#' ht_deg <- hmap(tmm_avg, title = "PB vs ALL - DEGs", deg = deg_genes,
#'                heatmap_colors = c("#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"),
#'                top_annotation_df = data.frame(Group = metadata$samples),
#'                top_annotation_colors = list(Group = c("PB" = "blue", "BM" = "red", "solid_tissues" = "green")))
#' }
#'
#' @import ComplexHeatmap
#' @import grid
#' @export
hmap <- function(tmm, title, deg = NULL,
                 heatmap_colors = c("#5a78b3","#6381b8","#93b2d3", "#c6dbea", "#e7f4f2",
                                    "#faf0ab","#f6e298","#db7d54", "#c14536", "#b9332c"),
                 top_annotation_df = NULL,
                 top_annotation_colors = NULL,
                 annotation_width = unit(c(1, 4), 'cm')) {
  #subset tmm if a DEG list is provided
  if (!is.null(deg)) {
    tmm <- tmm[rownames(tmm) %in% deg, ]
  }

  #compute z-scores for each gene
  cal_z_score <- function(x) { (x - mean(x)) / sd(x) }
  tmm_n <- apply(tmm, c(1, 2), as.numeric)
  tmm_norm <- t(apply(tmm_n, 1, cal_z_score))
  tmm_norm <- na.omit(tmm_norm)

  #if no top annotation is provided, try to use a global metadata if available
  if (is.null(top_annotation_df) && exists("metadata") && !is.null(metadata$samples)) {
    top_annotation_df <- data.frame(samples = metadata$samples)
  }

  #create top annotation if both annotation data and colors are provided
  topAnn <- NULL
  if (!is.null(top_annotation_df) && !is.null(top_annotation_colors)) {
    topAnn <- HeatmapAnnotation(df = top_annotation_df,
                                col = top_annotation_colors,
                                annotation_width = annotation_width)
  }

  #draw the heatmap using ComplexHeatmap
  ht <- Heatmap(tmm_norm,
                show_row_names = FALSE,
                show_row_dend = FALSE,
                show_column_names = FALSE,
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
