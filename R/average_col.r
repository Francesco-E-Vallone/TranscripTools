#' Average Columns by Group (e.g., Replicates or Conditions)
#'
#' Computes the mean expression for each gene across replicates that belong
#' to the same group, as defined in a metadata data frame.
#'
#' @param tmm Numeric matrix or data frame of expression values
#'   (rows = genes, columns = samples).
#' @param meta Data frame with sample metadata.
#'   Must contain a column named \code{group} and row names matching the
#'   column names of \code{tmm}.
#'
#' @return A data frame with the same row names as \code{tmm} and one column
#'   per group, containing the average expression values for that group.
#'
#' @details
#' This function is useful when visualizing summarized or “cumulated”
#' profiles, such as average expression per condition for heatmaps or PCA.
#' Missing values are ignored (\code{na.rm = TRUE}).
#'
#' @examples
#' \dontrun{
#' # Average normalized counts by condition
#' tmm_avg <- average_cols(tmm, meta)
#'
#' # Use averaged matrix for a heatmap
#' hmap(tmm_avg, title = "Group-Averaged Expression")
#' }
#'
#' @export
average_cols <- function(tmm, meta) {
  if (is.data.frame(tmm)) tmm <- as.matrix(tmm)
  if (!"group" %in% colnames(meta))
    stop("`meta` must contain a column named 'group'.")
  if (is.null(rownames(meta)) || any(!colnames(tmm) %in% rownames(meta)))
    stop("Row names of `meta` must match column names of `tmm`.")
  
  averaged_df <- data.frame(row.names = rownames(tmm))
  
  for (grp in levels(factor(meta$group))) {
    samples <- rownames(meta[meta$group == grp, , drop = FALSE])
    avg_values <- rowMeans(tmm[, samples, drop = FALSE], na.rm = TRUE)
    averaged_df[[grp]] <- avg_values
  }
  
  averaged_df
}