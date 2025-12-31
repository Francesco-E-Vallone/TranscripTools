#' Perform and Plot PCA with Customizable Aesthetics
#'
#' This function performs Principal Component Analysis (PCA) on a count matrix and generates a biplot
#' using \code{PCAtools}. It streamlines the workflow by computing the PCA and plotting in one step.
#'
#' \strong{Important note on shapes:} In \code{PCAtools::biplot()}, the argument \code{shape} specifies the
#' metadata column used to assign point shapes, while \code{shapekey} is an optional named vector mapping
#' factor levels to specific shape codes. For backward compatibility with older usage patterns in this package,
#' \code{pcafun()} uses \code{shapekey} to mean "the metadata column used for shapes", and exposes the optional
#' mapping via \code{shape_map}.
#'
#' @param counts A numeric matrix or data frame containing expression/count data (e.g., TMM, CPM, TPM),
#'   with rows as features (e.g., genes) and columns as samples.
#' @param metadata A data frame containing sample metadata. Row names of \code{metadata} must correspond to
#'   the column names of \code{counts} (sample IDs). The function will reorder \code{metadata} to match \code{counts}.
#' @param title A character string for the PCA plot title.
#' @param colby A character string specifying the column name in \code{metadata} used to color samples in the PCA plot.
#'   Default is \code{"samples"}.
#' @param shapekey Optional character string specifying the column name in \code{metadata} used to shape points in the PCA plot.
#'   Default is \code{NULL}, meaning no grouping by shape.
#' @param shape_map Optional named vector mapping levels of \code{shapekey} to point shape codes (integers).
#'   For example: \code{c("batch1" = 15, "batch2" = 17)}. Default is \code{NULL} (PCAtools chooses shapes automatically).
#' @param showLoadings Logical indicating whether to display loadings on the PCA biplot. Default is \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the PCA biplot.
#'
#' @details
#' The function computes PCA with \code{PCAtools::pca} and then plots using \code{PCAtools::biplot}.
#' It performs basic input validation:
#' \itemize{
#'   \item Ensures \code{counts} has column names and \code{metadata} has row names.
#'   \item Reorders \code{metadata} rows to match \code{counts} columns.
#'   \item Checks that \code{colby} and (if provided) \code{shapekey} exist in \code{metadata}.
#' }
#'
#' @examples
#' \dontrun{
#' # Assume counts is your count matrix and metadata contains corresponding sample information:
#' p1 <- pcafun(counts, metadata,
#'              title = "PCA",
#'              colby = "treatment",
#'              shapekey = "batch",
#'              showLoadings = FALSE)
#'
#' # With explicit mapping of batch levels to shapes:
#' p2 <- pcafun(counts, metadata,
#'              title = "PCA",
#'              colby = "treatment",
#'              shapekey = "batch",
#'              shape_map = c("batch1" = 15, "batch2" = 17, "batch3" = 8),
#'              showLoadings = FALSE)
#' }
#'
#' @import ggplot2
#' @export
pcafun <- function(counts,
                   metadata,
                   title,
                   colby = "samples",
                   shapekey = NULL,
                   shape_map = NULL,
                   showLoadings = TRUE) {
  
  #check dependencies eheh (can't fool me twice)
  if (!requireNamespace("PCAtools", quietly = TRUE)) {
    stop("PCAtools is required for pcafun(). Install via Bioconductor: BiocManager::install('PCAtools')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for pcafun(). Install via install.packages('ggplot2')")
  }
  
  #coerce counts to matrix and validate
  counts <- as.matrix(counts)
  if (!is.numeric(counts)) stop("`counts` must be numeric.")
  if (is.null(colnames(counts))) stop("`counts` must have colnames (sample IDs).")
  
  #validate metadata and align to counts
  if (!is.data.frame(metadata)) stop("`metadata` must be a data.frame.")
  if (is.null(rownames(metadata))) stop("`metadata` must have rownames (sample IDs).")
  
  #reorder metadata to match counts columns (will introduce NA rows if mismatch)
  missing_meta <- setdiff(colnames(counts), rownames(metadata))
  if (length(missing_meta) > 0) {
    stop("`metadata` is missing these samples (rownames) present in `counts` colnames: ",
         paste(missing_meta, collapse = ", "))
  }
  metadata <- metadata[colnames(counts), , drop = FALSE]
  
  #check aesthetics columns exist
  if (!(colby %in% colnames(metadata))) {
    stop("`colby` not found in metadata: ", colby)
  }
  if (!is.null(shapekey) && !(shapekey %in% colnames(metadata))) {
    stop("`shapekey` not found in metadata: ", shapekey)
  }
  
  #PCA
  p <- PCAtools::pca(counts, metadata = metadata)
  
  #Plot
  #NOTE: PCAtools expects `shape` = column name, and `shapekey` = mapping.
  plot <- PCAtools::biplot(
    p,
    colby = colby,
    shape = shapekey,
    shapekey = shape_map,
    legendPosition = "right",
    showLoadings = showLoadings,
    lab = NULL,
    title = title
  )
  
  print(plot)
  return(plot)
}