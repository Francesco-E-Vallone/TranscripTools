#' Perform and Plot PCA with Customizable Aesthetics
#'
#' This function performs Principal Component Analysis (PCA) on a count matrix and generates a biplot
#' using \code{PCAtools}. It wraps the native PCA tools functionality to streamline the workflow,
#' especially for pathway-specific PCAs, by performing the analysis and plotting in one step.
#'
#' @param counts A numeric matrix or data frame containing count data (e.g., TMM or TPM values) with rows as features (e.g., genes) and columns as samples.
#' @param metadata A data frame containing sample metadata. The row names of \code{metadata} should correspond to the column names of \code{counts}.
#' @param title A character string for the plot title.
#' @param colby A character string specifying the column name in \code{metadata} used to color the samples in the PCA plot.
#'   Default is \code{"samples"}.
#' @param shapeby An optional character string specifying the column name in \code{metadata} used to shape the points in the PCA plot.
#'   Default is \code{NULL}, meaning no shaping.
#' @param showLoadings A logical value indicating whether to display the loadings on the PCA biplot.
#'   Default is \code{TRUE}.
#'
#' @return A \code{ggplot} object representing the PCA biplot.
#'
#' @details
#' This function leverages \code{PCAtools::pca} to compute the principal components and \code{PCAtools::biplot}
#' to create a biplot. By including customizable parameters for color and shape aesthetics, along with an option
#' to toggle the display of loadings, it allows users to quickly generate tailored PCA visualizations without having
#' to perform the PCA calculation and plot generation as separate steps.
#'
#' @examples
#' \dontrun{
#' # Assume counts is your count matrix and metadata contains corresponding sample information:
#' pca_plot <- pcafun(counts, metadata, title = "Pathway-Specific PCA",
#'                    colby = "treatment", shapeby = "batch", showLoadings = FALSE)
#' }
#'
#' @import PCAtools
#' @import ggplot2
#' @export
pcafun <- function(counts, metadata, title, colby = "samples", shapeby = NULL, showLoadings = TRUE) {
  #perform PCA using PCAtools::pca
  p <- PCAtools::pca(counts, metadata = metadata)

  #create a biplot using PCAtools::biplot with customizable aesthetics
  plot <- PCAtools::biplot(p,
                           colby = colby,
                           shapeby = shapeby,
                           legendPosition = "right",
                           showLoadings = showLoadings,
                           lab = NULL,
                           title = title)

  print(plot)
  return(plot)
}
