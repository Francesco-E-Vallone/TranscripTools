#' Plot Gene Expression Statistics with Custom Labels and Statistical Comparisons
#'
#' This function creates a faceted boxplot (with overlayed jittered points) to display gene expression 
#' levels for a set of genes across samples. It is designed to work with normalized count matrices such as 
#' TMM or TPM and allows you to add pairwise statistical comparisons (e.g., via "t.test", "wilcox.test", "anova"). 
#' Additionally, users can customize plot title and axis labels.
#'
#' @param tmm A numeric matrix or data frame of normalized counts (TMM or TPM) with row names representing genes 
#'   and column names representing sample IDs.
#' @param genes A character vector of gene names to include in the plot.
#' @param metadata A data frame containing sample metadata. Its row names should correspond to the column names in \code{tmm}.
#' @param sample_col A character string specifying the column in \code{metadata} that identifies samples. Default is \code{"samples"}.
#' @param group_col A character string specifying the column in \code{metadata} that defines groups for comparison. Default is \code{"samples"}.
#' @param test A character string specifying the statistical test to use for comparisons (e.g., \code{"t.test"}, \code{"wilcox.test"}, \code{"anova"}).
#'   Default is \code{"t.test"}.
#' @param comparisons An optional list of pairwise comparisons (each a character vector of two group names) to be used in 
#'   \code{ggpubr::stat_compare_means()}. If \code{NULL} (default), no statistical comparisons are added.
#' @param plot_title A character string for the overall plot title. Default is \code{"Gene Expression"}.
#' @param x_lab A character string for the x-axis label. Default is the value of \code{group_col}.
#' @param y_lab A character string for the y-axis label. Default is \code{"Log2(TMM/TPM + 1)"}.
#'
#' @return A \code{ggplot} object representing a faceted boxplot with gene expression data and optional statistical annotations.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Subsets the count matrix \code{tmm} to include only the genes specified in \code{genes}.
#'   \item Reshapes the resulting subset into long format using \code{tibble::rownames_to_column()} and \code{tidyr::pivot_longer()}.
#'   \item Merges the reshaped data with \code{metadata} based on sample IDs.
#'   \item Converts the grouping column (specified by \code{group_col}) to a factor.
#'   \item Creates a faceted boxplot (one facet per gene) displaying expression levels along with jittered data points.
#'   \item Optionally adds statistical comparison annotations to each facet using \code{ggpubr::stat_compare_means()}.
#' }
#' The function is designed for use with both TMM and TPM normalized data.
#'
#' @examples
#' \dontrun{
#'   # Example usage:
#'   # Assume 'tmm' is your normalized count matrix,
#'   # 'gene_vector' is a vector containing the genes to be plotted,
#'   # and 'meta' is a data frame of sample metadata.
#'   # You want to compare groups "NM_IgM" vs "NWT_IgM" and "NM_NT" vs "NWT_NT"
#'   plot_genes_stats(
#'     tmm = tmm, 
#'     genes = gene_vector, 
#'     metadata = meta, 
#'     group_col = "sample_complete",
#'     comparisons = list(c("NM_IgM", "NWT_IgM"), c("NM_NT", "NWT_NT")),
#'     test = "wilcox.test",
#'     plot_title = "Expression Levels of Selected Genes",
#'     x_lab = "Experimental Group",
#'     y_lab = "Log2(TMM/TPM + 1)"
#'   )
#' }
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import ggplot2
#' @import ggpubr
#' @export
plot_genes_stats <- function(tmm, genes, metadata,
                             sample_col = "samples",
                             group_col = "samples",
                             test = "t.test",
                             comparisons = NULL,
                             plot_title = "Gene Expression",
                             x_lab = group_col,
                             y_lab = "Log2(TMM/TPM + 1)") {
  #subset the count matrix to keep only the specified genes
  tmm_subset <- tmm[rownames(tmm) %in% genes, , drop = FALSE]
  
  #reshape the matrix from wide to long format and add a Gene column
  df_long <- as.data.frame(tmm_subset) %>%
    tibble::rownames_to_column("Gene") %>%
    tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")
  
  #convert metadata: move row names into a column named "Sample"
  metadata <- metadata %>% 
    tibble::rownames_to_column("Sample")
  
  #merge expression data with metadata and filter out samples with missing group information
  df_long <- df_long %>%
    left_join(metadata, by = "Sample") %>%
    filter(!is.na(.data[[group_col]]))
  
  #convert the grouping column to a factor
  df_long[[group_col]] <- as.factor(df_long[[group_col]])
  
  #create the boxplot with jittered points
  p <- ggplot(df_long, aes(x = .data[[group_col]], y = Expression, fill = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~ Gene, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = plot_title,
         x = x_lab,
         y = y_lab)
  
  #add statistical comparisons if provided
  if (!is.null(comparisons)) {
    p <- p + ggpubr::stat_compare_means(
      comparisons = comparisons,
      method = test,
      label = "p.format",
      group.by = "Gene",
      hide.ns = FALSE
    )
  }
  
  return(p)
}
