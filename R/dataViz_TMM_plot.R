#' Plot Differentially Expressed Genes from a Count Matrix
#'
#' This function creates a bar plot of expression levels for the top differentially expressed genes (DEGs)
#' from a normalized count matrix (e.g., TMM or TPM). The function filters the DEG data frame based on a specified
#' adjusted p-value threshold and log fold-change, then selects and plots the top \code{top_n} genes. Gene names are italicized.
#'
#' @param deg A data frame containing differential expression results, including at least a log fold-change column
#'   and an adjusted p-value column. Row names should correspond to gene symbols.
#' @param tmm A normalized count matrix (e.g., TMM or TPM) with genes as rows and samples as columns.
#' @param title A character string specifying the title for the bar plot.
#' @param pval_col A character string specifying the column name in \code{deg} that contains adjusted p-values.
#'   Default is \code{"adj.P.Val"}.
#' @param logFC_col A character string specifying the column name in \code{deg} that contains log fold-changes.
#'   Default is \code{"logFC"}.
#' @param top_n An integer specifying the number of top DEGs to include in the plot. Default is \code{30}.
#' @param direction A character string indicating whether to visualize up-regulated genes ("UP") or down-regulated genes ("DOWN").
#'   Default is \code{"UP"}. For up-regulated genes, only genes with \code{logFC > 1} are selected; for down-regulated genes,
#'   only genes with \code{logFC < -1} are selected.
#'
#' @return A ggplot2 object representing a bar plot of the selected genes with their expression levels.
#'
#' @details
#' The function first filters the \code{deg} data frame to retain statistically significant genes (adjusted p-value below 0.05).
#' It then sorts the genes by log fold-change based on the specified \code{direction} ("UP" for up-regulated or "DOWN" for down-regulated)
#' and selects the top \code{top_n} genes. Next, the corresponding rows are extracted from the count matrix (\code{tmm}),
#' and the data is reshaped into a long format for plotting with ggplot2. Finally, a horizontal bar plot is created with gene names
#' displayed in italics.
#'
#' @examples
#' \dontrun{
#' # Example usage with TMM-normalized data:
#' # Assume deg_df is a DEG results data frame and tmm_avg is a TMM-normalized count matrix.
#' plot_up <- deg_barplot(deg_df, tmm_avg, "Top 30 UP Genes", direction = "UP")
#'
#' # For TPM-normalized data, the same function can be used:
#' plot_down <- deg_barplot(deg_df, tpm_avg, "Top 30 DOWN Genes", direction = "DOWN")
#' }
#'
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @export
deg_barplot <- function(deg, tmm, title, pval_col = "adj.P.Val", logFC_col = "logFC", top_n = 30, direction = "UP") {
  #filtering for statistically significant DEGs
  deg <- deg %>% filter(.data[[pval_col]] <= 0.05)

  #sorting DEGs based on logFC: UP-regulated or DOWN-regulated
  if (direction == "UP") {
    deg <- deg %>% filter(.data[[logFC_col]] > 1) %>% arrange(desc(.data[[logFC_col]]))  # Top UP genes
  } else if (direction == "DOWN") {
    deg <- deg %>% filter(.data[[logFC_col]] < -1) %>% arrange(.data[[logFC_col]])  # Top DOWN genes
  } else {
    stop("Invalid 'direction' argument. Use 'UP' or 'DOWN'.")
  }

  #selecting the top N genes
  deg <- head(deg, top_n)

  #selecting genes from the count matrix (TMM or TPM)
  selected_genes <- rownames(deg)
  tmm_filtered <- tmm[rownames(tmm) %in% selected_genes, , drop = FALSE]

  #reshaping data for ggplot2
  tmm_filtered <- data.frame(gene = rownames(tmm_filtered), tmm_filtered)
  tmm_long <- melt(tmm_filtered, id.vars = "gene", variable.name = "sample", value.name = "expression")

  #creating the bar plot with italicized gene names
  p <- ggplot(tmm_long, aes(x = gene, y = expression, fill = sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    labs(
      title = title,
      x = "Genes",
      y = "Expression level (TMM/TPM)",
      fill = "Samples"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 6, face = "italic")
    )

  print(p)
  return(p)
}
