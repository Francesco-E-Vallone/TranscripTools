#' Add Significance Categories to Differential Expression Results
#'
#' This function adds a significance category to a differential expression results data frame.
#' The significance is determined using an FDR (or adjusted p-value) column and a log2 fold change column.
#' By default, genes with FDR < 0.05 and absolute logFC > 1 are categorized as "Up" (if logFC > 1) or "Down" (if logFC < -1).
#'
#' @param res A data frame of differential expression results, with row names representing gene identifiers.
#' @param fdr_col A character string specifying the name of the column that contains FDR or adjusted p-values. Default is \code{"FDR"}.
#' @param logFC_col A character string specifying the name of the column that contains log2 fold changes. Default is \code{"logFC"}.
#' @param fdr_threshold Numeric, threshold for statistical significance. Default is \code{0.05}.
#' @param fc_threshold Numeric, threshold for absolute log2 fold change. Default is \code{1}.
#'
#' @return The original data frame with an added column \code{sign} that categorizes genes as "Up", "Down", or "Not Significant".
#'
#' @details
#' Genes are categorized as follows:
#' \itemize{
#'   \item \strong{Up}: genes with \code{logFC} > fc_threshold and FDR < fdr_threshold.
#'   \item \strong{Down}: genes with \code{logFC} < -fc_threshold and FDR < fdr_threshold.
#'   \item \strong{Not Significant}: all other genes.
#' }
#'
#' @examples
#' \dontrun{
#' # Assume deg_df is your differential expression results data frame:
#' deg_df <- add_significance(deg_df, fdr_col = "adj.P.Val", logFC_col = "logFC",
#'                            fdr_threshold = 0.05, fc_threshold = 1)
#' }
#'
#' @export
add_significance <- function(res, fdr_col = "FDR", logFC_col = "logFC", fdr_threshold = 0.05, fc_threshold = 1) {
  res$sign <- "Not Significant"
  res$sign[res[[fdr_col]] < fdr_threshold & res[[logFC_col]] > fc_threshold] <- "Up"
  res$sign[res[[fdr_col]] < fdr_threshold & res[[logFC_col]] < -fc_threshold] <- "Down"
  return(res)
}


#' Draw a Volcano Plot for Differential Expression Results
#'
#' This function creates a volcano plot from a differential expression results data frame.
#' It uses \code{ggplot2} and \code{ggrepel} to plot log2 fold changes versus -log10(p-value),
#' coloring points by their significance category. It also labels the top \code{top_n} up‑ and down‑regulated genes.
#'
#' @param res A data frame of differential expression results that includes columns for log fold change (\code{logFC}),
#'   p-values (\code{PValue}), and significance category (\code{sign}). It is recommended to run \code{add_significance()} first.
#' @param fdr_threshold Numeric, the FDR threshold used for drawing the horizontal reference line. Default is \code{0.05}.
#' @param fc_threshold Numeric, the log2 fold change threshold used for drawing vertical reference lines. Default is \code{1}.
#' @param colors A named vector of colors for the significance categories. Default is \code{c("Up" = "#FF8F8F", "Down" = "#30B3A9", "Not Significant" = "grey")}.
#' @param top_n An integer specifying the number of top up‑regulated and top down‑regulated genes to label on the plot. Default is \code{10}.
#'
#' @return A \code{ggplot} object representing the volcano plot.
#'
#' @details
#' The function identifies the top \code{top_n} up‑regulated genes (highest \code{logFC}) and the top \code{top_n} down‑regulated genes
#' (lowest \code{logFC}) from those with p-values less than the threshold. These genes are then labeled on the plot using \code{geom_text_repel}.
#' The x-axis displays log2 fold change and the y-axis shows -log10(p-value). Dashed lines indicate the significance thresholds.
#' Important: run the add_significance() function prior to this one, as it is essential to have the labels of logFC direction ("Up", "Down", "Not Significant")
#'
#' @examples
#' \dontrun{
#' # Assume deg_df is your differential expression results data frame and has been processed with add_significance()
#' volcano_plot <- volplot(deg_df, fdr_threshold = 0.05, fc_threshold = 1, top_n = 10)
#' }
#'
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#' @export
volplot <- function(res, fdr_threshold = 0.05, fc_threshold = 1,
                    colors = c("Up" = "#FF8F8F", "Down" = "#30B3A9", "Not Significant" = "grey"),
                    top_n = 10) {
  #select top up-regulated genes
  top_up <- res %>%
    filter(PValue < fdr_threshold, logFC > fc_threshold) %>%
    arrange(desc(logFC)) %>%
    slice_head(n = top_n)

  #select top down-regulated genes
  top_down <- res %>%
    filter(PValue < fdr_threshold, logFC < -fc_threshold) %>%
    arrange(logFC) %>%
    slice_head(n = top_n)

  #combine top genes for labeling
  top_genes <- bind_rows(top_up, top_down)

  #ensure gene names are available for labeling
  res$Gene <- rownames(res)
  top_genes$label <- rownames(top_genes)

  #create the volcano plot
  p <- ggplot(res, aes(x = logFC, y = -log10(PValue), color = sign)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed") +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    theme_minimal() +
    labs(title = "Volcano Plot",
         x = "Log2 Fold Change",
         y = "-log10(p-value)") +
    theme(legend.position = "right",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
    geom_text_repel(
      data = top_genes,
      aes(label = label),
      size = 2,
      color = "black",
      fontface = "italic",
      max.overlaps = Inf,
      box.padding = 0.5
    ) +
    scale_x_continuous(limits = c(min(res$logFC) - 1, max(res$logFC) + 1)) +
    scale_y_continuous(trans = "log1p")

  print(p)
  return(p)
}
