#' Plot Gene Expression Statistics with Custom Labels and Statistical Comparisons
#'
#' This function creates a faceted boxplot (with overlayed jittered points) to display gene expression 
#' levels for a set of genes across samples. It is designed to work with normalized count matrices such as 
#' TMM or TPM and allows you to add pairwise statistical comparisons (e.g., via "t.test", "wilcox.test", "anova"). 
#' Additionally, users can customize plot title and axis labels.
#'
#' @param tmm Numeric matrix/data.frame (rows = genes, cols = samples).
#' @param genes Character vector of genes to plot.
#' @param metadata Data frame with sample metadata.
#' @param sample_col Character. Column in `metadata` containing sample IDs matching `colnames(tmm)`. Default "samples".
#' @param group_col Character. Column in `metadata` defining groups. Default "samples".
#' @param test Character. Statistical test for pairwise comparisons (e.g., "t.test", "wilcox.test", "anova"). Default "t.test".
#' @param comparisons Optional list of length-2 character vectors with group names for pairwise tests.
#' @param log2p1 Logical. If TRUE, plot log2(Expression + 1). Default TRUE.
#' @param palette Optional named vector for group fill colors. Names must match group levels.
#' @param facet_scales Character, passed to `facet_wrap(scales=)`. Default "free_y".
#' @param plot_title,x_lab,y_lab,subtitle,caption Character plot labels. `x_lab` defaults to `group_col`; `y_lab` auto-set if `log2p1=TRUE`.
#'
#' @return A ggplot object (invisible).
#' @import dplyr tibble tidyr ggplot2 ggpubr
#' @export
plot_genes_stats <- function(tmm, genes, metadata,
                             sample_col = "samples",
                             group_col  = "samples",
                             test = "t.test",
                             comparisons = NULL,
                             log2p1 = TRUE,
                             palette = NULL,
                             facet_scales = "free_y",
                             plot_title = "Gene Expression",
                             x_lab = group_col,
                             y_lab = if (log2p1) "log2(Expression + 1)" else "Expression",
                             subtitle = NULL,
                             caption  = NULL) {
  
  #subset to requested genes
  tmm_sub <- tmm[rownames(tmm) %in% genes, , drop = FALSE]
  
  #long format
  df_long <- as.data.frame(tmm_sub) |>
    tibble::rownames_to_column("Gene") |>
    tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")
  
  #check metadata columns
  if (!all(c(sample_col, group_col) %in% colnames(metadata))) {
    stop("`metadata` must contain columns: '", sample_col, "' and '", group_col, "'.")
  }
  
  #join metadata
  df_long <- df_long |>
    dplyr::left_join(metadata |>
                       dplyr::select(dplyr::all_of(c(sample_col, group_col))) |>
                       dplyr::rename(Sample = !!sample_col, .group = !!group_col),
                     by = "Sample") |>
    dplyr::filter(!is.na(.group))
  
  #transform if requested
  if (log2p1) df_long$Expression <- log2(df_long$Expression + 1)
  
  #factor for consistent ordering
  df_long$.group <- as.factor(df_long$.group)
  
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .group, y = Expression, fill = .group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    ggplot2::facet_wrap(~ Gene, scales = facet_scales) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = plot_title, subtitle = subtitle, caption = caption,
                  x = x_lab, y = y_lab, fill = NULL)
  
  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }
  
  if (!is.null(comparisons)) {
    p <- p + ggpubr::stat_compare_means(
      comparisons = comparisons,
      method = test,
      label = "p.format",
      group.by = "Gene",
      hide.ns = FALSE
    )
  }
  
  print(p)
  return(p)
}