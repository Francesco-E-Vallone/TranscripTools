#' Add Significance Categories to Differential Expression Results
#'
#' This function adds a significance category to a differential expression results data frame.
#' The significance is determined using an FDR (or adjusted p-value) column and a log2 fold change column.
#' By default, genes with FDR < 0.05 and absolute logFC > 1 are categorized as "Up" (if logFC > 1) or "Down" (if logFC < -1).
#'
#' Flag differential expression direction by thresholds
#'
#' Classifies each row as "Up", "Down", or "Not Significant" using a chosen
#' adjusted p-value column and a log2 fold-change column.
#'
#' @param res data.frame of DE results (rows = genes).
#' @param fdr_col Character. Column name with adjusted p-values (e.g., FDR, padj). Default "FDR".
#' @param logFC_col Character. Column name with log2 fold change. Default "logFC".
#' @param fdr_threshold Numeric. Adjusted p-value cutoff. Default 0.05.
#' @param fc_threshold Numeric. Absolute log2FC cutoff. Default 1.
#'
#' @return `res` with an added factor column `sign` in c("Up","Down","Not Significant").
#' @details
#' Up: logFC >  fc_threshold & FDR < fdr_threshold  
#' Down: logFC < -fc_threshold & FDR < fdr_threshold
#' @examples
#' # df <- add_significance(df, fdr_col="adj.P.Val", logFC_col="logFC", fdr_threshold=0.05, fc_threshold=1)
#' @export
add_significance <- function(res,
                             fdr_col = "FDR",
                             logFC_col = "logFC",
                             fdr_threshold = 0.05,
                             fc_threshold = 1) {
  stopifnot(is.data.frame(res))
  need <- c(fdr_col, logFC_col)
  miss <- setdiff(need, names(res))
  if (length(miss)) stop("Missing required column(s): ", paste(miss, collapse = ", "))
  
  fdr <- suppressWarnings(as.numeric(res[[fdr_col]]))
  lfc <- suppressWarnings(as.numeric(res[[logFC_col]]))
  
  sign <- rep("Not Significant", nrow(res))
  ok <- !is.na(fdr) & !is.na(lfc)
  sign[ ok & fdr < fdr_threshold & lfc >  fc_threshold] <- "Up"
  sign[ ok & fdr < fdr_threshold & lfc < -fc_threshold] <- "Down"
  
  res$sign <- factor(sign, levels = c("Up","Down","Not Significant"))
  res
}


#' Draw a Volcano Plot for Differential Expression Results
#'
#' This function creates a volcano plot from a differential expression results data frame.
#' It uses \code{ggplot2} and \code{ggrepel} to plot log2 fold changes versus -log10(p-value),
#' coloring points by their significance category. It also labels the top \code{top_n} up‑ and down‑regulated genes.
#'
#' Volcano plot (custom p-value column + fully customizable labels)
#'
#' Draws log2FC (x) vs -log10(p) (y). Lets you declare which column is the p-value,
#' and customize all plot labels. If `sign` is absent, it will be computed using
#' `add_significance()` with `fdr_col` when available (else it falls back to `p_col`).
#'
#' @param res data.frame with at least `logFC_col` and `p_col`.
#' @param p_col Character. Column name for p-values used on y-axis. Default "PValue".
#' @param logFC_col Character. Column name for log2 fold change. Default "logFC".
#' @param fdr_col Character. Adjusted p-value column (used to derive `sign` if missing). Default "FDR".
#' @param fdr_threshold Numeric. Threshold for horizontal line and significance logic. Default 0.05.
#' @param fc_threshold Numeric. Threshold for vertical lines and significance logic. Default 1.
#' @param colors Named vector mapping c("Up","Down","Not Significant") to colors.
#' @param top_n Integer. Label top N up and top N down genes among significant rows. Default 10.
#' @param gene_col Optional character: column with gene names; if NULL uses rownames.
#' @param title,subtitle,caption,x_lab,y_lab,legend_title Character plot labels.
#'
#' @return A ggplot object (returned invisibly).
#' @examples
#' # p <- volplot(df, p_col="pval", logFC_col="log2FC", fdr_col="padj",
#' #              title="My volcano", x_lab="log2FC", y_lab="-log10(p)")
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#' @importFrom rlang .data
#' @export
volplot <- function(res,
                    p_col = "PValue",
                    logFC_col = "logFC",
                    fdr_col = "FDR",
                    fdr_threshold = 0.05,
                    fc_threshold = 1,
                    colors = c("Up" = "#FF8F8F", "Down" = "#30B3A9", "Not Significant" = "grey"),
                    top_n = 10,
                    gene_col = NULL,
                    title = "Volcano plot",
                    subtitle = NULL,
                    caption = NULL,
                    x_lab = "Log2 fold change",
                    y_lab = "-log10(p-value)",
                    legend_title = NULL) {
  
  stopifnot(is.data.frame(res))
  need <- c(p_col, logFC_col)
  miss <- setdiff(need, names(res))
  if (length(miss)) stop("Missing required column(s): ", paste(miss, collapse = ", "))
  
  res[[p_col]]     <- suppressWarnings(as.numeric(res[[p_col]]))
  res[[logFC_col]] <- suppressWarnings(as.numeric(res[[logFC_col]]))
  
  if (!("sign" %in% names(res))) {
    if (fdr_col %in% names(res)) {
      res <- add_significance(res, fdr_col = fdr_col, logFC_col = logFC_col,
                              fdr_threshold = fdr_threshold, fc_threshold = fc_threshold)
    } else {
      res <- add_significance(res, fdr_col = p_col, logFC_col = logFC_col,
                              fdr_threshold = fdr_threshold, fc_threshold = fc_threshold)
    }
  }
  
  if (is.null(gene_col)) {
    res$..gene_label.. <- if (!is.null(rownames(res))) rownames(res) else as.character(seq_len(nrow(res)))
    gene_col <- "..gene_label.."
  } else if (!gene_col %in% names(res)) {
    stop("gene_col '", gene_col, "' not found in data.")
  }
  
  #direction
  sig <- !is.na(res[[p_col]]) & res[[p_col]] < fdr_threshold
  up_df   <- dplyr::filter(res, sig, .data[[logFC_col]] >  fc_threshold) |> dplyr::arrange(dplyr::desc(.data[[logFC_col]])) |> dplyr::slice_head(n = top_n)
  down_df <- dplyr::filter(res, sig, .data[[logFC_col]] < -fc_threshold) |> dplyr::arrange(.data[[logFC_col]]) |> dplyr::slice_head(n = top_n)
  top_df  <- dplyr::bind_rows(up_df, down_df)
  
  #plot
  p <- ggplot2::ggplot(res, ggplot2::aes(x = .data[[logFC_col]], y = -log10(.data[[p_col]]), color = .data$sign)) +
    ggplot2::geom_point(alpha = 0.8, size = 2) +
    ggplot2::scale_color_manual(values = colors, name = legend_title) +
    ggplot2::geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = caption, x = x_lab, y = y_lab) +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  if (nrow(top_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top_df,
      ggplot2::aes(label = .data[[gene_col]]),
      size = 3,
      color = "black",
      max.overlaps = Inf,
      box.padding = 0.5
    )
  }
  
  print(p)
  return(p)
}