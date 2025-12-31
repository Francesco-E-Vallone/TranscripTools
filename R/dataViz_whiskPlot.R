#' Whiskyplot: Faceted Boxplots with Gene Expression and Pairwise Statistics
#'
#' Subsets a normalized count matrix (TMM or TPM) to a set of genes, merges with sample metadata,
#' and draws faceted boxplots (“whisker plots”) with jittered points. Pairwise comparisons are
#' computed via **rstatix** and annotated by **ggpubr**.
#'
#' @param tmm Numeric matrix/data.frame (rows = genes, cols = samples).
#' @param genes Character vector of gene names (must match rownames(tmm)).
#' @param metadata Data frame with sample metadata.
#' @param sample_col Character. Column in `metadata` containing sample IDs (matches colnames(tmm)). Default "samples".
#' @param group_col Character. Column in `metadata` defining groups. Default "samples".
#' @param test One of c("wilcox","t"). Default "wilcox".
#' @param p.adjust.method P-value adjustment for rstatix ("none","BH",...). Default "none".
#' @param comparisons Optional list of length-2 character vectors (group names). If NULL, all pairwise combos.
#' @param palette Optional named vector for group fill colors.
#' @param plot_title,x_lab,y_lab,subtitle,caption Character plot labels. `x_lab` defaults to `group_col`.
#'
#' @return A ggplot object.
#' 
#' @importFrom rlang .data
#' @importFrom rstatix pairwise_wilcox_test pairwise_t_test
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom dplyr filter group_by summarise mutate left_join n_distinct select rename
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter facet_wrap
#' @importFrom ggplot2 scale_y_continuous expansion theme_bw theme element_text labs
#' @export
whiskyplot <- function(tmm, genes, metadata,
                       sample_col      = "samples",
                       group_col       = "samples",
                       test            = c("wilcox", "t"),
                       p.adjust.method = "none",
                       comparisons     = NULL,
                       palette         = NULL,
                       plot_title      = "Gene expression",
                       x_lab           = group_col,
                       y_lab           = "log(TMM/TPM + 1)",
                       subtitle        = NULL,
                       caption         = NULL) {
  #dependencies check
  pkgs <- c("tibble","tidyr","dplyr","rstatix","ggplot2","ggpubr")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "))
  
  #match arg
  test <- match.arg(test)
  
  #subset & long format
  tmm_sub <- tmm[rownames(tmm) %in% genes, , drop = FALSE]
  df_long <- as.data.frame(tmm_sub) |>
    tibble::rownames_to_column("Gene") |>
    tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") |>
    dplyr::mutate(Expression = log1p(Expression)) #log-transformation
  
  #metadata join
  if (!all(c(sample_col, group_col) %in% colnames(metadata))) {
    stop("`metadata` must contain columns: '", sample_col, "' and '", group_col, "'.")
  }
  md <- metadata |>
    dplyr::select(dplyr::all_of(c(sample_col, group_col))) |>
    dplyr::rename(Sample = !!sample_col, .group = !!group_col)
  
  df_long <- df_long |>
    dplyr::left_join(md, by = "Sample") |>
    dplyr::filter(!is.na(.group)) |>
    dplyr::mutate(.group = as.factor(.group))
  
  #require ≥2 groups per gene
  keep_genes <- df_long |>
    dplyr::group_by(Gene) |>
    dplyr::filter(dplyr::n_distinct(.group) >= 2) |>
    dplyr::pull(Gene) |> unique()
  df_long <- dplyr::filter(df_long, Gene %in% keep_genes)
  if (nrow(df_long) == 0) {
    warning("No genes have data in two or more groups.")
    return(NULL)
  }
  
  #pairwise tests (per gene)
  stat_fun <- switch(test,
                     wilcox = rstatix::pairwise_wilcox_test,
                     t      = rstatix::pairwise_t_test
  )
  
  #if comparisons not provided, rstatix will compute all pairwise by default
  stat_df <- df_long |>
    dplyr::group_by(Gene) |>
    stat_fun(formula = stats::as.formula("Expression ~ .group"),
             p.adjust.method = p.adjust.method,
             comparisons = comparisons)
  
  #pick the right label column
  label_col <- if (identical(p.adjust.method, "none")) "p" else "p.adj"
  
  #y-position for p-value labels (per gene)
  if (nrow(stat_df) > 0) {
    y_df <- df_long |>
      dplyr::group_by(Gene) |>
      dplyr::summarise(y.position = max(Expression, na.rm = TRUE) * 1.08, .groups = "drop")
    stat_df <- dplyr::left_join(stat_df, y_df, by = "Gene")
  }
  
  #plot
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .group, y = Expression, fill = .group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    ggplot2::geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    ggplot2::facet_wrap(~ Gene, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(face = "bold")) +
    ggplot2::labs(title = plot_title, subtitle = subtitle, caption = caption,
                  x = x_lab, y = y_lab, fill = NULL)
  
  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_manual(values = palette)
  }
  
  if (nrow(stat_df) > 0) {
    #ggpubr expects columns: group1, group2, y.position, and a label column
    p <- p + ggpubr::stat_pvalue_manual(stat_df,
                                        label = label_col,
                                        size = 2.7,
                                        inherit.aes = FALSE)
  }
  
  print(p)
  return(p)
}