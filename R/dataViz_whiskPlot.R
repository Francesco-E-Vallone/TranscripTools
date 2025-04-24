#' Whiskyplot: Faceted Boxplots with Gene Expression and Pairwise Statistics
#'
#' Subsets a normalized count matrix (TMM or TPM) to a set of genes, merges with sample metadata,
#' and draws faceted boxplots (“whisker plots”) with jittered points. Pairwise comparisons are
#' computed via **rstatix** and annotated by **ggpubr**.
#'
#' @param tmm A numeric matrix or data frame (rows = genes, columns = samples).
#' @param genes Character vector of gene names (must match \code{rownames(tmm)}).
#' @param metadata A data frame of sample metadata; must contain a column matching \code{sample_col},
#'   and row names matching \code{colnames(tmm)}.
#' @param sample_col Character; metadata column containing sample IDs (default \code{"samples"}).
#' @param group_col Character; metadata column defining comparison groups (default \code{"samples"}).
#' @param test Character; which test to run—either \code{"wilcox"} or \code{"t"}. Default \code{"wilcox"}.
#' @param p.adjust.method Character; p-value adjustment method for **rstatix** (default \code{"none"}).
#' @param comparisons Optional list of pairwise comparisons (each a length-2 character vector).  
#'   If \code{NULL}, all pairwise combinations are tested.
#' @param plot_title Character; overall plot title (default \code{"Gene expression"}).
#' @param x_lab Character; x-axis label (default = \code{group_col}).
#' @param y_lab Character; y-axis label (default \code{"log(TMM/TPM + 1)"}).
#'
#' @return A \code{ggplot} object: faceted whisker plots annotated with significant p-values.
#'
#' @details
#' 1. Verifies \code{group_col} exists in \code{metadata}.  
#' 2. Subsets \code{tmm} to \code{genes}, reshapes to long form.  
#' 3. Merges in \code{metadata}, drops samples missing \code{group_col}.  
#' 4. Keeps only genes present in ≥2 groups.  
#' 5. Runs pairwise tests per gene via **rstatix**, filters p ≤ 0.05.  
#' 6. Computes safe \code{y.position} for each gene’s annotation.  
#' 7. Draws boxplots + jitter; adds p-value annotations with **ggpubr**.
#'
#' @examples
#' \dontrun{
#'   comps <- list(c("NM_IgM","NWT_IgM"), c("NM_NT","NWT_NT"))
#'   whiskyplot(
#'     tmm             = tmm,
#'     genes           = rownames(tmm_allint),
#'     metadata        = meta,
#'     sample_col      = "SampleID",
#'     group_col       = "sample_complete",
#'     test            = "wilcox",
#'     p.adjust.method = "BH",
#'     comparisons     = comps,
#'     plot_title      = "Whiskyplot: Gene Expression",
#'     x_lab           = "Experimental Group",
#'     y_lab           = "Log(TPM + 1)"
#'   )
#' }
#'
#' @importFrom rlang .data
#' @importFrom rstatix pairwise_wilcox_test pairwise_t_test
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom dplyr filter group_by summarise mutate left_join n_distinct
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter facet_wrap 
#'   scale_y_continuous expansion theme_bw theme labs
#' @export
whiskyplot <- function(tmm, genes, metadata,
                       sample_col      = "samples",
                       group_col       = "samples",
                       test            = c("wilcox", "t"),
                       p.adjust.method = "none",
                       comparisons     = NULL,
                       plot_title      = "Gene expression",
                       x_lab           = group_col,
                       y_lab           = "log(TPM + 1)") {
  test <- match.arg(test)
  
  # 1. Validate group_col
  if (!group_col %in% colnames(metadata)) {
    stop("Column '", group_col, "' not found in metadata.")
  }
  
  # 2. Subset & reshape
  tmm_subset <- tmm[rownames(tmm) %in% genes, , drop = FALSE]
  df_long <- as.data.frame(tmm_subset) %>%
    tibble::rownames_to_column("Gene") %>%
    tidyr::pivot_longer(-Gene, names_to = "Sample", values_to = "Expression")
  
  # 3. Prepare metadata & merge
  metadata2 <- metadata %>% dplyr::mutate(Sample = .data[[sample_col]])
  df_long <- df_long %>%
    dplyr::left_join(metadata2, by = "Sample") %>%
    dplyr::filter(!is.na(.data[[group_col]])) %>%
    dplyr::mutate(!!group_col := as.factor(.data[[group_col]]))
  
  # 4. Keep only genes present in ≥2 groups
  valid_genes <- df_long %>%
    dplyr::group_by(Gene) %>%
    dplyr::filter(dplyr::n_distinct(.data[[group_col]]) >= 2) %>%
    dplyr::pull(Gene) %>%
    unique()
  df_long <- df_long %>% dplyr::filter(Gene %in% valid_genes)
  if (nrow(df_long) == 0) {
    warning("No genes have expression in two or more groups.")
    return(NULL)
  }
  
  # 5. Pairwise testing
  stat_fun <- switch(test,
                     wilcox = rstatix::pairwise_wilcox_test,
                     t      = rstatix::pairwise_t_test)
  stat_df <- df_long %>%
    dplyr::group_by(Gene) %>%
    stat_fun(
      formula         = as.formula(paste("Expression ~", group_col)),
      p.adjust.method = p.adjust.method,
      comparisons     = comparisons
    ) %>%
    dplyr::filter(p <= 0.05)
  
  # 6. Compute y.position
  if (nrow(stat_df) > 0) {
    y_df <- df_long %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(y.position = max(Expression, na.rm = TRUE) * 1.1, .groups="drop")
    stat_df <- stat_df %>% dplyr::left_join(y_df, by="Gene")
  }
  
  # 7. Plot
  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[[group_col]], y = Expression, fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    ggplot2::geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    ggplot2::facet_wrap(~ Gene, scales = "free_y") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = plot_title, x = x_lab, y = y_lab)
  
  if (nrow(stat_df) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(stat_df, label = "p", size = 2.5, inherit.aes = FALSE)
  }
  
  p
}