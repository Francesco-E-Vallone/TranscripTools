#' Generate Enrichment Plots from GO Data (from the enrichR package)
#'
#' This function generates a list of enrichment plots for different databases
#' contained in the input list. Each element of the list \code{GO} is processed
#' to create a plot using the \code{plotEnrich} function.
#'
#' @param GO A named list where each element contains data (e.g., a data frame)
#'   representing GO (Gene Ontology) enrichment results for a specific database.
#'
#' @return A named list of plots. Each plot corresponds to one of the databases
#'   in the input \code{GO} list.
#'
#' @import enrichR
#'
#' @details
#' The function loops through each database in the \code{GO} list, converts
#' the corresponding element to a data frame (if it is not already one), and
#' generates a plot using the \code{plotEnrich} function.
#'
#' @examples
#' # Example usage with a hypothetical GO list:
#' # GO <- list(
#' #   KEGG = kegg_data,
#' #   Reactome = reactome_data
#' # )
#' # plots <- enrichPlot(GO)
#' @export
enrichPlot <- function(GO) {
  plots <- list()
  for (db in names(GO)) {
    df <- as.data.frame(GO[[db]])
    p <- enrichR::plotEnrich(df)
    plots[[db]] <- p
  }
  return(plots)
}

#' Perform GO Analysis on Up-regulated Genes
#'
#' This function performs Gene Ontology (GO) analysis on up-regulated genes from a differential
#' expression data frame. By default, it uses a preset list of GO and pathway databases. However,
#' you can supply your own database vector via the \code{databases} parameter. The function also
#' creates bar plots for the top 20 GO terms (based on p-value) for each database.
#'
#' @param df Data frame of DE results with rownames = symbols and a `logFC` column.
#' @param samples Character used in plot titles.
#' @param databases Character vector of Enrichr databases to use (any length).
#'
#' @return Named list of ggplot objects (one per database).
#' @import enrichR
#' @import ggplot2
#' @export
up_go <- function(df, samples,
                  databases = c("MSigDB_Hallmark_2020",
                                "GO_Biological_Process_2023",
                                "BioPlanet_2019",
                                "GO_Cellular_Component_2023",
                                "Reactome_Pathways_2024")) {
  if (!requireNamespace("enrichR", quietly = T)) stop("Package 'enrichR' is required")
  if (!requireNamespace("ggplot2", quietly = T)) stop("Package 'ggplot2' is required")
  
  #select UP genes
  df <- df[df$logFC >= 1, , drop = F]
  genes <- rownames(df)
  if (!length(genes)) {
    warning("No up-regulated genes (logFC >= 1).")
    return(list())
  }
  
  #run Enrichr
  go_list <- enrichR::enrichr(genes, databases = databases)
  
  plots <- list()
  for (db in databases) {
    if (!db %in% names(go_list) || is.null(go_list[[db]]) || nrow(go_list[[db]]) == 0) {
      warning("No results returned for database: ", db)
      next
    }
    
    df_db <- as.data.frame(go_list[[db]])
    
    #use column names, not fixed positions (robust to enrichR changes)
    if (!all(c("Term","P.value") %in% names(df_db))) next
    df_db <- df_db[df_db$`P.value` < 0.05, intersect(c("Term","P.value","Adjusted.P.value","Genes"), names(df_db)), drop = FALSE]
    if (!nrow(df_db)) next
    
    #top 20 by smallest p
    df_db <- df_db[order(df_db$`P.value`, decreasing = FALSE), , drop = FALSE]
    df_db <- head(df_db, 20)
    
    p <- ggplot2::ggplot(df_db, ggplot2::aes(x = stats::reorder(.data$Term, -.data$`P.value`),
                                             y = -log10(.data$`P.value`))) +
      ggplot2::geom_col(fill = "#8FE0BC", color = "black") +
      ggplot2::coord_flip() +
      ggplot2::theme_light() +
      ggplot2::labs(title = paste(db, "| Up-regulated genes:", samples),
                    x = "Terms", y = "-log10(P-value)") +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
    
    plots[[db]] <- p
  }
  
  return(plots)
}

#' Perform GO Analysis on Down-regulated Genes
#'
#' This function performs Gene Ontology (GO) analysis on down-regulated genes from a differential
#' expression data frame. It works similarly to \code{up_go()} but filters for genes with \code{logFC} <= -1.
#' By default, a preset list of databases is used, but you can supply your own via the \code{databases} parameter.
#' Bar plots are generated for the top 20 GO terms.
#'
#' Perform GO Analysis on Down-regulated Genes (flexible databases)
#'
#' @param df Data frame of DE results with rownames = symbols and a `logFC` column.
#' @param samples Character used in plot titles.
#' @param databases Character vector of Enrichr databases to use (any length).
#'
#' @return Named list of ggplot objects (one per database).
#' @import enrichR
#' @import ggplot2
#' @export
down_go <- function(df, samples,
                    databases = c("MSigDB_Hallmark_2020",
                                  "GO_Biological_Process_2023",
                                  "BioPlanet_2019",
                                  "GO_Cellular_Component_2023",
                                  "Reactome_Pathways_2024")) {
  if (!requireNamespace("enrichR", quietly = TRUE)) stop("Package 'enrichR' is required")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")
  
  #select DOWN genes
  df <- df[df$logFC <= -1, , drop = FALSE]
  genes <- rownames(df)
  if (!length(genes)) {
    warning("No down-regulated genes (logFC <= -1).")
    return(list())
  }
  
  #run Enrichr
  go_list <- enrichR::enrichr(genes, databases = databases)
  
  plots <- list()
  for (db in databases) {
    if (!db %in% names(go_list) || is.null(go_list[[db]]) || nrow(go_list[[db]]) == 0) {
      warning("No results returned for database: ", db)
      next
    }
    
    df_db <- as.data.frame(go_list[[db]])
    
    if (!all(c("Term","P.value") %in% names(df_db))) next
    df_db <- df_db[df_db$`P.value` < 0.05, intersect(c("Term","P.value","Adjusted.P.value","Genes"), names(df_db)), drop = FALSE]
    if (!nrow(df_db)) next
    
    df_db <- df_db[order(df_db$`P.value`, decreasing = FALSE), , drop = FALSE]
    df_db <- head(df_db, 20)
    
    p <- ggplot2::ggplot(df_db, ggplot2::aes(x = stats::reorder(.data$Term, -.data$`P.value`),
                                             y = -log10(.data$`P.value`))) +
      ggplot2::geom_col(fill = "#706F9B", color = "black") +
      ggplot2::coord_flip() +
      ggplot2::theme_light() +
      ggplot2::labs(title = paste(db, "| Down-regulated genes:", samples),
                    x = "Terms", y = "-log10(P-value)") +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
    
    plots[[db]] <- p
  }
  
  return(plots)
}

#' Save GO Analysis Plots to Files
#'
#' This function saves a list of ggplot2 objects to files in a chosen format. The file names are constructed from the provided
#' \code{list_name} and the names of each plot in the list.
#'
#' @param plot_list A named list of ggplot2 plot objects.
#' @param list_name A character string to prefix the file names.
#' @param path A character string specifying the directory path where the plots will be saved. Default is \code{"results/"}.
#' @param device A character string specifying the output device/format. Default is \code{"svg"}. Other possible values include \code{"png"}, \code{"pdf"}, etc.
#' @param height A numeric value specifying the height of the saved plot (in inches). Default is 7.
#' @param width A numeric value specifying the width of the saved plot (in inches). Default is 15.
#'
#' @return None. Files are written to disk.
#'
#' @details
#' The function iterates over each plot in \code{plot_list} and uses \code{ggsave} to save each plot in the specified format.
#'
#' @examples
#' \dontrun{
#'   # Assume up_plots is a list of ggplot objects from up_go()
#'   save_plot(up_plots, list_name = "up_GO", device = "png")
#' }
#'
#' @import ggplot2
#' @import svglite
#' @export
save_plot <- function(plot_list, list_name, path = "results/", device = "svg", height = 7, width = 15) {
  #check if dir exists
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  
  #looping for each name of the plot list
  for (i in names(plot_list)) {
    ggplot2::ggsave(filename = paste0(path, list_name, "_", i, ".", device),
                    plot = plot_list[[i]],
                    device = device,
                    height = height,
                    width = width)
  }
}

#' Build a Nested List of GO Analysis Results for Export
#'
#' Combines up- and down-regulated GO/enrichment results into a nested list suitable for export.
#' Missing comparisons (e.g. only UP available) and missing databases are allowed.
#'
#' Input can be either:
#' \itemize{
#'   \item Multi-comparison: \code{list(comparison -> list(database -> result))}
#'   \item Single-comparison: \code{list(database -> result)}
#' }
#'
#' Each database result can be:
#' \itemize{
#'   \item a data.frame (already extracted), or
#'   \item a list-like object containing \code{$data} (your original assumption).
#' }
#'
#' @param up_list Up-regulated results (single-comparison or multi-comparison list).
#' @param down_list Down-regulated results (single-comparison or multi-comparison list).
#' @param databases Character vector of database names to extract.
#' @param comparison_name Name to use if a single-comparison structure is supplied.
#' @param keep_empty Logical. If TRUE (default), missing db/direction becomes empty data.frame
#'   so structure is stable for Excel export. If FALSE, missing/empty entries are omitted.
#' @param warn Logical. If TRUE (default), warn when results are missing/malformed.
#'
#' @return A nested list where each top-level element is a \code{comparison} and
#' contains names like \code{paste0(db, "_UP")} and \code{paste0(db, "_Down")}.
#'
#' @export
build_go_list <- function(
    up_list,
    down_list,
    databases,
    comparison_name = "comparison1",
    keep_empty = TRUE,
    warn = TRUE
) {
  
  if (!is.list(up_list) || !is.list(down_list)) {
    stop("`up_list` and `down_list` must be lists.")
  }
  
  empty_df <- function() data.frame()
  
  #detect "single comparison" db-level lists by name overlap with databases
  is_single_comparison <- function(x, dbs) {
    is.list(x) && !is.null(names(x)) && any(names(x) %in% dbs)
  }
  
  #extract a data.frame from a db entry:
  # - if already a data.frame -> use it
  # - else if has $data -> as.data.frame($data)
  # - else -> missing
  extract_df <- function(obj, comp, db, suffix) {
    label <- paste0("'", comp, "' / '", db, "_", suffix, "'")
    
    if (is.null(obj)) {
      if (warn) warning("Missing result for ", label)
      return(if (keep_empty) empty_df() else NULL)
    }
    
    if (is.data.frame(obj)) {
      return(if (nrow(obj) == 0 && keep_empty) empty_df() else if (nrow(obj) == 0) NULL else obj)
    }
    
    #accept tibbles etc.
    if (inherits(obj, "tbl_df")) {
      obj <- as.data.frame(obj)
      return(if (nrow(obj) == 0 && keep_empty) empty_df() else if (nrow(obj) == 0) NULL else obj)
    }
    
    #obj$data
    if (is.list(obj) && !is.null(obj$data)) {
      df <- as.data.frame(obj$data)
      if (is.null(df) || nrow(df) == 0) return(if (keep_empty) empty_df() else NULL)
      return(df)
    }
    
    if (warn) warning("Unrecognized object (no data.frame, no $data) for ", label)
    if (keep_empty) empty_df() else NULL
  }
  
  #normalize single-comparison inputs to multi-comparison form
  if (is_single_comparison(up_list, databases)) {
    up_list <- setNames(list(up_list), comparison_name)
  }
  if (is_single_comparison(down_list, databases)) {
    down_list <- setNames(list(down_list), comparison_name)
  }
  
  up_names <- names(up_list); if (is.null(up_names)) up_names <- character(0)
  dn_names <- names(down_list); if (is.null(dn_names)) dn_names <- character(0)
  
  #KEY CHANGE: UNION to avoid loosing UP-only or DOWN-only comparisons
  comparisons <- union(up_names, dn_names)
  
  #if both are totally unnamed at comparison level, the only safe rule is index pairing.
  if (length(comparisons) == 0) {
    n <- max(length(up_list), length(down_list))
    if (n == 0) stop("Empty `up_list` and `down_list`.")
    
    if (warn) warning("Both `up_list` and `down_list` are unnamed: pairing comparisons by index.")
    
    comparisons <- paste0("comparison", seq_len(n))
    out <- setNames(vector("list", n), comparisons)
    
    for (i in seq_len(n)) {
      comp <- comparisons[i]
      up_comp <- if (i <= length(up_list)) up_list[[i]] else NULL
      dn_comp <- if (i <= length(down_list)) down_list[[i]] else NULL
      
      comp_list <- list()
      for (db in databases) {
        up_df <- extract_df(if (!is.null(up_comp)) up_comp[[db]] else NULL, comp, db, "UP")
        if (!is.null(up_df)) comp_list[[paste0(db, "_UP")]] <- up_df
        
        dn_df <- extract_df(if (!is.null(dn_comp)) dn_comp[[db]] else NULL, comp, db, "Down")
        if (!is.null(dn_df)) comp_list[[paste0(db, "_Down")]] <- dn_df
      }
      
      if (!keep_empty) comp_list <- comp_list[vapply(comp_list, nrow, integer(1)) > 0]
      out[[comp]] <- comp_list
    }
    
    return(out)
  }
  
  out <- setNames(vector("list", length(comparisons)), comparisons)
  
  for (comp in comparisons) {
    up_comp <- up_list[[comp]]
    dn_comp <- down_list[[comp]]
    
    if (is.null(up_comp) && warn) warning("Comparison '", comp, "' missing in `up_list`")
    if (is.null(dn_comp) && warn) warning("Comparison '", comp, "' missing in `down_list`")
    
    comp_list <- list()
    
    for (db in databases) {
      up_df <- extract_df(if (!is.null(up_comp)) up_comp[[db]] else NULL, comp, db, "UP")
      if (!is.null(up_df)) comp_list[[paste0(db, "_UP")]] <- up_df
      
      dn_df <- extract_df(if (!is.null(dn_comp)) dn_comp[[db]] else NULL, comp, db, "Down")
      if (!is.null(dn_df)) comp_list[[paste0(db, "_Down")]] <- dn_df
    }
    
    if (!keep_empty) comp_list <- comp_list[vapply(comp_list, nrow, integer(1)) > 0]
    out[[comp]] <- comp_list
  }
  
  out
}


#' Combine and Prepare GO Analysis Results for Export
#'
#' This function takes a nested list of GO analysis results and, for each comparison/category,
#' combines the contained data frames into a single data frame by row-binding them.
#' A new column \code{Database} is added (from the names of the nested list elements) to indicate
#' the source database and direction (e.g. \code{"GO_Biological_Process_2025_UP"}).
#'
#' Unlike earlier versions, this function does \strong{not} overwrite names with a fixed template.
#' This makes it robust to missing databases (very common in GO analyses) and prevents silent mislabeling.
#'
#' @param go_list A named list of GO analysis results. Each element of \code{go_list} should itself be a named list
#'   of data frames (or tibbles). Names should indicate database and direction (e.g. \code{"Reactome_Pathways_2024_Down"}).
#'
#' @return A named list of data frames. Each data frame contains an additional column \code{Database} indicating
#'   the source of each row (taken from the nested list names).
#'
#' @examples
#' \dontrun{
#' combined_go_list <- prep_go_exp(go_list)
#' writexl::write_xlsx(combined_go_list, "GO_Analysis_Results.xlsx")
#' }
#'
#' @export
prep_go_exp <- function(go_list) {
  
  #dependency check
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for prep_go_exp(). Install it with: install.packages('dplyr')")
  }
  
  if (!is.list(go_list) || is.null(names(go_list))) {
    stop("`go_list` must be a named list (e.g., comparisons/categories as names).")
  }
  
  #combined list
  combined_go_list <- lapply(go_list, function(sublist) {
    
    if (!is.list(sublist)) {
      return(data.frame())
    }
    
    #if sublist has no names, create stable fallback names
    if (is.null(names(sublist))) {
      names(sublist) <- paste0("db", seq_along(sublist))
    } else {
      # Replace empty names with fallbacks (rare but happens)
      empty <- which(is.na(names(sublist)) | names(sublist) == "")
      if (length(empty) > 0) names(sublist)[empty] <- paste0("db", empty)
    }
    
    #coerce each element to data.frame when possible
    sublist <- lapply(sublist, function(x) {
      if (is.null(x)) return(data.frame())
      if (inherits(x, "tbl_df")) return(as.data.frame(x))
      if (is.data.frame(x)) return(x)
      data.frame()
    })
    
    dplyr::bind_rows(sublist, .id = "Database")
  })
  
  return(combined_go_list)
}
