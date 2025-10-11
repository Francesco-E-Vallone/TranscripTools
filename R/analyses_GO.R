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
#' This function creates a nested list of GO analysis results for multiple comparisons by combining
#' up-regulated and down-regulated results. It expects two named lists: one with up-regulated GO results and one
#' with down-regulated GO results. For each comparison, it extracts the specified databases from both lists and
#' converts the results into data frames. The output is a nested list where each top-level element corresponds to a
#' comparison and contains a named list of data frames, with names indicating both the database and the regulation direction
#' (e.g., \code{"MSigDB_Hallmarks_UP"}).
#'
#' @param up_list A named list of GO analysis results for up-regulated genes. Each element should itself be a list
#'   with one entry per database. Each database entry must contain a \code{$data} element.
#' @param down_list A named list of GO analysis results for down-regulated genes. Each element should have the same
#'   structure as \code{up_list}.
#' @param databases A character vector of database names to extract from each result. Default is:
#'   \code{c("MSigDB_Hallmarks", "GO_Biological_Process_2023", "BioPlanet_2019",
#'   "GO_Cellular_Component_2023", "Reactome_Pathways_2024")} for continuity with up_go() and down_go() functions.
#'
#' @return A named list of data frames. Each element corresponds to a comparison and is itself a list
#'   of data frames. Each data frame is named to indicate the source database and whether it corresponds
#'   to up- or down-regulated genes.
#'
#' @details
#' For each comparison (the names that are common between \code{up_list} and \code{down_list}),
#' the function performs the following steps:
#' \itemize{
#'   \item Iterates over the specified \code{databases} for the up-regulated results and converts the result (stored in \code{$data})
#'         into a data frame, naming it with a suffix \code{"_UP"}.
#'   \item Iterates over the specified \code{databases} for the down-regulated results and converts these similarly,
#'         naming them with a suffix \code{"_Down"}.
#'   \item Combines these into a single list for that comparison.
#' }
#' This nested list can be further processed or exported (for example, with \code{writexl::write_xlsx}, where each top-level
#' element becomes a separate sheet).
#'
#' @examples
#' \dontrun{
#'   # Assume you have two named lists, one for up-regulated and one for down-regulated GO results:
#'   # Each result object (e.g., up_PB_vs_all) contains entries like up_PB_vs_all[["MSigDB_Hallmarks"]]
#'   # with a 'data' field.
#'   up_list <- list(
#'     PB_vs_all = list(
#'       MSigDB_Hallmarks = up_PB_vs_all$MSigDB_Hallmarks,
#'       GO_Biological_Process_2023 = up_PB_vs_all$GO_Biological_Process_2023,
#'       BioPlanet_2019 = up_PB_vs_all$BioPlanet_2019,
#'       GO_Cellular_Component_2023 = up_PB_vs_all$GO_Cellular_Component_2023,
#'       Reactome_Pathways_2024 = up_PB_vs_all$Reactome_Pathways_2024
#'     ),
#'     PB_vs_KLS = list(
#'       MSigDB_Hallmarks = up_PB_vs_OT$MSigDB_Hallmarks,
#'       GO_Biological_Process_2023 = up_PB_vs_OT$GO_Biological_Process_2023,
#'       BioPlanet_2019 = up_PB_vs_OT$BioPlanet_2019,
#'       GO_Cellular_Component_2023 = up_PB_vs_OT$GO_Cellular_Component_2023,
#'       Reactome_Pathways_2024 = up_PB_vs_OT$Reactome_Pathways_2024
#'     )
#'   )
#'
#'   down_list <- list(
#'     PB_vs_all = list(
#'       MSigDB_Hallmarks = down_PB_vs_all$MSigDB_Hallmarks,
#'       GO_Biological_Process_2023 = down_PB_vs_all$GO_Biological_Process_2023,
#'       BioPlanet_2019 = down_PB_vs_all$BioPlanet_2019,
#'       GO_Cellular_Component_2023 = down_PB_vs_all$GO_Cellular_Component_2023,
#'       Reactome_Pathways_2024 = down_PB_vs_all$Reactome_Pathways_2024
#'     ),
#'     PB_vs_KLS = list(
#'       MSigDB_Hallmarks = down_PB_vs_OT$MSigDB_Hallmarks,
#'       GO_Biological_Process_2023 = down_PB_vs_OT$GO_Biological_Process_2023,
#'       BioPlanet_2019 = down_PB_vs_OT$BioPlanet_2019,
#'       GO_Cellular_Component_2023 = down_PB_vs_OT$GO_Cellular_Component_2023,
#'       Reactome_Pathways_2024 = down_PB_vs_OT$Reactome_Pathways_2024
#'     )
#'   )
#'
#'   # Build the nested GO list:
#'   go_list <- build_go_list(up_list, down_list)
#' }
#'
#' @export
build_go_list <- function(up_list, down_list,
                          databases = c("MSigDB_Hallmarks", "GO_Biological_Process_2023",
                                        "BioPlanet_2019", "GO_Cellular_Component_2023", "Reactome_Pathways_2024")) {
  #identify comparisons common to both lists
  comparisons <- intersect(names(up_list), names(down_list))
  go_list <- setNames(vector("list", length(comparisons)), comparisons)

  for (comp in comparisons) {
    up_comp <- up_list[[comp]]
    down_comp <- down_list[[comp]]

    comp_list <- list()

    #loop over each database for up-regulated results
    for (db in databases) {
      # Extract the result, convert to data frame, and name with "_UP"
      comp_list[[paste0(db, "_UP")]] <- as.data.frame(up_comp[[db]]$data)
    }
    #loop over each database for down-regulated results
    for (db in databases) {
      #extract the result, convert to data frame, and name with "_Down"
      comp_list[[paste0(db, "_Down")]] <- as.data.frame(down_comp[[db]]$data)
    }
    go_list[[comp]] <- comp_list
  }

  return(go_list)
}

#' Combine and Prepare GO Analysis Results for Export
#'
#' This function takes a nested list of GO analysis results and, for each category, combines the contained data frames
#' into a single data frame by row-binding them. A new column \code{Database} is added to indicate the source of each record.
#' The resulting list is suitable for export to an Excel workbook (using packages like \code{writexl}), where each element becomes a separate sheet.
#'
#' @param go_list A named list of GO analysis results. Each element of \code{go_list} should itself be a list of data frames,
#'   each corresponding to a different database or comparison.
#'
#' @return A named list of data frames. Each data frame contains an additional column \code{Database} that indicates
#'   the source of the data.
#'
#' @details
#' For each category in \code{go_list}, the function assigns standardized names to the nested data frames,
#' then uses \code{dplyr::bind_rows} (within an \code{lapply} loop) to combine them into a single data frame.
#' This approach avoids explicit loops and improves efficiency over more verbose iterative methods.
#'
#' @examples
#' \dontrun{
#'   # Assume go_list is a nested list of GO results from various comparisons:
#'   combined_go_list <- prep_go_exp(go_list)
#'
#'   # Export to an Excel file where each list element becomes a separate sheet:
#'   writexl::write_xlsx(combined_go_list, "GO_Analysis_Results.xlsx")
#' }
#'
#' @import dplyr
#' @export
prep_go_exp <- function(go_list) {
  #standardized names for the nested list elements
  new_names <- c("MSigDB_Hallmarks_UP", "GO_Biological_Process_2023_UP", "BioPlanet_2019_UP",
                 "GO_Cellular_Component_2023_UP", "Reactome_Pathways_2024_UP",
                 "MSigDB_Hallmarks_Down", "GO_Biological_Process_2023_Down", "BioPlanet_2019_Down",
                 "GO_Cellular_Component_2023_Down", "Reactome_Pathways_2024_DOWN")

  #use lapply to process each category in the GO list efficiently
  combined_go_list <- lapply(go_list, function(sublist) {
    #assign standardized names to the sublist
    names(sublist) <- new_names[seq_along(sublist)]
    #combine the data frames with an added 'Database' column
    dplyr::bind_rows(sublist, .id = "Database")
  })

  return(combined_go_list)
}
