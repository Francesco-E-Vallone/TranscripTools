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
#' @param df A data frame of differential expression results that must contain at least the columns
#'   \code{logFC} and \code{P.value}. Row names should represent gene symbols.
#' @param samples A character string that describes the sample group, used in plot titles.
#' @param databases A character vector specifying the databases to use for GO analysis. Default is:
#'   \code{c("MSigDB_Hallmark_2020", "GO_Biological_Process_2023", "BioPlanet_2019",
#'   "GO_Cellular_Component_2023", "Reactome_Pathways_2024")}.
#'   Importantly, only 5 databases are allowed at a time.
#'
#'
#' @return A named list of ggplot2 objects, each representing a bar plot of the top 20 significant GO terms for a database.
#'
#' @details
#' The function filters the input data frame to retain only those genes with \code{logFC} >= 1 (up-regulated genes),
#' extracts the gene symbols, and then uses the \code{enrichR} package to perform GO analysis on the provided databases.
#' For each database, only terms with \code{P.value} < 0.05 are retained and selected columns (columns 1, 3, 4, and 9, corresponding to "Terms", "P.value", "Adj. P.value", and "Genes", respectively)
#' are used for plotting. The function returns a list of ggplot objects, one for each database.
#'
#' @examples
#' \dontrun{
#'   # Assume df is your differential expression data frame with row names as gene symbols
#'   up_plots <- up_go(df, samples = "SampleGroup1")
#' }
#'
#' @import enrichR
#' @import ggplot2
#' @export
up_go <- function(df, samples, databases = c("MSigDB_Hallmark_2020",
                                             "GO_Biological_Process_2023",
                                             "BioPlanet_2019",
                                             "GO_Cellular_Component_2023",
                                             "Reactome_Pathways_2024")) {
  if (!requireNamespace("enrichR", quietly = TRUE)) stop("Package 'enrichR' is required")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")

  #select up-regulated genes (logFC >= 1) and extract gene symbols from row names
  df <- df[df$logFC >= 1, ]
  genes <- rownames(df)

  #perform GO analysis using enrichR on the specified databases
  go <- enrichR::enrichr(genes, databases = databases)

  #for each database, filter significant results (P.value < 0.05) and select specific columns
  msig <- go[[databases[1]]]
  msig <- msig[msig$P.value < 0.05, c(1, 3, 4, 9)]

  biop <- go[[databases[2]]]
  biop <- biop[biop$P.value < 0.05, c(1, 3, 4, 9)]

  bioplanet <- go[[databases[3]]]
  bioplanet <- bioplanet[bioplanet$P.value < 0.05, c(1, 3, 4, 9)]

  cellcomp <- go[[databases[4]]]
  cellcomp <- cellcomp[cellcomp$P.value < 0.05, c(1, 3, 4, 9)]

  react <- go[[databases[5]]]
  react <- react[react$P.value < 0.05, c(1, 3, 4, 9)]

  #combine all GO results into a named list
  all_go <- list(MSigDB_Hallmarks = msig,
                 GO_Biological_Process_2023 = biop,
                 BioPlanet_2019 = bioplanet,
                 GO_Cellular_Component_2023 = cellcomp,
                 Reactome_Pathways_2024 = react)

  #create bar plots for the top 20 GO terms from each database
  plots <- list()
  for (name in names(all_go)) {
    plot <- ggplot2::ggplot(na.omit(all_go[[name]][1:20, ])) +
      ggplot2::aes(x = reorder(Term, -P.value),
                   y = -log10(P.value)) +
      ggplot2::geom_col(fill = "#8FE0BC", color = "black") +
      ggplot2::coord_flip() +
      ggplot2::theme_light() +
      ggplot2::labs(title = paste(name, "| Up-regulated genes:", samples),
                    x = "Terms",
                    y = "-log10(P-value)") +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
    print(plot)
    plots[[name]] <- plot
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
#' @param df A data frame of differential expression results that must contain the columns \code{logFC} and \code{P.value}.
#'   Row names should represent gene symbols.
#' @param samples A character string indicating the sample group, which is used in the plot titles.
#' @param databases A character vector specifying the databases to use for GO analysis. Default is:
#'   \code{c("MSigDB_Hallmark_2020", "GO_Biological_Process_2023", "BioPlanet_2019",
#'   "GO_Cellular_Component_2023", "Reactome_Pathways_2024")}.
#'   Importantly, only 5 databases are allowed at a time.
#'
#' @return A named list of ggplot2 objects, each representing a bar plot for the top 20 significant GO terms for a database.
#'
#' @details
#' The function filters the input data frame to include only genes with \code{logFC} <= -1 (down-regulated),
#' then uses the \code{enrichR} package to perform GO analysis. For each database, terms with \code{P.value} < 0.05 are
#' selected and plotted. The resulting plots display the top 20 GO terms based on p-value.
#'
#' @examples
#' \dontrun{
#'   # Assume df is your differential expression data frame with gene symbols as row names
#'   down_plots <- down_go(df, samples = "SampleGroup2")
#' }
#'
#' @import enrichR
#' @import ggplot2
#' @export
down_go <- function(df, samples, databases = c("MSigDB_Hallmark_2020",
                                               "GO_Biological_Process_2023",
                                               "BioPlanet_2019",
                                               "GO_Cellular_Component_2023",
                                               "Reactome_Pathways_2024")) {
  if (!requireNamespace("enrichR", quietly = TRUE)) stop("Package 'enrichR' is required")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")

  #select down-regulated genes (logFC <= -1) and extract gene symbols
  df <- df[df$logFC <= -1, ]
  genes <- rownames(df)

  #perform GO analysis using the specified databases
  go <- enrichR::enrichr(genes, databases = databases)

  #filter and subset results for each database
  msig <- go[[databases[1]]]
  msig <- msig[msig$P.value < 0.05, c(1, 3, 4, 9)]

  biop <- go[[databases[2]]]
  biop <- biop[biop$P.value < 0.05, c(1, 3, 4, 9)]

  bioplanet <- go[[databases[3]]]
  bioplanet <- bioplanet[bioplanet$P.value < 0.05, c(1, 3, 4, 9)]

  cellcomp <- go[[databases[4]]]
  cellcomp <- cellcomp[cellcomp$P.value < 0.05, c(1, 3, 4, 9)]

  react <- go[[databases[5]]]
  react <- react[react$P.value < 0.05, c(1, 3, 4, 9)]

  #combine all filtered results into a named list
  all_go <- list(MSigDB_Hallmarks = msig,
                 GO_Biological_Process_2023 = biop,
                 BioPlanet_2019 = bioplanet,
                 GO_Cellular_Component_2023 = cellcomp,
                 Reactome_Pathways_2024 = react)

  #create bar plots for the top 20 GO terms for each database
  plots <- list()
  for (name in names(all_go)) {
    plot <- ggplot2::ggplot(na.omit(all_go[[name]][1:20, ])) +
      ggplot2::aes(x = reorder(Term, -P.value),
                   y = -log10(P.value)) +
      ggplot2::geom_col(fill = "#706F9B", color = "black") +
      ggplot2::coord_flip() +
      ggplot2::theme_light() +
      ggplot2::labs(title = paste(name, "| Down-regulated genes:", samples),
                    x = "Terms",
                    y = "-log10(P-value)") +
      ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
    print(plot)
    plots[[name]] <- plot
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
