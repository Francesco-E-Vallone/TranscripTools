#' Create a Consolidated Data Frame from Multiple Raw Count Files
#'
#' This function reads and merges raw count data from several separate files into one data frame.
#' It processes each file by reading the data, assigning descriptive column names based on the file name,
#' and appending the count columns to the existing data. After processing, only columns with count values
#' are retained, and rows with labels such as "no_feature" or "ambiguous" (if present) are removed.
#'
#' @param raw_counts A data frame containing an initial set of raw counts. This can be an empty data frame
#'   or one that is pre-populated with some counts. It serves as the base to which new count data will be appended.
#' @param files A character vector of file paths to the raw count files. Each file is expected to be in a format
#'   where it has no header and contains two columns.
#'
#' @return A data frame that consolidates the count columns from all the provided files, with non-count columns and
#'   extra rows (such as "no_feature" or "ambiguous") removed.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Loops over each file in \code{files} and reads the data using \code{\link[utils]{read.delim}} with \code{header = FALSE}.
#'   \item Assigns column names to the data based on the file's base name, creating names like \code{genes_<filename>} and
#'   \code{counts_<filename>}.
#'   \item Sets the row names of the data frame using the first column, so that the count data can be merged by matching gene names.
#'   \item Appends the new data to \code{raw_counts} using \code{\link[dplyr]{bind_cols}}.
#'   \item Retains only those columns whose names start with "counts", effectively discarding the gene identifier columns.
#'   \item Removes any rows where the row name matches \code{"no_feature"} or \code{"ambiguous"}, if present.
#' }
#'
#' @note
#' Before running this function, ensure you have defined the paths to your raw count files and that each file adheres
#' to the expected format (i.e., two columns with no header). Adjust the filtering of extra rows as necessary based
#' on your data.
#'
#' @examples
#' \dontrun{
#' # Create an initial raw counts data frame (this can be empty or contain initial data)
#' raw_counts <- data.frame(row.names = c("Gene1", "Gene2", "Gene3"))
#'
#' # Define a list of file paths to the raw count files (modify the path and pattern as needed)
#' files <- list.files(path = "path/to/counts", pattern = "*.txt", full.names = TRUE)
#'
#' # Consolidate the raw counts from all files into one data frame
#' consolidated_counts <- create_df(raw_counts, files)
#' }
#'
#' @import dplyr
#' @importFrom utils head read.delim

#' @export
create_df <- function(raw_counts, files) {
  #looping for each file
  for (file in files) {
    print(paste("processing: ", file)) # Debug message
    data <- read.delim(file, header = FALSE)

    #setting proper names based on the file's base name
    file_name <- basename(file)
    colnames(data) <- c(paste0("genes_", file_name), paste0("counts_", file_name))

    #setting row names so that files can be appended by matching gene names
    rownames(data) <- data[,1]

    #appending the new data columns without overwriting existing ones
    raw_counts <- bind_cols(raw_counts, data)
  }

  #retaining only columns that start with "counts"
  to_retain <- grep("^counts", colnames(raw_counts), value = TRUE)
  raw_counts <- raw_counts[, colnames(raw_counts) %in% to_retain]

  #removing rows for genes labeled as "no_feature" or "ambiguous", if present
  to_remove <- c("no_feature", "ambiguous")
  raw_counts <- raw_counts[!(rownames(raw_counts) %in% to_remove), ]

  return(raw_counts)
}

#' Select Statistically Significant Differentially Expressed Genes
#'
#' This function processes a named list of differential gene expression (DGE) data frames and filters each
#' data frame to retain only those genes that are statistically significant based on a specified significance 
#' column and a log fold-change threshold. This is useful for selecting the most differentially expressed genes 
#' across multiple analyses.
#'
#' @param dge_list A named list of data frames, each containing DGE results. Each data frame should have a numeric column 
#'   for significance (e.g., p-value or adjusted p-value) and a \code{logFC} column indicating the log fold-change.
#' @param stat_sign A character string specifying the name of the significance column (e.g., `"p.value"` or `"adj.P.Val"`)
#'   used to determine statistical significance.
#' @param logFC_threshold A numeric value specifying the threshold for the absolute log fold-change.
#'   Only genes with \code{logFC} greater than or equal to this threshold (or less than or equal to the negative 
#'   of this threshold) will be retained. Default is \code{1}.
#'
#' @return A named list of data frames, where each data frame contains only the statistically significant genes from 
#'   the corresponding element in \code{dge_list}. Any element of \code{dge_list} that is empty or lacks the specified 
#'   columns is skipped.
#'
#' @details
#' The function loops over each element in \code{dge_list} and for each one:
#' \itemize{
#'   \item Checks if the data frame is empty; if so, it moves on to the next element.
#'   \item Verifies that both the significance column (given by \code{stat_sign}) and the \code{logFC} column are present.
#'   \item Filters the data frame to retain only those genes for which the significance value is less than 0.05 and 
#'         the absolute value of \code{logFC} is greater than or equal to \code{logFC_threshold}.
#' }
#'
#' @examples
#' \dontrun{
#' # Example: Suppose you have a list of DGE results:
#' dge_list <- list(
#'   condition1 = data.frame(p.value = runif(100), logFC = rnorm(100)),
#'   condition2 = data.frame(p.value = runif(100), logFC = rnorm(100))
#' )
#'
#' # Filter the list using "p.value" as the significance column and a logFC threshold of 1
#' filtered_list <- select_stat_sign(dge_list, stat_sign = "p.value", logFC_threshold = 1)
#' }
#'
#' @export
select_stat_sign <- function(dge_list, stat_sign, logFC_threshold = 1) {
  #create an empty list to store filtered results
  filt_res <- list()
  
  #loop over each element (by name) in dge_list
  for (res_name in names(dge_list)) {
    df <- dge_list[[res_name]]
    
    #skip if the data frame is empty
    if (nrow(df) == 0) {
      next
    }
    
    #check that the specified significance column and the "logFC" column exist
    if (!(stat_sign %in% colnames(df))) {
      warning(sprintf("Skipping '%s': Column '%s' not found", res_name, stat_sign))
      next
    }
    if (!("logFC" %in% colnames(df))) {
      warning(sprintf("Skipping '%s': Column 'logFC' not found", res_name))
      next
    }
    
    #filter the data frame based on the significance threshold and log fold-change threshold
    df_filtered <- df[df[[stat_sign]] < 0.05 & (df$logFC >= logFC_threshold | df$logFC <= -logFC_threshold), ]
    
    #sdd the filtered data frame to the result list under the same name
    filt_res[[res_name]] <- df_filtered
  }
  
  return(filt_res)
}



#' Filter TPM Data by Pathway and Differential Expression
#'
#' This function filters a TPM data matrix to include only the genes that belong
#' to a specified pathway and are also identified as differentially expressed.
#' This is particularly useful for conducting pathway-specific analyses.
#'
#' @param tpm A matrix or data frame containing TPM (Transcripts Per Million)
#'   values, where row names correspond to gene names.
#' @param deg A data frame containing differential expression results. The row
#'   names of this data frame should correspond to the genes identified as
#'   differentially expressed.
#' @param pathway A vector of gene names representing the pathway of interest.
#'
#' @return A subset of the TPM data containing only the rows (genes) that are
#'   present in both the provided \code{pathway} vector and in the row names of
#'   \code{deg}.
#'
#' @details
#' This function is designed for pathway-specific analyses. It intersects the set
#' of genes from the specified pathway with those identified as differentially
#' expressed (from \code{deg}). Ensure that the \code{deg} data frame is formatted
#' correctly with gene names as row names.
#'
#' @examples
#' \dontrun{
#' # Create an example TPM matrix with gene names as row names
#' tpm <- matrix(runif(100), nrow = 10,
#'               dimnames = list(paste0("Gene", 1:10), paste0("Sample", 1:10)))
#'
#' # Create an example differential expression results data frame
#' deg <- data.frame(logFC = rnorm(10), p.value = runif(10))
#' rownames(deg) <- paste0("Gene", sample(1:10))
#'
#' # Define a pathway of interest (vector of gene names)
#' pathway_genes <- c("Gene1", "Gene3", "Gene5")
#'
#' # Filter the TPM matrix for genes in the pathway that are also differentially expressed
#' filtered_tpm <- deg_in_pathway(tpm, deg, pathway_genes)
#' }
#'
#' @export
deg_in_pathway <- function(tpm, deg, pathway) {
  tpm <- tpm[rownames(tpm) %in% pathway, ]
  tpm <- tpm[rownames(tpm) %in% rownames(deg), ]
  return(tpm)
}


#' Prepare Differential Expression Results for Export
#'
#' This function takes a list of DEG result data frames and adds the row names as a new column named \code{Genes}.
#' This is useful for exporting DEG results to file formats (such as XLSX, CSV, or TXT) where preserving gene identifiers
#' is necessary. When used with functions like \code{writexl::write_xlsx}, each element of the list can be saved as a separate
#' sheet in an Excel workbook.
#'
#' @param res_list A named list of data frames containing DEG results. Each data frame should have row names that represent gene identifiers.
#'
#' @return A named list of data frames, where each data frame includes a new column \code{Genes} containing the original row names.
#'
#' @details
#' The function iterates over each element in \code{res_list}, converts the row names into a new column named \code{Genes},
#' and returns the updated list. This transformation is helpful when exporting results, ensuring that gene identifiers
#' are not lost.
#'
#' @examples
#' \dontrun{
#' # Assume res_list is a list of DEG data frames obtained from an analysis:
#' res_list <- list(
#'   condition1 = data.frame(logFC = rnorm(50), p.value = runif(50), row.names = paste0("Gene", 1:50)),
#'   condition2 = data.frame(logFC = rnorm(50), p.value = runif(50), row.names = paste0("Gene", 51:100))
#' )
#'
#' # Prepare the results for export by adding the gene identifiers as a column
#' res_deg <- prep_deg_exp(res_list)
#'
#' # You can now export res_deg using writexl::write_xlsx(res_deg, "DEG_results.xlsx")
#' }
#'
#' @export
prep_deg_export <- function(res_list) {
  #create an empty list to store modified data frames
  res_deg <- list()

  #loop over each element in res_list
  for (nested in names(res_list)) {
    #extract the current data frame
    df <- res_list[[nested]]
    #add row names as a new column named 'Genes'
    df <- data.frame(Genes = rownames(df), df, check.names = FALSE, stringsAsFactors = FALSE)
    #store the modified data frame in the result list
    res_deg[[nested]] <- df
  }

  return(res_deg)
}
