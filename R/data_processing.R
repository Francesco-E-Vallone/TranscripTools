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
  if (!is.data.frame(raw_counts)) raw_counts <- as.data.frame(raw_counts)
  
  for (file in files) {
    message("processing: ", file)
    data <- utils::read.delim(file, header = FALSE)
    
    file_name <- basename(file)
    colnames(data) <- c("gene", paste0("counts_", file_name))
    data <- data[!is.na(data$gene), , drop = FALSE]
    rownames(data) <- data$gene
    data$gene <- NULL
    
    #merge by rownames (gene IDs)
    raw_counts <- merge(
      raw_counts,
      data,
      by = "row.names",
      all = TRUE,
      sort = FALSE
    )
    rownames(raw_counts) <- raw_counts$Row.names
    raw_counts$Row.names <- NULL
  }
  
  #keep only counts_*
  raw_counts <- raw_counts[, grep("^counts_", colnames(raw_counts)), drop = FALSE]
  
  #common featureCounts extra rows (optional)
  raw_counts <- raw_counts[!(rownames(raw_counts) %in% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")), , drop = FALSE]
  
  raw_counts
}

#' Select Statistically Significant Differentially Expressed Genes
#'
#' Filters differential expression (DE) results by statistical significance and
#' absolute log fold-change. The input can be either:
#' \itemize{
#'   \item a single data frame of DE results, or
#'   \item a named (or unnamed) list of DE result data frames (e.g. one per contrast).
#' }
#'
#' @param dge A data frame of DE results or a list of data frames.
#'   Each data frame must contain:
#'   \itemize{
#'     \item a significance column (specified by \code{stat_sign})
#'     \item a \code{logFC} column (log fold-change)
#'   }
#' @param stat_sign Character string. Name of the significance column
#'   (e.g., \code{"p.value"}, \code{"adj.P.Val"}, \code{"padj"}).
#' @param p_threshold Numeric. Significance cutoff applied to \code{stat_sign}.
#'   Default is \code{0.05}.
#' @param logFC_threshold Numeric. Absolute log fold-change cutoff.
#'   Genes are kept when \code{abs(logFC) >= logFC_threshold}. Default is \code{1}.
#' @param keep_empty Logical. If \code{TRUE}, keeps empty filtered data frames in the output
#'   (useful when you want to preserve the full set of contrasts). If \code{FALSE} (default),
#'   drops empty results.
#'
#' @return If \code{dge} is a data frame, returns a filtered data frame.
#' If \code{dge} is a list, returns a list of filtered data frames.
#' Elements that are not data frames or that lack required columns are skipped with a warning.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Removes rows with \code{NA} in either \code{stat_sign} or \code{logFC}.
#'   \item Filters rows with \code{dge[[stat_sign]] < p_threshold}.
#'   \item Filters rows with \code{abs(logFC) >= logFC_threshold}.
#' }
#'
#' @examples
#' \dontrun{
#' dge_list <- list(
#'   condition1 = data.frame(p.value = runif(100), logFC = rnorm(100)),
#'   condition2 = data.frame(p.value = runif(100), logFC = rnorm(100))
#' )
#'
#' filtered_list <- select_stat_sign(dge_list, stat_sign = "p.value", logFC_threshold = 1)
#'
#' # Single data.frame also works:
#' filtered_df <- select_stat_sign(dge_list$condition1, stat_sign = "p.value", logFC_threshold = 1)
#' }
#'
#' @export
select_stat_sign <- function(
    dge,
    stat_sign,
    p_threshold = 0.05,
    logFC_threshold = 1,
    keep_empty = FALSE
) {
  
  # ---- internal worker: filter ONE data.frame ----
  filter_one <- function(df, label = NULL) {
    
    label_safe <- if (!is.null(label) && !is.na(label) && nzchar(label)) label else "<unnamed>"
    
    if (!is.data.frame(df)) {
      warning(sprintf("Skipping '%s': element is not a data.frame", label_safe))
      return(NULL)
    }
    
    # If empty input, return it as-is; keep/drop decision handled by caller
    if (nrow(df) == 0) {
      return(df)
    }
    
    # Required columns
    if (!(stat_sign %in% colnames(df))) {
      warning(sprintf("Skipping '%s': column '%s' not found", label_safe, stat_sign))
      return(NULL)
    }
    if (!("logFC" %in% colnames(df))) {
      warning(sprintf("Skipping '%s': column 'logFC' not found", label_safe))
      return(NULL)
    }
    
    # Drop NA rows in required fields (keeps filtering deterministic)
    ok <- !is.na(df[[stat_sign]]) & !is.na(df[["logFC"]])
    df <- df[ok, , drop = FALSE]
    
    # Apply thresholds
    df_filtered <- df[df[[stat_sign]] < p_threshold & abs(df[["logFC"]]) >= logFC_threshold, , drop = FALSE]
    
    df_filtered
  }
  
  # ---- Case 1: list input (note: data.frame is also a list, so check !is.data.frame) ----
  if (is.list(dge) && !is.data.frame(dge)) {
    
    nms <- names(dge)
    if (is.null(nms)) {
      # keep structure even if unnamed
      nms <- rep("", length(dge))
    }
    
    out <- vector("list", length(dge))
    
    for (i in seq_along(dge)) {
      out[[i]] <- filter_one(dge[[i]], label = nms[[i]])
    }
    
    # restore names only if there were any meaningful names
    if (!all(nms == "")) names(out) <- nms
    
    # drop NULLs (invalid inputs / missing columns)
    out <- out[!vapply(out, is.null, logical(1))]
    
    # optionally drop empty filtered results
    if (!keep_empty) {
      out <- out[vapply(out, nrow, integer(1)) > 0]
    }
    
    return(out)
  }
  
  # ---- Case 2: single data.frame input ----
  if (is.data.frame(dge)) {
    res <- filter_one(dge, label = "dge")
    if (is.null(res)) {
      stop("Input data.frame is missing required columns; see warnings.")
    }
    if (!keep_empty && nrow(res) == 0) {
      # For single df, returning empty is fine; keep_empty mainly matters for list
      return(res)
    }
    return(res)
  }
  
  stop("`dge` must be a data.frame or a list of data.frames.")
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
