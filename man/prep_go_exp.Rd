% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyses_GO.R
\name{prep_go_exp}
\alias{prep_go_exp}
\title{Combine and Prepare GO Analysis Results for Export}
\usage{
prep_go_exp(go_list)
}
\arguments{
\item{go_list}{A named list of GO analysis results. Each element of \code{go_list} should itself be a list of data frames,
each corresponding to a different database or comparison.}
}
\value{
A named list of data frames. Each data frame contains an additional column \code{Database} that indicates
the source of the data.
}
\description{
This function takes a nested list of GO analysis results and, for each category, combines the contained data frames
into a single data frame by row-binding them. A new column \code{Database} is added to indicate the source of each record.
The resulting list is suitable for export to an Excel workbook (using packages like \code{writexl}), where each element becomes a separate sheet.
}
\details{
For each category in \code{go_list}, the function assigns standardized names to the nested data frames,
then uses \code{dplyr::bind_rows} (within an \code{lapply} loop) to combine them into a single data frame.
This approach avoids explicit loops and improves efficiency over more verbose iterative methods.
}
\examples{
\dontrun{
  # Assume go_list is a nested list of GO results from various comparisons:
  combined_go_list <- prep_go_exp(go_list)

  # Export to an Excel file where each list element becomes a separate sheet:
  writexl::write_xlsx(combined_go_list, "GO_Analysis_Results.xlsx")
}

}
