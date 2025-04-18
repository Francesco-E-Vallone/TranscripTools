% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyses_GO.R
\name{save_plot}
\alias{save_plot}
\title{Save GO Analysis Plots to Files}
\usage{
save_plot(
  plot_list,
  list_name,
  path = "results/",
  device = "svg",
  height = 7,
  width = 15
)
}
\arguments{
\item{plot_list}{A named list of ggplot2 plot objects.}

\item{list_name}{A character string to prefix the file names.}

\item{path}{A character string specifying the directory path where the plots will be saved. Default is \code{"results/"}.}

\item{device}{A character string specifying the output device/format. Default is \code{"svg"}. Other possible values include \code{"png"}, \code{"pdf"}, etc.}

\item{height}{A numeric value specifying the height of the saved plot (in inches). Default is 7.}

\item{width}{A numeric value specifying the width of the saved plot (in inches). Default is 15.}
}
\value{
None. Files are written to disk.
}
\description{
This function saves a list of ggplot2 objects to files in a chosen format. The file names are constructed from the provided
\code{list_name} and the names of each plot in the list.
}
\details{
The function iterates over each plot in \code{plot_list} and uses \code{ggsave} to save each plot in the specified format.
}
\examples{
\dontrun{
  # Assume up_plots is a list of ggplot objects from up_go()
  save_plot(up_plots, list_name = "up_GO", device = "png")
}

}
