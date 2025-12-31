# Create a Consolidated Data Frame from Multiple Raw Count Files

This function reads and merges raw count data from several separate
files into one data frame. It processes each file by reading the data,
assigning descriptive column names based on the file name, and appending
the count columns to the existing data. After processing, only columns
with count values are retained, and rows with labels such as
"no_feature" or "ambiguous" (if present) are removed.

## Usage

``` r
create_df(raw_counts, files)
```

## Arguments

- raw_counts:

  A data frame containing an initial set of raw counts. This can be an
  empty data frame or one that is pre-populated with some counts. It
  serves as the base to which new count data will be appended.

- files:

  A character vector of file paths to the raw count files. Each file is
  expected to be in a format where it has no header and contains two
  columns.

## Value

A data frame that consolidates the count columns from all the provided
files, with non-count columns and extra rows (such as "no_feature" or
"ambiguous") removed.

## Details

The function performs the following steps:

- Loops over each file in `files` and reads the data using
  [`read.delim`](https://rdrr.io/r/utils/read.table.html) with
  `header = FALSE`.

- Assigns column names to the data based on the file's base name,
  creating names like `genes_<filename>` and `counts_<filename>`.

- Sets the row names of the data frame using the first column, so that
  the count data can be merged by matching gene names.

- Appends the new data to `raw_counts` using
  [`bind_cols`](https://dplyr.tidyverse.org/reference/bind_cols.html).

- Retains only those columns whose names start with "counts",
  effectively discarding the gene identifier columns.

- Removes any rows where the row name matches `"no_feature"` or
  `"ambiguous"`, if present.

## Note

Before running this function, ensure you have defined the paths to your
raw count files and that each file adheres to the expected format (i.e.,
two columns with no header). Adjust the filtering of extra rows as
necessary based on your data.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create an initial raw counts data frame (this can be empty or contain initial data)
raw_counts <- data.frame(row.names = c("Gene1", "Gene2", "Gene3"))

# Define a list of file paths to the raw count files (modify the path and pattern as needed)
files <- list.files(path = "path/to/counts", pattern = "*.txt", full.names = TRUE)

# Consolidate the raw counts from all files into one data frame
consolidated_counts <- create_df(raw_counts, files)
} # }
```
