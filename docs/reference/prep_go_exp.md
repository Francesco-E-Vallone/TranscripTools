# Combine and Prepare GO Analysis Results for Export

This function takes a nested list of GO analysis results and, for each
comparison/category, combines the contained data frames into a single
data frame by row-binding them. A new column `Database` is added (from
the names of the nested list elements) to indicate the source database
and direction (e.g. `"GO_Biological_Process_2025_UP"`).

## Usage

``` r
prep_go_exp(go_list)
```

## Arguments

- go_list:

  A named list of GO analysis results. Each element of `go_list` should
  itself be a named list of data frames (or tibbles). Names should
  indicate database and direction (e.g.
  `"Reactome_Pathways_2024_Down"`).

## Value

A named list of data frames. Each data frame contains an additional
column `Database` indicating the source of each row (taken from the
nested list names).

## Details

Unlike earlier versions, this function does **not** overwrite names with
a fixed template. This makes it robust to missing databases (very common
in GO analyses) and prevents silent mislabeling.

## Examples

``` r
if (FALSE) { # \dontrun{
combined_go_list <- prep_go_exp(go_list)
writexl::write_xlsx(combined_go_list, "GO_Analysis_Results.xlsx")
} # }
```
