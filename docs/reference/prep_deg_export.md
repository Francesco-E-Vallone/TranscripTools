# Prepare Differential Expression Results for Export

This function takes a list of DEG result data frames and adds the row
names as a new column named `Genes`. This is useful for exporting DEG
results to file formats (such as XLSX, CSV, or TXT) where preserving
gene identifiers is necessary. When used with functions like
[`writexl::write_xlsx`](https://docs.ropensci.org/writexl//reference/write_xlsx.html),
each element of the list can be saved as a separate sheet in an Excel
workbook.

## Usage

``` r
prep_deg_export(res_list)
```

## Arguments

- res_list:

  A named list of data frames containing DEG results. Each data frame
  should have row names that represent gene identifiers.

## Value

A named list of data frames, where each data frame includes a new column
`Genes` containing the original row names.

## Details

The function iterates over each element in `res_list`, converts the row
names into a new column named `Genes`, and returns the updated list.
This transformation is helpful when exporting results, ensuring that
gene identifiers are not lost.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assume res_list is a list of DEG data frames obtained from an analysis:
res_list <- list(
  condition1 = data.frame(logFC = rnorm(50), p.value = runif(50), row.names = paste0("Gene", 1:50)),
  condition2 = data.frame(logFC = rnorm(50), p.value = runif(50), row.names = paste0("Gene", 51:100))
)

# Prepare the results for export by adding the gene identifiers as a column
res_deg <- prep_deg_exp(res_list)

# You can now export res_deg using writexl::write_xlsx(res_deg, "DEG_results.xlsx")
} # }
```
