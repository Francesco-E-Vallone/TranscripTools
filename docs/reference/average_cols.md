# Average Columns by Group (e.g., Replicates or Conditions)

Computes the mean expression for each gene across replicates that belong
to the same group, as defined in a metadata data frame.

## Usage

``` r
average_cols(tmm, meta)
```

## Arguments

- tmm:

  Numeric matrix or data frame of expression values (rows = genes,
  columns = samples).

- meta:

  Data frame with sample metadata. Must contain a column named `group`
  and row names matching the column names of `tmm`.

## Value

A data frame with the same row names as `tmm` and one column per group,
containing the average expression values for that group.

## Details

This function is useful when visualizing summarized or “cumulated”
profiles, such as average expression per condition for heatmaps or PCA.
Missing values are ignored (`na.rm = TRUE`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Average normalized counts by condition
tmm_avg <- average_cols(tmm, meta)

# Use averaged matrix for a heatmap
hmap(tmm_avg, title = "Group-Averaged Expression")
} # }
```
