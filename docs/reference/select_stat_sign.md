# Select Statistically Significant Differentially Expressed Genes

Filters differential expression (DE) results by statistical significance
and absolute log fold-change. The input can be either:

- a single data frame of DE results, or

- a named (or unnamed) list of DE result data frames (e.g. one per
  contrast).

## Usage

``` r
select_stat_sign(
  dge,
  stat_sign,
  p_threshold = 0.05,
  logFC_threshold = 1,
  keep_empty = FALSE
)
```

## Arguments

- dge:

  A data frame of DE results or a list of data frames. Each data frame
  must contain:

  - a significance column (specified by `stat_sign`)

  - a `logFC` column (log fold-change)

- stat_sign:

  Character string. Name of the significance column (e.g., `"p.value"`,
  `"adj.P.Val"`, `"padj"`).

- p_threshold:

  Numeric. Significance cutoff applied to `stat_sign`. Default is
  `0.05`.

- logFC_threshold:

  Numeric. Absolute log fold-change cutoff. Genes are kept when
  `abs(logFC) >= logFC_threshold`. Default is `1`.

- keep_empty:

  Logical. If `TRUE`, keeps empty filtered data frames in the output
  (useful when you want to preserve the full set of contrasts). If
  `FALSE` (default), drops empty results.

## Value

If `dge` is a data frame, returns a filtered data frame. If `dge` is a
list, returns a list of filtered data frames. Elements that are not data
frames or that lack required columns are skipped with a warning.

## Details

The function:

- Removes rows with `NA` in either `stat_sign` or `logFC`.

- Filters rows with `dge[[stat_sign]] < p_threshold`.

- Filters rows with `abs(logFC) >= logFC_threshold`.

## Examples

``` r
if (FALSE) { # \dontrun{
dge_list <- list(
  condition1 = data.frame(p.value = runif(100), logFC = rnorm(100)),
  condition2 = data.frame(p.value = runif(100), logFC = rnorm(100))
)

filtered_list <- select_stat_sign(dge_list, stat_sign = "p.value", logFC_threshold = 1)

# Single data.frame also works:
filtered_df <- select_stat_sign(dge_list$condition1, stat_sign = "p.value", logFC_threshold = 1)
} # }
```
