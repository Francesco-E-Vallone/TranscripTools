# Plot Differentially Expressed Genes from a Count Matrix

This function creates a bar plot of expression levels for the top
differentially expressed genes (DEGs) from a normalized count matrix
(e.g., TMM or TPM). The function filters the DEG data frame based on a
specified adjusted p-value threshold and log fold-change, then selects
and plots the top `top_n` genes. Gene names are italicized.

## Usage

``` r
deg_barplot(
  deg,
  tmm,
  title,
  pval_col = "adj.P.Val",
  logFC_col = "logFC",
  top_n = 30,
  direction = "UP"
)
```

## Arguments

- deg:

  A data frame containing differential expression results, including at
  least a log fold-change column and an adjusted p-value column. Row
  names should correspond to gene symbols.

- tmm:

  A normalized count matrix (e.g., TMM or TPM) with genes as rows and
  samples as columns.

- title:

  A character string specifying the title for the bar plot.

- pval_col:

  A character string specifying the column name in `deg` that contains
  adjusted p-values. Default is `"adj.P.Val"`.

- logFC_col:

  A character string specifying the column name in `deg` that contains
  log fold-changes. Default is `"logFC"`.

- top_n:

  An integer specifying the number of top DEGs to include in the plot.
  Default is `30`.

- direction:

  A character string indicating whether to visualize up-regulated genes
  ("UP") or down-regulated genes ("DOWN"). Default is `"UP"`. For
  up-regulated genes, only genes with `logFC > 1` are selected; for
  down-regulated genes, only genes with `logFC < -1` are selected.

## Value

A ggplot2 object representing a bar plot of the selected genes with
their expression levels.

## Details

The function first filters the `deg` data frame to retain statistically
significant genes (adjusted p-value below 0.05). It then sorts the genes
by log fold-change based on the specified `direction` ("UP" for
up-regulated or "DOWN" for down-regulated) and selects the top `top_n`
genes. Next, the corresponding rows are extracted from the count matrix
(`tmm`), and the data is reshaped into a long format for plotting with
ggplot2. Finally, a horizontal bar plot is created with gene names
displayed in italics.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage with TMM-normalized data:
# Assume deg_df is a DEG results data frame and tmm_avg is a TMM-normalized count matrix.
plot_up <- deg_barplot(deg_df, tmm_avg, "Top 30 UP Genes", direction = "UP")

# For TPM-normalized data, the same function can be used:
plot_down <- deg_barplot(deg_df, tpm_avg, "Top 30 DOWN Genes", direction = "DOWN")
} # }
```
