# Draw a Volcano Plot for Differential Expression Results

This function creates a volcano plot from a differential expression
results data frame. It uses `ggplot2` and `ggrepel` to plot log2 fold
changes versus -log10(p-value), coloring points by their significance
category. It also labels the top `top_n` up‑ and down‑regulated genes.

## Usage

``` r
volplot(
  res,
  p_col = "PValue",
  logFC_col = "logFC",
  fdr_col = "FDR",
  fdr_threshold = 0.05,
  fc_threshold = 1,
  colors = c(Up = "#FF8F8F", Down = "#30B3A9", `Not Significant` = "grey"),
  top_n = 10,
  gene_col = NULL,
  title = "Volcano plot",
  subtitle = NULL,
  caption = NULL,
  x_lab = "Log2 fold change",
  y_lab = "-log10(p-value)",
  legend_title = NULL
)
```

## Arguments

- res:

  data.frame with at least `logFC_col` and `p_col`.

- p_col:

  Character. Column name for p-values used on y-axis. Default "PValue".

- logFC_col:

  Character. Column name for log2 fold change. Default "logFC".

- fdr_col:

  Character. Adjusted p-value column (used to derive `sign` if missing).
  Default "FDR".

- fdr_threshold:

  Numeric. Threshold for horizontal line and significance logic. Default
  0.05.

- fc_threshold:

  Numeric. Threshold for vertical lines and significance logic. Default
  1.

- colors:

  Named vector mapping c("Up","Down","Not Significant") to colors.

- top_n:

  Integer. Label top N up and top N down genes among significant rows.
  Default 10.

- gene_col:

  Optional character: column with gene names; if NULL uses rownames.

- title, subtitle, caption, x_lab, y_lab, legend_title:

  Character plot labels.

## Value

A ggplot object (returned invisibly).

## Details

Volcano plot (custom p-value column + fully customizable labels)

Draws log2FC (x) vs -log10(p) (y). Lets you declare which column is the
p-value, and customize all plot labels. If `sign` is absent, it will be
computed using
[`add_significance()`](https://francesco-e-vallone.github.io/TranscripTools/reference/add_significance.md)
with `fdr_col` when available (else it falls back to `p_col`).

## Examples

``` r
# p <- volplot(df, p_col="pval", logFC_col="log2FC", fdr_col="padj",
#              title="My volcano", x_lab="log2FC", y_lab="-log10(p)")
```
