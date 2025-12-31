# Plot Gene Expression Statistics with Custom Labels and Statistical Comparisons

This function creates a faceted boxplot (with overlayed jittered points)
to display gene expression levels for a set of genes across samples. It
is designed to work with normalized count matrices such as TMM or TPM
and allows you to add pairwise statistical comparisons (e.g., via
"t.test", "wilcox.test", "anova"). Additionally, users can customize
plot title and axis labels.

## Usage

``` r
plot_genes_stats(
  tmm,
  genes,
  metadata,
  sample_col = "samples",
  group_col = "samples",
  test = "t.test",
  comparisons = NULL,
  log2p1 = TRUE,
  palette = NULL,
  facet_scales = "free_y",
  plot_title = "Gene Expression",
  x_lab = group_col,
  y_lab = if (log2p1) "log2(Expression + 1)" else "Expression",
  subtitle = NULL,
  caption = NULL
)
```

## Arguments

- tmm:

  Numeric matrix/data.frame (rows = genes, cols = samples).

- genes:

  Character vector of genes to plot.

- metadata:

  Data frame with sample metadata.

- sample_col:

  Character. Column in `metadata` containing sample IDs matching
  `colnames(tmm)`. Default "samples".

- group_col:

  Character. Column in `metadata` defining groups. Default "samples".

- test:

  Character. Statistical test for pairwise comparisons (e.g., "t.test",
  "wilcox.test", "anova"). Default "t.test".

- comparisons:

  Optional list of length-2 character vectors with group names for
  pairwise tests.

- log2p1:

  Logical. If TRUE, plot log2(Expression + 1). Default TRUE.

- palette:

  Optional named vector for group fill colors. Names must match group
  levels.

- facet_scales:

  Character, passed to `facet_wrap(scales=)`. Default "free_y".

- plot_title, x_lab, y_lab, subtitle, caption:

  Character plot labels. `x_lab` defaults to `group_col`; `y_lab`
  auto-set if `log2p1=TRUE`.

## Value

A ggplot object (invisible).
