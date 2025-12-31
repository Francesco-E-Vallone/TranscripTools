# Whiskyplot: Faceted Boxplots with Gene Expression and Pairwise Statistics

Subsets a normalized count matrix (TMM or TPM) to a set of genes, merges
with sample metadata, and draws faceted boxplots (“whisker plots”) with
jittered points. Pairwise comparisons are computed via **rstatix** and
annotated by **ggpubr**.

## Usage

``` r
whiskyplot(
  tmm,
  genes,
  metadata,
  sample_col = "samples",
  group_col = "samples",
  test = c("wilcox", "t"),
  p.adjust.method = "none",
  comparisons = NULL,
  palette = NULL,
  plot_title = "Gene expression",
  x_lab = group_col,
  y_lab = "log(TMM/TPM + 1)",
  subtitle = NULL,
  caption = NULL
)
```

## Arguments

- tmm:

  Numeric matrix/data.frame (rows = genes, cols = samples).

- genes:

  Character vector of gene names (must match rownames(tmm)).

- metadata:

  Data frame with sample metadata.

- sample_col:

  Character. Column in `metadata` containing sample IDs (matches
  colnames(tmm)). Default "samples".

- group_col:

  Character. Column in `metadata` defining groups. Default "samples".

- test:

  One of c("wilcox","t"). Default "wilcox".

- p.adjust.method:

  P-value adjustment for rstatix ("none","BH",...). Default "none".

- comparisons:

  Optional list of length-2 character vectors (group names). If NULL,
  all pairwise combos.

- palette:

  Optional named vector for group fill colors.

- plot_title, x_lab, y_lab, subtitle, caption:

  Character plot labels. `x_lab` defaults to `group_col`.

## Value

A ggplot object.
