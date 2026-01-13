# Create a Customized Heatmap for Count Data

This function generates a heatmap from a normalized count matrix (e.g.,
TMM or TPM data), with the option to restrict the visualization to a set
of differentially expressed genes (DEGs). It supports a top annotation
built from sample metadata and lets you control whether row/column names
are shown, plus their font size and font face (via ComplexHeatmap's
`*_names_gp`).

## Usage

``` r
hmap(
  tmm,
  title,
  deg = NULL,
  heatmap_colors = c("#5a78b3", "#93b2d3", "#e7f4f2", "#faf0ab", "#db7d54", "#b9332c"),
  metadata = NULL,
  metadata_col = NULL,
  annotation_colors = NULL,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_names_fontsize = NULL,
  row_names_fontface = NULL,
  column_names_fontsize = NULL,
  column_names_fontface = NULL
)
```

## Arguments

- tmm:

  Numeric matrix or data.frame (rows = genes, cols = samples).

- title:

  Character. Heatmap title.

- deg:

  Optional character vector of genes to include (must match
  `rownames(tmm)`).

- heatmap_colors:

  Color vector for the heatmap gradient.

- metadata:

  Optional data.frame with sample metadata (rows = samples).

- metadata_col:

  Optional character. Column in `metadata` to use for top annotation.

- annotation_colors:

  Optional named vector mapping annotation levels to colors.

- show_row_names, show_column_names:

  Logical. Show gene/sample names. Default FALSE.

- row_names_fontsize:

  Numeric or NULL. Font size (points) for row names. NULL = default.

- row_names_fontface:

  Character/integer or NULL. Font face for row names (e.g., `"plain"`,
  `"bold"`, `"italic"`, `"bold.italic"`). NULL = default.

- column_names_fontsize:

  Numeric or NULL. Font size (points) for column names. NULL = default.

- column_names_fontface:

  Character/integer or NULL. Font face for column names (e.g.,
  `"plain"`, `"bold"`, `"italic"`, `"bold.italic"`). NULL = default.

## Value

A
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object (drawn and returned).
