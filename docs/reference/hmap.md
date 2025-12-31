# Create a Customized Heatmap for Count Data

This function generates a heatmap from a normalized count matrix (e.g.,
TMM or TPM data), with the option to restrict the visualization to a set
of differentially expressed genes (DEGs). In addition, it allows full
customization of the heatmap colors, the top annotation, and the display
of row (gene names) and column names.

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
  show_column_names = FALSE
)
```

## Arguments

- tmm:

  Numeric matrix or data.frame (rows = genes, cols = samples).

- title:

  Character. Heatmap title.

- deg:

  Optional character vector of genes to include.

- heatmap_colors:

  Color vector for the heatmap gradient.

- metadata:

  Data frame with sample metadata (rows = samples).

- metadata_col:

  Character. Column in `metadata` to use for top annotation.

- annotation_colors:

  Optional named vector for annotation levels.

- show_row_names, show_column_names:

  Logical. Show gene/sample names. Default FALSE.

## Value

A
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object (drawn and returned invisibly).
