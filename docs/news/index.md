# Changelog

## TranscripTools 0.4.3

- Fix
  [`hmap()`](https://francesco-e-vallone.github.io/TranscripTools/reference/hmap.md)
  crash when ComplexHeatmap expects `gpar()` (row/column name gp now
  always valid).
- Improve top-annotation color mapping: missing levels are filled with
  `grey80` with a warning instead of erroring.
- Add optional control of row/column name font size and font face via
  `*_names_gp`.

## TranscripTools 0.4.2

- hmap(): fix annotation_colors/gpar error; add optional row/col name
  font controls; warn+fill missing annotation colors

## TranscripTools 0.4.1

- hmap(): added row/column name font size + font face options.

## TranscripTools 0.4.0

- pcafun: fixed shape handling (shapekey = metadata column; optional
  shape_map).
- prep_go_exp: robust GO export without hard-coded renaming; safe with
  missing databases.
- hmap: aligns metadata to sample order and checks ComplexHeatmap
  dependency.
- whiskyplot: expression is now log1p-transformed to match y-axis label.

## TranscripTools 0.3.0

- volplot: custom p_col + full label control
- hmap: explicit metadata + column
- plot_genes_stats / whiskyplot: custom labels
- deg_barplot: label controls
- up_go / down_go: flexible Enrichr databases
- new “average_cols” function
