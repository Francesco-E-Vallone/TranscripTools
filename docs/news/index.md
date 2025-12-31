# Changelog

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
