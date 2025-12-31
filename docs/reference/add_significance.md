# Add Significance Categories to Differential Expression Results

This function adds a significance category to a differential expression
results data frame. The significance is determined using an FDR (or
adjusted p-value) column and a log2 fold change column. By default,
genes with FDR \< 0.05 and absolute logFC \> 1 are categorized as "Up"
(if logFC \> 1) or "Down" (if logFC \< -1).

## Usage

``` r
add_significance(
  res,
  fdr_col = "FDR",
  logFC_col = "logFC",
  fdr_threshold = 0.05,
  fc_threshold = 1
)
```

## Arguments

- res:

  data.frame of DE results (rows = genes).

- fdr_col:

  Character. Column name with adjusted p-values (e.g., FDR, padj).
  Default "FDR".

- logFC_col:

  Character. Column name with log2 fold change. Default "logFC".

- fdr_threshold:

  Numeric. Adjusted p-value cutoff. Default 0.05.

- fc_threshold:

  Numeric. Absolute log2FC cutoff. Default 1.

## Value

`res` with an added factor column `sign` in c("Up","Down","Not
Significant").

## Details

Flag differential expression direction by thresholds

Classifies each row as "Up", "Down", or "Not Significant" using a chosen
adjusted p-value column and a log2 fold-change column.

Up: logFC \> fc_threshold & FDR \< fdr_threshold Down: logFC \<
-fc_threshold & FDR \< fdr_threshold

## Examples

``` r
# df <- add_significance(df, fdr_col="adj.P.Val", logFC_col="logFC", fdr_threshold=0.05, fc_threshold=1)
```
