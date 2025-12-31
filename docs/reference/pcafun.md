# Perform and Plot PCA with Customizable Aesthetics

This function performs Principal Component Analysis (PCA) on a count
matrix and generates a biplot using `PCAtools`. It streamlines the
workflow by computing the PCA and plotting in one step.

## Usage

``` r
pcafun(
  counts,
  metadata,
  title,
  colby = "samples",
  shapekey = NULL,
  shape_map = NULL,
  showLoadings = TRUE
)
```

## Arguments

- counts:

  A numeric matrix or data frame containing expression/count data (e.g.,
  TMM, CPM, TPM), with rows as features (e.g., genes) and columns as
  samples.

- metadata:

  A data frame containing sample metadata. Row names of `metadata` must
  correspond to the column names of `counts` (sample IDs). The function
  will reorder `metadata` to match `counts`.

- title:

  A character string for the PCA plot title.

- colby:

  A character string specifying the column name in `metadata` used to
  color samples in the PCA plot. Default is `"samples"`.

- shapekey:

  Optional character string specifying the column name in `metadata`
  used to shape points in the PCA plot. Default is `NULL`, meaning no
  grouping by shape.

- shape_map:

  Optional named vector mapping levels of `shapekey` to point shape
  codes (integers). For example: `c("batch1" = 15, "batch2" = 17)`.
  Default is `NULL` (PCAtools chooses shapes automatically).

- showLoadings:

  Logical indicating whether to display loadings on the PCA biplot.
  Default is `TRUE`.

## Value

A `ggplot` object representing the PCA biplot.

## Details

**Important note on shapes:** In
[`PCAtools::biplot()`](https://rdrr.io/pkg/PCAtools/man/biplot.html),
the argument `shape` specifies the metadata column used to assign point
shapes, while `shapekey` is an optional named vector mapping factor
levels to specific shape codes. For backward compatibility with older
usage patterns in this package, `pcafun()` uses `shapekey` to mean "the
metadata column used for shapes", and exposes the optional mapping via
`shape_map`.

The function computes PCA with
[`PCAtools::pca`](https://rdrr.io/pkg/PCAtools/man/pca.html) and then
plots using
[`PCAtools::biplot`](https://rdrr.io/pkg/PCAtools/man/biplot.html). It
performs basic input validation:

- Ensures `counts` has column names and `metadata` has row names.

- Reorders `metadata` rows to match `counts` columns.

- Checks that `colby` and (if provided) `shapekey` exist in `metadata`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assume counts is your count matrix and metadata contains corresponding sample information:
p1 <- pcafun(counts, metadata,
             title = "PCA",
             colby = "treatment",
             shapekey = "batch",
             showLoadings = FALSE)

# With explicit mapping of batch levels to shapes:
p2 <- pcafun(counts, metadata,
             title = "PCA",
             colby = "treatment",
             shapekey = "batch",
             shape_map = c("batch1" = 15, "batch2" = 17, "batch3" = 8),
             showLoadings = FALSE)
} # }
```
