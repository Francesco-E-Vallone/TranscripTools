# Build a Nested List of GO Analysis Results for Export

Combines up- and down-regulated GO/enrichment results into a nested list
suitable for export. Missing comparisons (e.g. only UP available) and
missing databases are allowed.

## Usage

``` r
build_go_list(
  up_list,
  down_list,
  databases,
  comparison_name = "comparison1",
  keep_empty = TRUE,
  warn = TRUE
)
```

## Arguments

- up_list:

  Up-regulated results (single-comparison or multi-comparison list).

- down_list:

  Down-regulated results (single-comparison or multi-comparison list).

- databases:

  Character vector of database names to extract.

- comparison_name:

  Name to use if a single-comparison structure is supplied.

- keep_empty:

  Logical. If TRUE (default), missing db/direction becomes empty
  data.frame so structure is stable for Excel export. If FALSE,
  missing/empty entries are omitted.

- warn:

  Logical. If TRUE (default), warn when results are missing/malformed.

## Value

A nested list where each top-level element is a `comparison` and
contains names like `paste0(db, "_UP")` and `paste0(db, "_Down")`.

## Details

Input can be either:

- Multi-comparison: `list(comparison -> list(database -> result))`

- Single-comparison: `list(database -> result)`

Each database result can be:

- a data.frame (already extracted), or

- a list-like object containing `$data` (your original assumption).
