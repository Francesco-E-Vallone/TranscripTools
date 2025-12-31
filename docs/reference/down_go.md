# Perform GO Analysis on Down-regulated Genes

This function performs Gene Ontology (GO) analysis on down-regulated
genes from a differential expression data frame. It works similarly to
[`up_go()`](https://francesco-e-vallone.github.io/TranscripTools/reference/up_go.md)
but filters for genes with `logFC` \<= -1. By default, a preset list of
databases is used, but you can supply your own via the `databases`
parameter. Bar plots are generated for the top 20 GO terms.

## Usage

``` r
down_go(
  df,
  samples,
  databases = c("MSigDB_Hallmark_2020", "GO_Biological_Process_2023", "BioPlanet_2019",
    "GO_Cellular_Component_2023", "Reactome_Pathways_2024")
)
```

## Arguments

- df:

  Data frame of DE results with rownames = symbols and a `logFC` column.

- samples:

  Character used in plot titles.

- databases:

  Character vector of Enrichr databases to use (any length).

## Value

Named list of ggplot objects (one per database).

## Details

Perform GO Analysis on Down-regulated Genes (flexible databases)
