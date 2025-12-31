# Perform GO Analysis on Up-regulated Genes

This function performs Gene Ontology (GO) analysis on up-regulated genes
from a differential expression data frame. By default, it uses a preset
list of GO and pathway databases. However, you can supply your own
database vector via the `databases` parameter. The function also creates
bar plots for the top 20 GO terms (based on p-value) for each database.

## Usage

``` r
up_go(
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
