# Generate Enrichment Plots from GO Data (from the enrichR package)

This function generates a list of enrichment plots for different
databases contained in the input list. Each element of the list `GO` is
processed to create a plot using the `plotEnrich` function.

## Usage

``` r
enrichPlot(GO)
```

## Arguments

- GO:

  A named list where each element contains data (e.g., a data frame)
  representing GO (Gene Ontology) enrichment results for a specific
  database.

## Value

A named list of plots. Each plot corresponds to one of the databases in
the input `GO` list.

## Details

The function loops through each database in the `GO` list, converts the
corresponding element to a data frame (if it is not already one), and
generates a plot using the `plotEnrich` function.

## Examples

``` r
# Example usage with a hypothetical GO list:
# GO <- list(
#   KEGG = kegg_data,
#   Reactome = reactome_data
# )
# plots <- enrichPlot(GO)
```
