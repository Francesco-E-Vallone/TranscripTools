# Filter TPM Data by Pathway and Differential Expression

This function filters a TPM data matrix to include only the genes that
belong to a specified pathway and are also identified as differentially
expressed. This is particularly useful for conducting pathway-specific
analyses.

## Usage

``` r
deg_in_pathway(tpm, deg, pathway)
```

## Arguments

- tpm:

  A matrix or data frame containing TPM (Transcripts Per Million)
  values, where row names correspond to gene names.

- deg:

  A data frame containing differential expression results. The row names
  of this data frame should correspond to the genes identified as
  differentially expressed.

- pathway:

  A vector of gene names representing the pathway of interest.

## Value

A subset of the TPM data containing only the rows (genes) that are
present in both the provided `pathway` vector and in the row names of
`deg`.

## Details

This function is designed for pathway-specific analyses. It intersects
the set of genes from the specified pathway with those identified as
differentially expressed (from `deg`). Ensure that the `deg` data frame
is formatted correctly with gene names as row names.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create an example TPM matrix with gene names as row names
tpm <- matrix(runif(100), nrow = 10,
              dimnames = list(paste0("Gene", 1:10), paste0("Sample", 1:10)))

# Create an example differential expression results data frame
deg <- data.frame(logFC = rnorm(10), p.value = runif(10))
rownames(deg) <- paste0("Gene", sample(1:10))

# Define a pathway of interest (vector of gene names)
pathway_genes <- c("Gene1", "Gene3", "Gene5")

# Filter the TPM matrix for genes in the pathway that are also differentially expressed
filtered_tpm <- deg_in_pathway(tpm, deg, pathway_genes)
} # }
```
