<h1 align="center">

TranscripTools

</h1>

<p align="center">

<img src="man/figures/Tlogo.png" alt="TranscripTools logo" width="400"/>

</p>

TranscripTools is an R package born out of practical necessity. It collects a set of lightweight helper functions for bulk RNA-seq analyses, developed while working on multiple projects with many differential expression comparisons per project.

The package focuses on common but repetitive tasks: handling and filtering differential expression results, generating quick and consistent visualizations (volcano plots, heatmaps, PCA), running basic enrichment analyses, and exporting results in analysis-ready formats for reporting. While it was developed for the author’s projects, the same workflow pattern shows up in many labs—especially when you need to inspect and report results across multiple comparisons.

This is not intended to be a general-purpose framework. It is a personal research toolbox, built to reduce repetitive code, keep analyses readable, and make it easier to iterate over many contrasts without rewriting the same logic every time. Improving reproducibility comes as a consequence of having fewer ad-hoc scripts and more consistent interfaces.

TranscripTools is under active development; while it was built around personal workflows, it is meant to be useful to anyone running bulk RNA-seq analyses, especially when working through many comparisons.

## Installation

```r
# Install devtools if needed
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Install TranscripTools from GitHub
devtools::install_github("Francesco-E-Vallone/TranscripTools")

# Load the package
library(TranscripTools)

# Optional (Bioconductor):
# - PCA plotting: PCAtools
# - Heatmaps: ComplexHeatmap
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("PCAtools", "ComplexHeatmap"))
```

## Documentation

Function reference and examples are available via the pkgdown site:

[TranscripTools documentation](https://francesco-e-vallone.github.io/TranscripTools/reference/index.html)

## Workflow example (minimal)

If you want to see how the functions in TranscripTools fit together, start from this **minimal workflow**.

It’s written as a practical “how I run it” guide (no toy dataset included): it walks through a typical bulk RNA-seq downstream loop with multiple contrasts; from DE results to PCA, DEG filtering, enrichment/pathway steps, focused gene exploration, visualization, and export.

[Practical analysis workflow (minimal)](https://francesco-e-vallone.github.io/TranscripTools/articles/TranscripTools_workflow.html)

## Collaboration

This package is openly developed. If you find it useful, have suggestions, or want to collaborate on extending or refining parts of it, feel free to get in touch! :)

## Note on AI assistance

Some documentation text and examples were drafted with the help of a large language model (LLM) and then reviewed and edited by the author. The code, design choices, and final responsibility for the package remain the author’s.
