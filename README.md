# Epigenomics &amp; Chromatin Accessibility
# ATAC-seq Differential Accessibility Analysis Pipeline

An end-to-end R pipeline for identifying **Differentially Accessible Regions (DARs)** from ATAC-seq data, featuring robust linear modeling, quality controls, genomic feature annotation, and transcription factor motif enrichment.

## Pipeline Overview

1. **Preprocessing & Quality Control:** Count matrix filtering via `edgeR::filterByExpr` and TMM normalization.
2. **Differential Accessibility Testing:** Log-CPM transformation and precision weighting using **Limma-Voom**, followed by empirical Bayes moderation (`eBayes`).
3. **Visualization:** Volcano plots, MA plots, sample MDS/boxplots, and hierarchical clustering heatmaps of significant peaks.
4. **Genomic Annotation:** Integration with ENCODE Regulatory Build (`.gff`) using custom `GenomicRanges` overlap functions.
5. **Motif Enrichment:** Transcription factor binding site (TFBS) enrichment analysis using **monaLisa** and the **JASPAR2020** database against a genomic background.

## Dependencies

- R >= 4.2
- Bioconductor Packages: `GenomicRanges`, `GenomicAlignments`, `edgeR`, `limma`, `BSgenome.Hsapiens.UCSC.hg38`, `monaLisa`, `JASPAR2020`, `TFBSTools`
- CRAN Packages: `ggplot2`, `ggrepel`, `pheatmap`, `tidyr`, `stringr`, `magrittr`
