# ==============================================================================
# Title: ATAC-seq Peak Annotation & Functional Enrichment Analysis
# Description: Detailed genomic annotation of Differentially Accessible Regions
#              (DARs) using ChIPseeker, followed by GO, KEGG, and Reactome
#              pathway enrichment analysis.
# Source: YuLab-SMU ChIPseeker (https://github.com/YuLab-SMU/ChIPseeker)
# ==============================================================================

# --- 0. Package Loading ---
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)

# Load hg38 transcript database (consistent with BSgenome.Hsapiens.UCSC.hg38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Note: This script assumes 'sig_peaks_up' and 'sig_peaks_down' (GRanges objects
# with UCSC seqlevelsStyle) have been generated from Script 01.

# --- 1. Genomic Annotation & TSS Distribution ---
# Annotate up-regulated peaks relative to known gene features (Promoter, Exon, Intron, etc.)
peak_anno_up <- annotatePeak(sig_peaks_up,
                             tssRegion = c(-3000, 3000), # Promoter region: -3kb to +3kb from TSS
                             TxDb = txdb,
                             annoDb = "org.Hs.eg.db")

# Visualize genomic feature distribution and distance to Transcription Start Site (TSS)
plotAnnoBar(peak_anno_up)
plotAnnoPie(peak_anno_up)
plotDistToTSS(peak_anno_up, title = "Distribution of ATAC-seq Peaks Relative to TSS")

# Optional QC: Average accessibility profile around TSS (requires tagMatrix from peak/BAM files)
# tagMatrix <- getTagMatrix(sig_peaks_up, windows = getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000))
# plotAvgProf(tagMatrix, xlim = c(-3000, 3000), xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency")

# --- 2. Functional Enrichment Analysis (GO, KEGG, Reactome) ---
# Extract unique Entrez Gene IDs from annotated up-regulated peaks
genes_up <- unique(na.omit(as.data.frame(peak_anno_up)$geneId))

# 2.1 Gene Ontology (Biological Process)
enrich_go <- enrichGO(gene          = genes_up,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

dotplot(enrich_go, showCategory = 15, title = "GO Biological Process Enrichment")

# 2.2 Reactome Pathway Enrichment
enrich_reactome <- enrichPathway(gene         = genes_up,
                                 organism     = "human",
                                 pvalueCutoff = 0.05,
                                 readable     = TRUE)

dotplot(enrich_reactome, showCategory = 15, title = "Reactome Pathway Enrichment")

# 2.3 Many-to-Many Peak-to-Gene Mapping (seq2gene) & KEGG Enrichment
# Associates genomic regions with host genes, proximal genes, and flanking genes within 3kb
gene_flex <- seq2gene(sig_peaks_up,
                      tssRegion     = c(-3000, 3000),
                      flankDistance = 3000,
                      TxDb          = txdb)

enrich_kegg_flex <- enrichKEGG(gene          = gene_flex,
                               organism      = "hsa",
                               pvalueCutoff  = 0.05,
                               pAdjustMethod = "BH")

dotplot(enrich_kegg_flex, showCategory = 15, title = "KEGG Pathway Enrichment (seq2gene)")

# --- 3. Comparative Enrichment & Overlap Analysis (Up vs. Down DARs) ---
# Annotate down-regulated peaks to enable multi-group functional comparison
peak_anno_down <- annotatePeak(sig_peaks_down,
                               tssRegion = c(-3000, 3000),
                               TxDb = txdb,
                               annoDb = "org.Hs.eg.db")

genes_down <- unique(na.omit(as.data.frame(peak_anno_down)$geneId))

# Create a named list of gene sets for comparative clustering
genes_list <- list(UpRegulated = genes_up, DownRegulated = genes_down)

# Visualize overlapping target genes between Up and Down accessible regions
vennplot(genes_list)

# Compare KEGG pathways across Up- and Down-regulated DARs
comp_kegg <- compareCluster(geneCluster   = genes_list,
                            fun           = "enrichKEGG",
                            organism      = "hsa",
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH")

dotplot(comp_kegg, showCategory = 10, title = "Comparative KEGG Enrichment: Up vs. Down DARs")

# --- 4. Peak Overlap Permutation Test (Optional Reference) ---
# Statistically tests whether DARs overlap with external ChIP-seq/BED peak sets more than expected by chance
# enrichPeakOverlap(queryPeak     = files[[5]],
#                   targetPeak    = unlist(files[1:4]),
#                   TxDb          = txdb,
#                   pAdjustMethod = "BH",
#                   nShuffle      = 1000,
#                   verbose       = FALSE)
