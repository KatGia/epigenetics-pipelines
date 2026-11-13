# ==============================================================================
# Title: Single-Cell Multi-Omic Integration (scRNA-seq + scATAC-seq)
# Description: Preprocessing of scATAC-seq data (TF-IDF, LSI/SVD), calculation of
#              Gene Activity Scores, and Canonical Correlation Analysis (CCA)
#              anchor-based label transfer from a reference scRNA-seq dataset.
# Sources:
#   - Satija Lab Seurat v5 ATAC-seq Integration Vignette
#     (https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette)
#

# --- 0. Package Loading ---
library(Seurat)
library(SeuratData)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)

# --- 1. Load Reference scRNA-seq & Query scATAC-seq Data ---
# Load pre-processed and normalized scRNA-seq reference object (contains 'pbmc_normalized')
load("data/Blood_Seurat.RData")
# https://www.arraygen.com/scATASeq.php#:~:text=scATAC%2Dseq%20(Single%2Dcell,gene%20regulation%20and%20cellular%20heterogeneity.&text=Why%20Use%20scATAC%2Dseq
# https://satijalab.org/seurat/articles/seurat5_atacseq_integration_vignette

# Load 10x Genomics PBMC Multiome scATAC-seq dataset
pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")

# Remove low-quality / filtered cells
pbmc.atac <- subset(pbmc.atac, seurat_annotations != "filtered")

# --- 2. scATAC-seq Preprocessing & Dimensionality Reduction ---
# Extract hg38 gene annotations (required to compute Gene Activity Scores)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(pbmc.atac) <- annotations

# Normalization (TF-IDF) and feature selection
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")

# Latent Semantic Indexing (LSI) via Singular Value Decomposition (SVD)
pbmc.atac <- RunSVD(pbmc.atac)

# UMAP projection using LSI components 2:30 (component 1 typically correlates with sequencing depth)
pbmc.atac <- RunUMAP(pbmc.atac,
                     reduction = "lsi",
                     dims = 2:30,
                     reduction.name = "umap.atac",
                     reduction.key = "atacUMAP_")

# --- 3. Quantify Gene Activity Scores ---
# Assign pre-normalized scRNA-seq object as reference
pbmc.rna <- pbmc_normalized

# Identify variable features from the scRNA-seq reference
genes.to.use <- VariableFeatures(pbmc.rna)

# Estimate gene transcriptional activity from ATAC-seq chromatin accessibility
gene.activities <- GeneActivity(pbmc.atac, features = genes.to.use)

# Store Gene Activity Scores as a new assay and normalize/scale
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))

# --- 4. Multi-Omic Anchor Identification & Label Transfer ---
# Find transfer anchors between scRNA-seq (Reference) and scATAC-seq ACTIVITY assay (Query) via CCA
transfer.anchors <- FindTransferAnchors(
    reference = pbmc.rna,
    query = pbmc.atac,
    features = VariableFeatures(object = pbmc.rna),
    reference.assay = "RNA",
    query.assay = "ACTIVITY",
    reduction = "cca"
)

# Define reference cluster annotations to transfer
ref.labels <- pbmc.rna$seurat_clusters

# Transfer cell-type/cluster labels onto scATAC-seq cells weighted by LSI reduction
celltype.predictions <- TransferData(
    anchorset = transfer.anchors,
    refdata = ref.labels,
    weight.reduction = pbmc.atac[["lsi"]],
    dims = 2:30
)

# Add transferred metadata to the scATAC-seq Seurat object
pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)

# Visualize predicted cluster annotations on the ATAC-seq UMAP
DimPlot(pbmc.atac, group.by = "predicted.id", label = TRUE, reduction = "umap.atac") +
    labs(title = "Transferred Cell Annotations on scATAC-seq UMAP")

# Save integrated workspace (optional)
# save.image("results/scATAC_seq_integrated.RData")
