# ==============================================================================
# Title: ATAC-seq Differential Accessibility & Motif Enrichment Analysis
# Description: Identification of Differentially Accessible Regions (DARs), 
#              genomic feature annotation, and transcription factor motif enrichment.
# Sources: 
#   - Nature Communications / ENCODE regulatory workflows
#   - Sam Buckberry's ATAC-seq Workshop (https://github.com/SamBuckberry/ATACseq-workshop)
# ==============================================================================

# --- 0. Package Loading ---
library(GenomicRanges)
library(GenomicAlignments)
library(edgeR)

library(BSgenome.Hsapiens.UCSC.hg38) 
# Contains the complete human genome sequence (hg38 version), which is essential for extracting the DNA sequences required for motif enrichment analysis.

library(monaLisa)
library(JASPAR2020)
library(TFBSTools)
# The monaLisa, JASPAR2020, and TFBSTools packages are used for transcription factor (TF) motif enrichment analysis, in order to identify which TFs bind to differential peaks.

library(stringr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(tidyr)

# --- 1. Setup & Data Loading ---
# Import peak counts matrix: rows = genomic peaks/regions, columns = samples (Stimulated/Unstimulated)
peak_counts <- read.table("atac_peak_counts.tsv", header = TRUE, row.names = 1)
head(peak_counts)

# --- 2. Experimental Design & Normalization ---
# Define treatment and control groups based on sample column names
groups <- ifelse(test = grepl(pattern = "_T.U", x = colnames(peak_counts)),
                 yes = "Unstimulated", no = "Stimulated")
samples <- str_sub(string = colnames(peak_counts), start = 2, end = 5)
samples <- str_c(samples, groups, sep = "_")

# Create DGEList object
y <- DGEList(counts = peak_counts, group = groups, samples = samples)

# Design matrix for linear modeling (~ 0 + groups estimates mean accessibility per group)
design <- model.matrix(~ 0 + groups)
colnames(design) <- c("Stimulated", "Unstimulated")

# Filter out lowly accessible peaks with low counts across samples
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes = FALSE]

# Calculate TMM normalization factors for sequencing depth and compositional bias
y <- calcNormFactors(y)

# --- 3. Differential Analysis (Limma-Voom) ---
# Transform count data to log2-CPM and calculate precision weights
v <- voom(counts = y, design = design, plot = TRUE)

#  What does voom do?
# It converts ATAC counts into log2 Counts Per Million (log-CPM) and calculates specific weights for each peak. This adjusts count data (where variance depends on the mean) for use in linear models (which require homogeneous variance).
# A linear model is fitted to the log2 CPM for each peak, and the residuals are calculated.
# A smoothed curve is fitted to the square root of the residual standard deviation as a function of the mean counts for the peaks (see the red line in the plot above).
# This smoothed curve is used to derive weights for each peak and sample, which are then passed to limma along with the log2 CPM values.
# Further details on the application of this method to RNA-seq data can be found here:  https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29

# Quality checks
boxplot(v$E, names = v$targets$samples, main = "Log-CPM Boxplot per Sample")
plotMDS(v$E, main = "MDS Plot of ATAC-seq Samples")

# If the boxplots appear notably different across the samples (the ones above seem fine), quantile normalization might be useful.

# Fit linear model
fit <- lmFit(object = v, design = design)

# Contrast matrix: comparing Stimulated vs Unstimulated
contr <- makeContrasts(Stimulated - Unstimulated, levels = colnames(coef(fit)))
fit2 <- contrasts.fit(fit, contr)

# Empirical Bayes moderation of standard errors
eb <- eBayes(fit2)

# Empirical Bayesian smoothing of standard errors (see https://www.degruyter.com/doi/10.2202/1544-6115.1027)
# This improves variance (standard error) estimates for individual peaks by "shrinking" them toward a common mean. This is crucial for stabilizing estimates and increasing statistical power, especially when the number of replicates is limited.

# Extract full results table: check for the highest differential peaks
# Extract the complete test results (Log Fold Change, P-value, adjusted P-value) and define a peak as significant
tt2 <- topTable(eb, number = nrow(v))
tt2$Significant <- ((abs(tt2$logFC) > 1) & (tt2$adj.P.Val < 0.05))

# How many DE peaks are there?
cat("Number of significant differentially accessible peaks:\n")
print(table(tt2$Significant))

# --- 4. Results Visualization ---
# Volcano Plot
gg_volcano <- ggplot(tt2, aes(x = logFC, y = -log10(P.Value), fill = Significant, colour = Significant)) +
    geom_point(alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    theme_bw() +
    ggtitle("Volcano Plot: ATAC-seq DARs")
gg_volcano

# MA Plot
gg_ma <- ggplot(tt2, aes(x = AveExpr, y = logFC, colour = Significant)) +
    geom_point(alpha = 0.5) +
    geom_vline(xintercept = c(1), linetype = "dashed") +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
    theme_bw() +
    ggtitle("MA Plot: Accessibility vs Fold Change")
gg_ma

# * Volcano Plot: Displays LogFC (difference in accessibility) versus P-value (significance), with differential peaks highlighted in red.
# * MA Plot: Displays LogFC versus AveExpr (average accessibility). It should show a pattern centered around zero, indicating that the differential results are not influenced by average accessibility.

# Heatmap of significant DE peaks
# Visualize normalized accessibility (log-CPM) for significant peaks only in a heatmap, showing how peaks turn on or off consistently in stimulated versus non-stimulated samples.
top_de_cpm <- v$E[rownames(v$E) %in% rownames(tt2)[tt2$Significant == TRUE], ]
colnames(top_de_cpm) <- v$targets$samples
pheatmap(top_de_cpm, scale = "row", show_rownames = FALSE, main = "Significant DARs Heatmap")

# Helper function to plot normalized counts for top peaks
plot_de <- function(ids){
    dat <- v$E[rownames(v$E) %in% ids, ] %>% data.frame()
    colnames(dat) <- v$targets$samples
    dat <- tibble::rownames_to_column(dat, var = "id") %>%
        tidyr::gather(sample, value, -id)
    dat$group <- v$targets$group[match(dat$sample, v$targets$samples)]
    
    ggplot(data = dat, aes(x = group, y = value)) + 
        geom_point(size = 3, alpha = 0.5, position = position_jitter(width = 0.1)) +
        facet_wrap(~id, scales = "free_y") +
        ylab("Normalised log-CPM") +
        theme_bw()
}

my_top_peaks2 <- rownames(tt2)[tt2$Significant == TRUE & tt2$logFC < 0][1:6] # down
# my_top_peaks2 <- rownames(tt2)[tt2$Significant == TRUE & tt2$logFC > 0][1:6] # up
plot_de(ids = my_top_peaks2)

# --- 5. Functional Interpretation (Annotation & Motif Enrichment) ---
# Convert significant peaks to GRanges  (up-regulated and down-regulated) and set UCSC seqnames style
sig_peaks_up <- rownames(tt2)[tt2$Significant == TRUE & tt2$logFC > 0] %>% GRanges()
sig_peaks_down <- rownames(tt2)[tt2$Significant == TRUE & tt2$logFC < 0] %>% GRanges()

seqlevelsStyle(sig_peaks_up) <- "UCSC"
seqlevelsStyle(sig_peaks_down) <- "UCSC"

# Import ENCODE regulatory build features
reg_gr <- rtracklayer::import("homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20210107.gff.gz", format = "gff")

annotate_peaks <- function(query_gr, subject_gr, subject_annotate_col = "feature_type"){
    olaps <- findOverlaps(query = query_gr, subject = subject_gr)
    gr_no_hit <- query_gr[!overlapsAny(query_gr, subject_gr)]    
    gr_no_hit$feature <- "No feature"
    
    gr_hit <- query_gr[olaps@from]
    gr_hit$feature <- mcols(subject_gr)[olaps@to, subject_annotate_col]
    
    c(gr_no_hit, gr_hit)
}

peaks_up_annot <- annotate_peaks(query_gr = sig_peaks_up, subject_gr = reg_gr)
peaks_up_annot$direction <- "Up"

peaks_down_annot <- annotate_peaks(query_gr = sig_peaks_down, subject_gr = reg_gr)
peaks_down_annot$direction <- "Down"

all_peaks_annot <- c(peaks_up_annot, peaks_down_annot) %>% as.data.frame()

# Plot proportional genomic feature annotations
gg_annot_prop <- ggplot(all_peaks_annot, aes(x = direction, fill = feature)) +
    geom_bar(position = "fill") + 
    geom_text(aes(label = signif(..count.. / tapply(..count.., ..x.., sum)[as.character(..x..)], digits = 2)),
              stat = "count", position = position_fill(vjust = 0.5)) +
    labs(y = "Proportion", title = "Genomic Feature Annotation of DARs") +
    theme_bw()
gg_annot_prop

# Motif Enrichment Analysis using monaLisa & JASPAR2020
# Load the motif database (JASPAR), which defines the sequences to which transcription factors bind.
pwms <- getMatrixSet(JASPAR2020, opts = list(matrixtype = "PWM", tax_group = "vertebrates"))

# Resize the peaks to a uniform size for motif calling.
# It is best to divide them into windows of equal size, while ensuring the process remains computationally efficient.
# Therefore, the peaks are resized to a standard length (e.g., 100 bp) centered on the peak, as transcription factors (TFs) bind to short sequences.
peaks_up_reduce <- resize(sig_peaks_up, width = 100, fix = 'center')

# Retrieve the DNA sequences for the peaks
peaks_up_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peaks_up_reduce)

# Obtain the DNA sequence of chr22 as a genomic background for enrichment tests.
chr22_seq <- readDNAStringSet(filepath = "Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz")

# Calculate motif enrichments
# Perform a statistical test to determine whether motifs for specific TFs are significantly overrepresented in differential peaks (e.g., open peaks contain a specific NF-kappaB TF motif more frequently than would be expected based on the genome).
se <- calcBinnedMotifEnrR(seqs = peaks_up_seqs, pwmL = pwms,
                         genome.oversample = 1,
                         background = "genome", genome = chr22_seq)

# Narrow down to the most significant results for the pattern
sel <- apply(assay(se, "negLog10Padj"), 1, function(x) max(abs(x), 0, na.rm = TRUE)) > 1
seSel <- se[sel, ]

# Motif Enrichment Heatmap
plotMotifHeatmaps(x = seSel, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 2.0, cluster = TRUE, maxEnr = 2, maxSig = 10, 
                  show_motif_GC = TRUE)
