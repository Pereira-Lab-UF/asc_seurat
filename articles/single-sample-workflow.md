# Single-sample workflow

The single-sample workflow follows the published Asc-Seurat structure:

1.  Load data or demo data.
2.  Run QC and optionally remove doublets.
3.  Normalize, score cell cycle if needed, and inspect the PCA elbow
    plot.
4.  Cluster the data, optionally compute tSNE, and rename clusters.
5.  Run differential expression and marker discovery.
6.  Generate feature, violin, and heatmap visualizations.
