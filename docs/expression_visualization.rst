.. _expression_visualization:

************************
Expression visualization
************************

Step 6 visualizes marker genes on the clustered object. Genes can come from an
uploaded marker list or from the de novo marker genes detected in Step 5.

Gene sources
============

Use :guilabel:`Upload file` for a custom list. The file may be CSV or TSV and
should contain one gene per row. A second column is optional and is used as a
gene group selector.

Use :guilabel:`Use de novo marker genes detected in Step 5` to pull genes
directly from the differential-expression table. Choose the cluster and the
number of top genes to visualize; genes are ranked by adjusted p-value.

Expression values
=================

The :guilabel:`Expression values` selector controls the layer used for the
clustered dot plot, feature plots, violin plots, and per-gene bundle:

* :guilabel:`Counts` uses raw counts.
* :guilabel:`Normalized` uses normalized expression.
* :guilabel:`Scaled` uses scaled expression. If a selected gene was not already
  present in ``scale.data``, Asc-Seurat scales the selected genes on demand for
  plotting.

Plots
=====

Click :guilabel:`Show Expression Plots` to generate:

* a clustered dot plot showing mean expression by cluster and percent of cells
  expressing each gene;
* UMAP feature plots for the selected genes;
* violin plots by cluster;
* an optional per-gene zip bundle.

Download controls are collapsed by default. Expand
:guilabel:`Download options` to choose plot size, resolution, and file type.
