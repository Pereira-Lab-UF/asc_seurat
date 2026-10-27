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
* :guilabel:`Scaled` uses scaled expression. Selected genes that are not already
  in ``scale.data`` are scaled on demand, regressing out the same variables as
  during normalization, so genes that were already scaled keep their values.
  If ``scale.data`` was computed on other cells (for example, before a subset
  that kept the current clusters), all selected genes are scaled again on the
  current cells.

Plots
=====

Click :guilabel:`Show Expression Plots` to generate:

* a clustered dot plot showing mean expression by cluster and percent of cells
  expressing each gene;
* UMAP feature plots for the selected genes;
* violin plots by cluster;
* an optional per-gene zip bundle.

The plots are cleared when the clustered object changes (re-clustering,
renaming clusters, or subsetting) or when the DE results they were built
from change; click :guilabel:`Show Expression Plots` again to plot the
current object. If no genes can be selected, a warning is shown.

Download controls are collapsed by default. Expand
:guilabel:`Download options` to choose plot size, resolution, and file type.
