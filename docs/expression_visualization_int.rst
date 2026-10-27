.. _expression_visualization_int:

***********************************
Integrated expression visualization
***********************************

The integrated workflow uses the same visualization controls as the
single-sample workflow, but the clustered object contains all samples.

Gene sources
============

Use :guilabel:`Upload file` for a custom CSV or TSV marker list with one gene
per row and an optional second column for gene groups.

Use :guilabel:`Use de novo marker genes detected in Step 3` to visualize genes
directly from the integrated differential-expression results. When cluster
markers are available, choose the cluster and the number of top genes to plot.

Expression values
=================

The :guilabel:`Expression values` selector controls all generated plots:

* :guilabel:`Counts` uses raw counts.
* :guilabel:`Normalized` uses normalized expression.
* :guilabel:`Scaled` uses scaled expression. Selected genes that are not already
  in ``scale.data`` are scaled on demand, so genes that were already scaled keep
  their values. If ``scale.data`` was computed on other cells (for example,
  before a subset that kept the current clusters), all selected genes are
  scaled again on the current cells.

Plots
=====

Click :guilabel:`Show Expression Plots` to generate the clustered dot plot,
feature plots, violin plots, and optional per-gene zip bundle. For integrated
objects, the plots should be interpreted together with the sample-colored UMAPs
and sample-aware cluster tables from the clustering step.

The plots are cleared when the clustered object changes (re-clustering,
renaming clusters, or subsetting) or when the DE results they were built
from change; click :guilabel:`Show Expression Plots` again to plot the
current object. If no genes can be selected, a warning is shown.

Download controls are collapsed by default. Expand
:guilabel:`Download options` to choose plot size, resolution, and file type.
