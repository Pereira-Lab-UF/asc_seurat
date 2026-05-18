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

Use :guilabel:`Use de novo marker genes detected in Step 5` to visualize genes
directly from the integrated differential-expression results. When cluster
markers are available, choose the cluster and the number of top genes to plot.

Expression values
=================

The :guilabel:`Expression values` selector controls all generated plots:

* :guilabel:`Counts` uses raw counts.
* :guilabel:`Normalized` uses normalized expression.
* :guilabel:`Scaled` uses scaled expression and scales selected genes on demand
  if they were not already present in ``scale.data``.

Plots
=====

Click :guilabel:`Show Expression Plots` to generate the clustered dot plot,
feature plots, violin plots, and optional per-gene zip bundle. For integrated
objects, the plots should be interpreted together with the sample-colored UMAPs
and sample-aware cluster tables from the clustering step.

Download controls are collapsed by default. Expand
:guilabel:`Download options` to choose plot size, resolution, and file type.
