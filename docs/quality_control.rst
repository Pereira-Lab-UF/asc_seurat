.. _quality_control:

***************
Quality control
***************

After loading a single-sample dataset, Asc-Seurat shows QC violin plots for:

* ``nFeature_RNA``: genes detected in each cell.
* ``nCount_RNA``: molecules detected in each cell.
* ``percent.mt``: percentage of transcripts matching the mitochondrial pattern.

Filtering mode
==============

The filtering panel has two modes:

* :guilabel:`No filtering` keeps all loaded cells for normalization.
* :guilabel:`Use filtering thresholds` applies minimum genes, maximum genes, and
  maximum mitochondrial percentage cutoffs.

The suggested threshold values are calculated from the loaded dataset as median
plus or minus three median absolute deviations (MADs), with safeguards for
valid ranges. These defaults are meant as a starting point, not a substitute for
reviewing the violin plots.

Doublet removal
===============

Doublet detection is independent from metric filtering. The
:guilabel:`Detect and remove doublets (scDblFinder)` checkbox is available in
both filtering modes:

* With :guilabel:`No filtering`, Asc-Seurat keeps all cells except predicted
  doublets when the checkbox is enabled and :guilabel:`Apply Filters` is run.
* With :guilabel:`Use filtering thresholds`, Asc-Seurat first applies the QC
  thresholds and then removes predicted doublets if requested.

If doublet removal is not selected and :guilabel:`No filtering` is active, you
can proceed directly to :guilabel:`Run Normalization + PCA`.

Downloads
=========

Plot download controls are collapsed by default. Expand
:guilabel:`Download options` to choose size, resolution, and file type.
