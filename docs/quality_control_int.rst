.. _quality_control_int:

******************************
Integrated quality control
******************************

After multi-sample integration, Asc-Seurat shows integration QC metrics and an
optional additional filtering panel. This step is intended for post-integration
cleanup after you have already reviewed each sample individually.

The integrated QC violin plot shows:

* ``nFeature_RNA``: genes detected in each cell.
* ``nCount_RNA``: molecules detected in each cell.
* ``percent.mt``: percentage of transcripts matching the mitochondrial pattern.

Choose :guilabel:`No additional filtering` to cluster the integrated object as
created. Choose :guilabel:`Use filtering thresholds` to apply another round of
minimum genes, maximum genes, and maximum mitochondrial percentage cutoffs.

The integration tab also shows an elbow plot after integration. Use it to choose
the number of PCs for downstream clustering.

Download controls are collapsed by default. Expand
:guilabel:`Download options` to choose plot size, resolution, and file type.
