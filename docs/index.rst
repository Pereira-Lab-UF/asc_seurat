Asc-Seurat documentation
========================

**Asc-Seurat** is a web application for single-cell RNA-seq analysis,
built around `Seurat v5 <https://satijalab.org/seurat/>`_. It walks you
through the full scRNA-seq workflow — quality control, normalization,
clustering, differential expression, trajectory inference, cell-type
annotation, and publication-ready visualization — without writing any
R code.

This documentation covers Asc-Seurat **v3**, the current release. v3
builds on the original Asc-Seurat published in *BMC Bioinformatics*
(2021) and keeps the same end-to-end click-driven workflow.

.. figure:: images/v3/v3_home.png
   :alt: Asc-Seurat home screen.
   :width: 100%
   :align: center

   Asc-Seurat home screen.

What you can do with Asc-Seurat
-------------------------------

- **Load data from many formats** — 10X Genomics directories, 10X HDF5, CSV/TSV count matrices, AnnData (``.h5ad``), and existing Seurat objects.
- **Quality control and filtering** — interactive violin plots, per-cell metric thresholds, and optional doublet removal.
- **Normalization and clustering** — ``LogNormalize`` or ``SCTransform``, UMAP and t-SNE embeddings, and an inline control to rename clusters with biological labels.
- **Differential expression** — cluster markers, between-cluster DE, and between-sample DE for integrated datasets.
- **Multi-sample integration** — declare samples inline, set per-sample QC, and integrate with RPCA (Seurat) or Harmony.
- **Trajectory inference** — three supported methods: Slingshot, PAGA, and Monocle 3, plus method-specific gene discovery.
- **Cell-type annotation** — automated labelling against reference atlases via SingleR.
- **Advanced plots** — stacked violin and multi-gene dot plots, with control over gene and cluster ordering.
- **Per-gene plot bundles** — download a zipped collection of every plot for a selected gene in one click.
- **Bookmarking and session reports** — capture an analysis state and reproduce or share it later.

Where to start
--------------

- **First time?** Install Asc-Seurat (:ref:`installation`) and then
  follow the :ref:`getting_started` walkthrough. The built-in **Demo**
  tab auto-loads a 2,000-cell PBMC dataset so you can try every step
  without supplying your own data.
- **Ready to load your data?** Jump to
  :ref:`Single-sample loading <loading_data>`,
  :ref:`Integration loading <loading_data_int>`, or
  :doc:`trajectory_inference`.

.. toctree::
   :hidden:
   :caption: General information
   :maxdepth: 2

   installation
   getting_started
   references
   license

.. toctree::
   :hidden:
   :caption: Analysis of individual sample
   :maxdepth: 4

   loading_data
   quality_control
   clustering
   differential_expression
   expression_visualization

.. toctree::
   :hidden:
   :caption: Analysis of multiple samples
   :maxdepth: 4

   loading_data_int
   quality_control_int
   clustering_int
   differential_expression_int
   expression_visualization_int

.. toctree::
   :hidden:
   :caption: Trajectory inference
   :maxdepth: 4

   trajectory_inference

.. toctree::
   :hidden:
   :caption: Cell-type annotation
   :maxdepth: 4

   annotation

.. toctree::
   :hidden:
   :caption: Advanced plots
   :maxdepth: 4

   Advanced_plots
