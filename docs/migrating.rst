.. _migrating:

*******************
Migrating from v2.x
*******************

Asc-Seurat v3 keeps the same high-level workflow structure but modernizes
the analytical backend and the way the app is installed. This page
summarizes what carries over, what has changed, and what you may need
to do to move an existing v2 project forward.

What stays the same
===================

- The published Asc-Seurat navigation — **Single Sample**, **Integration**,
  **Trajectory**, **Annotation**, **Advanced Plots** — is preserved.
- Input formats are unchanged: 10x-style MTX directories for raw data,
  CSV configuration files for integration, CSV/TSV gene lists for
  visualization.
- **Old Seurat RDS files still load.** Asc-Seurat v3 calls
  ``UpdateSeuratObject()`` automatically when you open an object saved
  with Seurat v3 or v4.

What changed
============

Installation and launch
-----------------------

Asc-Seurat is now a proper R package:

.. code-block:: r

   install.packages("pak")
   pak::pkg_install("Pereira-Lab-UF/asc_seurat", dependencies = TRUE)
   ascseurat::run_app()

The Docker image is still available and is the recommended option if you
don't want to manage R/Bioconductor dependencies yourself:

.. code-block:: bash

   docker pull kirstlab/asc_seurat:3
   docker run --rm -p 3838:3838 kirstlab/asc_seurat:3

The v3 Docker image **no longer requires mounting** ``/var/run/docker.sock`` —
that Docker-in-Docker requirement (needed by Dynverse in v2) is gone.

Seurat v5
---------

Asc-Seurat v3 uses Seurat v5 throughout. Objects created by older versions
are upgraded on load. There is nothing you need to do manually.

Trajectory inference
--------------------

- Dynverse / ``dyno`` have been removed. They required a Docker-in-Docker
  runtime, were slow, and frequently broke between releases.
- Three supported trajectory methods are now available in the Shiny
  workflow: **Slingshot** (default), **PAGA**, and **Monocle 3**.
- Trajectory gene discovery is method-specific. Slingshot uses
  trajectory differential-expression engines, PAGA reports connected-state
  markers and pseudotime-associated genes, and Monocle 3 reports
  graph-variable genes.
- The Docker image bundles all of these. R package users can install
  PAGA dependencies with ``ascseurat::setup_paga()``, install
  ``monocle3`` from CRAN/Bioconductor, and install ``PseudotimeDE`` from
  ``https://github.com/dsong-lab/PseudotimeDE.git``.

**Expect slightly different lineage assignments** compared with v2,
because the method stack has changed. You may need to re-choose
root/start and end clusters in the trajectory tab.

Integration
-----------

- **RPCA integration** (Seurat v5 ``IntegrateLayers``) is the default
  integration workflow. Harmony is also available.
- The older CCA-based anchor integration used by Asc-Seurat v2 is no
  longer the default, though existing integrated RDS objects will still
  load.

UI and workflow
---------------

- New **Demo** tab that auto-loads a 2,000-cell PBMC dataset. See
  :doc:`getting_started`.
- Top-level **Tools** menu groups the secondary workflows (Cell-type
  Annotation and Advanced Plots) and keeps the
  primary navigation focused on the analysis pipeline.
- **Integration** tab now uses an inline per-sample table (no separate
  config CSV) with per-sample QC parameters.
- **Trajectory** tab supports an upload-from-browser **and** a
  path-on-server option for very large RDS files.
- New **cluster renaming** UI in the clustering tab.
- **DE → visualization** direct hand-off: selected DE genes can be sent
  straight into the visualization module.
- **Session reports** and **bookmarking** are built in.
- **Doublet detection** via ``scDblFinder`` is available in the QC step.
- **SingleR** cell-type annotation is included.
- **Per-gene plot-bundle export**: download a zipped collection of all
  plots for a selected gene.
- H5AD (AnnData) input via ``anndataR``.

Removed
-------

- ``dynverse`` / ``dyno`` / ``dynplot`` / ``dynwrap`` / ``dynfeature``
  (trajectory inference).
- ``biomaRt`` (online annotation).
- Functional annotation / GO enrichment has been removed from the app;
  use a dedicated enrichment web service or standalone R workflow.
- ``rclipboard`` copy-to-clipboard buttons (unused in practice).

See :doc:`references` for the full list of packages that v3 directly
depends on.

Action items for existing users
===============================

1. **Uninstall the old Docker image** if disk space matters to you —
   ``docker image rm kirstlab/asc_seurat:2`` (keep it if you still need
   to refer back to the v2 behavior).
2. **Re-run trajectory inference** on existing projects where lineage
   assignments matter, since the method stack has changed.
3. Pick the trajectory method that fits your data: **Slingshot** for the
   default lineage analysis (with the PseudotimeDE-fast / scMaSigPro /
   tradeSeq DE engines), **PAGA** when working from a clustered neighbor
   graph, or **Monocle 3** when you want a learned principal graph on the
   existing UMAP. R package users should run ``ascseurat::setup_paga()``
   before using PAGA, and install ``monocle3`` and ``PseudotimeDE`` if
   they want those methods locally (the Docker image bundles all of
   them).
