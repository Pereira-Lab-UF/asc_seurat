.. _installation:

************
Installation
************

Asc-Seurat can be installed in two ways:

- **Docker** — recommended for most users. Zero dependency management;
  every optional feature is bundled.
- **R package from GitHub** — lighter and starts faster if you already
  work in R and are comfortable managing Bioconductor packages.

Both options run the same Shiny app and produce identical results.

Source code lives at https://github.com/Pereira-Lab-UF/asc_seurat.

Option 1 — Docker (recommended)
===============================

Requirements: a recent `Docker <https://docs.docker.com/get-docker/>`_
installation. Nothing else.

.. code-block:: bash

   docker pull pereiralabbio/asc-seurat:3 && docker run --rm -p 3838:3838 pereiralabbio/asc-seurat:3

Then open `http://localhost:3838 <http://localhost:3838>`_ in a browser.

Option 2 — R package from GitHub
================================

Requires **R ≥ 4.3.0** and Rstudio is recommended.

From inside an R session:

.. code-block:: r

   install.packages("pak")
   pak::pkg_install("Pereira-Lab-UF/asc_seurat", dependencies = TRUE)

Then launch the app:

.. code-block:: r

   ascseurat::run_app()

Or, from a terminal:

.. code-block:: bash

   Rscript -e 'if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak", repos = "https://cloud.r-project.org"); pak::pkg_install("Pereira-Lab-UF/asc_seurat", dependencies = TRUE)'

Python dependencies for PAGA
----------------------------

Docker users can skip this step — Scanpy is already in the image.

PAGA uses the Python `Scanpy <https://scanpy.readthedocs.io/>`_ stack
via `reticulate <https://rstudio.github.io/reticulate/>`_. A one-liner
prepares and verifies a managed Python environment:

.. code-block:: r

   ascseurat::setup_paga()

Resource notes
==============

Single-cell analysis is memory-intensive. Larger datasets will require
more RAM whether you use Docker or the R package. we recommend at least 8 GB for datasets of ~20,000 cells.

Verifying the installation
==========================

After launching the app, open the **Demo** tab in the top navigation.
A 2,000-cell subset of the PBMC 3k reference dataset is loaded
automatically and you can walk through the full workflow in a few
minutes. See :ref:`getting_started` for the step-by-step tour.
