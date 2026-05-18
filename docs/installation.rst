.. _installation:

************
Installation
************

Asc-Seurat v3 can be installed in two ways:

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

   docker pull kirstlab/asc_seurat:3
   docker run --rm -p 3838:3838 kirstlab/asc_seurat:3

Then open `http://localhost:3838 <http://localhost:3838>`_ in a browser.

The Docker image bundles every optional dependency — Scanpy/PAGA,
Monocle 3, ``PseudotimeDE``, ``tradeSeq``, ``SingleR``, ``celldex``,
``scDblFinder``, and AnnData support — so every Asc-Seurat feature is
available out of the box.

Mounting your own data
----------------------

To make a local folder visible inside the container, mount it onto
``/home/ascseurat/data`` (and optionally ``/home/ascseurat/RDS_files``):

.. code-block:: bash

   docker run --rm -p 3838:3838 \
     -v "$(pwd)/data":/home/ascseurat/data \
     -v "$(pwd)/RDS_files":/home/ascseurat/RDS_files \
     kirstlab/asc_seurat:3

Files placed in ``./data`` and ``./RDS_files`` then appear inside the app.

.. note::

   The v3 Docker image **does not** require mounting
   ``/var/run/docker.sock``. The Docker-in-Docker requirement that v2
   needed for Dynverse trajectory inference is gone.

Option 2 — R package from GitHub
================================

Requires **R ≥ 4.3.0**.

From a terminal:

.. code-block:: bash

   Rscript -e 'if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak", repos = "https://cloud.r-project.org"); pak::pkg_install("Pereira-Lab-UF/asc_seurat", dependencies = TRUE)'

Or from inside an R session:

.. code-block:: r

   install.packages("pak")
   pak::pkg_install("Pereira-Lab-UF/asc_seurat", dependencies = TRUE)

Then launch the app:

.. code-block:: r

   ascseurat::run_app()

Using ``dependencies = TRUE`` installs the optional R packages used by
Asc-Seurat features such as SingleR annotation, doublet detection,
AnnData input, and the trajectory DE engines (``PseudotimeDE``,
``scMaSigPro``, ``tradeSeq``).

Python dependencies for PAGA
----------------------------

PAGA uses the Python `Scanpy <https://scanpy.readthedocs.io/>`_ stack
via `reticulate <https://rstudio.github.io/reticulate/>`_. A one-liner
prepares and verifies a managed Python environment:

.. code-block:: r

   ascseurat::setup_paga()

Docker users can skip this step — Scanpy is already in the image.

Trajectory dependencies
-----------------------

Slingshot, PAGA, and Monocle 3 are the three supported trajectory
methods in v3, with method-specific gene-discovery engines.

- **Docker image** — Scanpy/PAGA, Monocle 3, ``PseudotimeDE``,
  ``scMaSigPro``, and ``tradeSeq`` are all preinstalled.
- **R package** — installed when you run ``pak::pkg_install(...,
  dependencies = TRUE)`` (see above). Run ``ascseurat::setup_paga()``
  once to prepare the Python side of PAGA.

Resource notes
==============

Single-cell analysis is memory-intensive. Larger datasets will require
more RAM whether you use Docker or the R package.

The RPCA integration workflow raises ``future.globals.maxSize``
dynamically to just below detected available memory. To force a
specific cap, set the ``ASCSEURAT_FUTURE_GLOBALS_MAXSIZE`` environment
variable to a byte value before launching the app.

Docker example — cap at 8 GB:

.. code-block:: bash

   docker run --rm -p 3838:3838 \
     -e ASCSEURAT_FUTURE_GLOBALS_MAXSIZE=8000000000 \
     kirstlab/asc_seurat:3

R example — cap at 8 GB:

.. code-block:: r

   Sys.setenv(ASCSEURAT_FUTURE_GLOBALS_MAXSIZE = "8000000000")
   ascseurat::run_app()

Verifying the installation
==========================

After launching the app, open the **Demo** tab in the top navigation.
A 2,000-cell subset of the PBMC 3k reference dataset is loaded
automatically and you can walk through the full workflow in a few
minutes. See :ref:`getting_started` for the step-by-step tour.
