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

Docker containers cannot see arbitrary host paths unless those paths are
mounted when the container starts. The easiest pattern is to launch the
container from the folder that contains your data and mount that current
folder into Asc-Seurat's ``data/`` directory:

.. code-block:: bash

   cd "/path/to/folder/that/contains/your/data"
   docker run --rm -p 3838:3838 \
     -v "$PWD:/home/ascseurat/data:ro" \
     pereiralabbio/asc-seurat:3

Then enter paths relative to the app workdir, for example
``data/sample/filtered_feature_bc_matrix``. If the mounted folder itself
is the 10X matrix directory, enter ``data/``.

If you prefer to paste normal absolute paths from anywhere under your
home directory on macOS or Linux, mount your home directory at the same
path inside the container:

.. code-block:: bash

   docker run --rm -p 3838:3838 \
     -v "$HOME:$HOME:ro" \
     pereiralabbio/asc-seurat:3

Then enter the normal absolute path in Asc-Seurat, for example
``/Users/you/project/sample/filtered_feature_bc_matrix``, without
wrapping the path in quotes. To expose external drives on macOS, also
mount ``/Volumes:/Volumes:ro``. On Windows, mount a folder to a Linux
container path such as ``/home/ascseurat/data`` and enter ``data/...``
in the app.

Option 2 — R package from GitHub
================================

Requires **R ≥ 4.3.0** and Rstudio is recommended.

From inside an R session:

Follow the upstream `BPCells R installation instructions <https://github.com/bnprks/BPCells#r-installation>`_
first. In particular, install the HDF5 system dependency before installing
``BPCells``:

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt-get install -y libhdf5-dev pkg-config

   # macOS
   brew install hdf5 pkg-config

Then install ``BPCells`` before installing Asc-Seurat. ``BPCells`` is
required by the trajectory stack, and installing it first avoids a known
``pak`` failure with GitHub sub-directory packages.

.. code-block:: r

   install.packages(c("pak", "remotes"))
   remotes::install_github("bnprks/BPCells/r", upgrade = "never")
   pak::pkg_install(
     "Pereira-Lab-UF/asc_seurat",
     dependencies = c("Depends", "Imports", "LinkingTo")
   )

The ``dependencies`` argument intentionally installs only required runtime
packages. Do not use ``dependencies = TRUE`` for the default install unless
you also want optional ``Suggests`` packages; on macOS that can force source
builds such as ``PseudotimeDE`` that need a working Fortran toolchain.

If the GitHub install of ``BPCells`` is rate-limited or unavailable, use
the BPCells R-universe repository first, then run the Asc-Seurat install:

.. code-block:: r

   install.packages("BPCells", repos = c("https://bnprks.r-universe.dev", "https://cloud.r-project.org"))
   pak::pkg_install(
     "Pereira-Lab-UF/asc_seurat",
     dependencies = c("Depends", "Imports", "LinkingTo")
   )

Then launch the app:

.. code-block:: r

   ascseurat::run_app()

Or, from a terminal:

.. code-block:: bash

   Rscript -e 'if (!requireNamespace("pak", quietly = TRUE) || !requireNamespace("remotes", quietly = TRUE)) install.packages(c("pak", "remotes"), repos = "https://cloud.r-project.org"); remotes::install_github("bnprks/BPCells/r", upgrade = "never"); pak::pkg_install("Pereira-Lab-UF/asc_seurat", dependencies = c("Depends", "Imports", "LinkingTo"))'

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
