:orphan:

.. _packages_version:

*******************
Session information
*******************

All R-package dependencies for the current Asc-Seurat release are
declared in ``DESCRIPTION``. Install ``BPCells`` first, then install
Asc-Seurat:

.. code-block:: r

   install.packages(c("pak", "remotes"))
   remotes::install_github("bnprks/BPCells/r", upgrade = "never")
   pak::pkg_install("Pereira-Lab-UF/asc_seurat", dependencies = TRUE)

This explicit ``BPCells`` step avoids a known ``pak`` failure with the
GitHub sub-directory package used by BPCells.

To capture the exact versions on your machine for a methods section or
a bug report, run ``sessionInfo()`` after launching the app:

.. code-block:: r

   ascseurat::run_app()
   # … in another R session, or after closing the app:
   sessionInfo()

See :doc:`references` for the list of packages Asc-Seurat directly
calls.
