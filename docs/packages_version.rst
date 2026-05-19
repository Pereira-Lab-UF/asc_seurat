:orphan:

.. _packages_version:

*******************
Session information
*******************

All R-package dependencies for the current Asc-Seurat release are
declared in ``DESCRIPTION``. Follow the upstream
`BPCells R installation instructions <https://github.com/bnprks/BPCells#r-installation>`_
first so the HDF5 system dependency is available, then install
``BPCells`` before installing Asc-Seurat:

.. code-block:: r

   install.packages(c("pak", "remotes"))
   remotes::install_github("bnprks/BPCells/r", upgrade = "never")
   pak::pkg_install(
     "Pereira-Lab-UF/asc_seurat",
     dependencies = c("Depends", "Imports", "LinkingTo")
   )

This explicit ``BPCells`` step avoids a known ``pak`` failure with the
GitHub sub-directory package used by BPCells. The explicit ``dependencies``
value installs every package declared as an app runtime dependency in
``DESCRIPTION``, including ``PseudotimeDE``; it avoids only developer,
documentation, and test-only packages. On macOS, ``PseudotimeDE`` requires
the official R GNU Fortran toolchain from
`R for macOS tools <https://mac.r-project.org/tools/>`_.

To capture the exact versions on your machine for a methods section or
a bug report, run ``sessionInfo()`` after launching the app:

.. code-block:: r

   ascseurat::run_app()
   # … in another R session, or after closing the app:
   sessionInfo()

See :doc:`references` for the list of packages Asc-Seurat directly
calls.
