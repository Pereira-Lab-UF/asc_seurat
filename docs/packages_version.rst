:orphan:

.. _packages_version:

*******************
Session information
*******************

All R-package dependencies for the current Asc-Seurat release are
declared in ``DESCRIPTION`` and installed automatically when you run

.. code-block:: r

   pak::pkg_install("Pereira-Lab-UF/asc_seurat", dependencies = TRUE)

To capture the exact versions on your machine for a methods section or
a bug report, run ``sessionInfo()`` after launching the app:

.. code-block:: r

   ascseurat::run_app()
   # … in another R session, or after closing the app:
   sessionInfo()

See :doc:`references` for the list of packages Asc-Seurat directly
calls.
