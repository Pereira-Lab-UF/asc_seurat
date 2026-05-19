.. _getting_started:

***************
Getting started
***************

The fastest way to evaluate Asc-Seurat is to use the built-in **Demo**
tab, which auto-loads a small PBMC dataset and lets you walk through the
whole single-sample workflow without supplying any data of your own.

Launching the app
=================

If you installed the R package:

.. code-block:: r

   ascseurat::run_app()

If you are using the Docker image:

.. code-block:: bash

   docker pull pereiralabbio/asc-seurat:3 && \\
   docker run --rm -d --name asc-seurat -p 3838:3838 pereiralabbio/asc-seurat:3 && open http://localhost:3838

then open ``http://localhost:3838`` in a browser.

For your own data, launch Docker from the folder that contains the data
and mount that folder into the app:

.. code-block:: bash

   cd "/path/to/folder/that/contains/your/data"
   docker run --rm -p 3838:3838 \
     -v "$PWD:/home/ascseurat/data:ro" \
     pereiralabbio/asc-seurat:3

Then use paths such as ``data/sample/filtered_feature_bc_matrix`` in the
app. See :ref:`installation` for additional mount examples.

.. figure:: images/v3/v3_home.png
   :alt: Asc-Seurat home screen.
   :width: 100%
   :align: center

   The home screen. Click **Try it with Demo Data** or use the top
   navigation to jump straight to a workflow.

The Demo tab
============

The **Demo** tab auto-loads a 2,000-cell subset of the PBMC 3k
reference dataset that ships with the package. It is the recommended
starting point for first-time users: the dataset is small enough to
finish every step in a few minutes while still exercising QC,
normalization, clustering, differential expression, and visualization.

Suggested first run:

1. Open the app.
2. Click **Try it with Demo Data** on the home page (or the **Demo**
   tab in the navigation bar). The PBMC object is loaded automatically.
3. Step through QC, Normalize & Cluster, and DE / Visualization with
   the default settings.
4. From the DE results table, send a gene set directly into the
   visualization module to inspect markers without re-importing a CSV.

.. figure:: images/v3/v3_demo.png
   :alt: Demo tab with example dataset loaded.
   :width: 100%
   :align: center

   The Demo tab. The example dataset is auto-loaded and every step
   (Load, QC, Normalize & Cluster, Explore) is available in a single
   scroll.

Top-level navigation
====================

The navigation bar groups the workflows into:

- **Home** — landing page with the Demo button and a quick-start guide.
- **Single Sample** — load, QC, cluster, and explore one sample.
- **Integration** — declare multiple samples inline, set per-sample QC,
  and integrate with RPCA or Harmony.
- **Trajectory Inference** — Slingshot, PAGA, and Monocle 3 in
  one tab.
- **Tools** (dropdown) — Cell-type Annotation (SingleR) and Advanced
  Plots (stacked violin / multi-gene dot plot).
- **Demo** — the auto-loaded PBMC walkthrough described above.

Working with your own data
==========================

Once you are comfortable with the demo, load your own dataset:

- For **a single sample**, see :ref:`loading_data`.
- For **multiple samples (integration)**, see :ref:`loading_data_int`.
- For **trajectory inference**, see :ref:`trajectory_inference`.
- For **cell-type annotation**, see :doc:`annotation`.
