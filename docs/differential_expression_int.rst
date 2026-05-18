.. _differental_expression_int:

***********************************************************
Markers identification and differential expression analysis
***********************************************************

After integrated clustering, Asc-Seurat v3 exposes marker testing in
**Step 3: Differential Expression / Marker Identification (Optional)**.
The integrated workflow supports the same core options as the
single-sample workflow, plus sample-aware tests for multi-sample
objects.

Selected marker genes can be sent directly from the results table to
the integrated visualization module (see
:ref:`expression_visualization_int`). You can also download the results
as a CSV file for external use.

Analysis modes
==============

The integrated DE panel can:

- find markers for all clusters;
- find markers for one selected cluster;
- compare two selected clusters;
- test differential expression between conditions within a selected
  cluster;
- find conserved markers across samples for a selected cluster.

All modes share the statistical test selector, adjusted p-value cutoff,
log2 fold-change threshold, minimum fraction of cells expressing the
gene, and the option to return only positive markers.

.. figure:: images/v3/v3_de_int_all_clusters.png
   :alt: Asc-Seurat v3 integrated DE settings for all clusters.
   :width: 100%
   :align: center

   Settings for finding marker genes for all clusters in the integrated
   object.

.. figure:: images/v3/v3_de_int_one_cluster.png
   :alt: Asc-Seurat v3 integrated DE settings for one cluster.
   :width: 100%
   :align: center

   Settings for finding markers for one selected cluster.

.. figure:: images/v3/v3_de_int_compare_clusters.png
   :alt: Asc-Seurat v3 integrated DE settings for comparing clusters.
   :width: 100%
   :align: center

   Settings for comparing two selected clusters.

.. figure:: images/v3/v3_de_int_conditions.png
   :alt: Asc-Seurat v3 integrated DE settings for condition testing.
   :width: 100%
   :align: center

   Settings for testing differential expression between sample
   conditions within a selected cluster.

.. figure:: images/v3/v3_de_int_conserved.png
   :alt: Asc-Seurat v3 integrated DE settings for conserved markers.
   :width: 100%
   :align: center

   Settings for finding conserved markers across samples.

Marker table
============

After the search runs, Asc-Seurat displays an interactive marker table.
The table can be searched, paged, and downloaded as CSV. The selected
genes can also be passed directly to **Step 4: Gene Expression
Visualization** without re-uploading a marker list.

.. figure:: images/v3/v3_de_int_marker_table.png
   :alt: Asc-Seurat v3 integrated marker table.
   :width: 100%
   :align: center

   Marker table generated for cluster 0 of the integrated WT and rhd6
   example samples.
