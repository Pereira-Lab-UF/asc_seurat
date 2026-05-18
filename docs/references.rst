.. _references:

**********
References
**********

Asc-Seurat is built on the work of many other people and relies on a
diversity of R packages. These packages, in turn, have many dependencies.
Here we list all packages that Asc-Seurat v3 directly calls.

Analytical core
===============

* **Seurat**

  * web page: https://satijalab.org/seurat/
  * Publications:

    * `Satija, R. et al. (2015) <https://www.nature.com/articles/nbt.3192>`_ Spatial reconstruction of single-cell gene expression data. Nat. Biotechnol., 33, 495–502.
    * `Stuart, T. et al. (2019) <https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8>`_ Comprehensive Integration of Single-Cell Data. Cell, 177, 1888–1902.e21.
    * `Hao, Y. et al. (2024) <https://doi.org/10.1038/s41592-024-02276-1>`_ Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nat. Methods, 21, 1–9.

* **Slingshot** (trajectory inference)

  * web page: https://bioconductor.org/packages/release/bioc/html/slingshot.html
  * Publication:

    * `Street, K. et al. (2018) <https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-4772-0>`_ Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics. BMC Genomics, 19, 477.

* **scDblFinder** (doublet detection)

  * web page: https://bioconductor.org/packages/release/bioc/html/scDblFinder.html
  * Publication:

    * `Germain, P.-L. et al. (2022) <https://doi.org/10.12688/f1000research.73600.2>`_ Doublet identification in single-cell sequencing data using scDblFinder. F1000Res., 10, 979.

* **SingleR** (cell-type annotation)

  * web page: https://bioconductor.org/packages/release/bioc/html/SingleR.html
  * Publication:

    * `Aran, D. et al. (2019) <https://www.nature.com/articles/s41590-018-0276-y>`_ Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nat. Immunol., 20, 163–172.

* **scMaSigPro** (trajectory differential expression)

  * web page: https://github.com/BioBam/scMaSigPro
  * Publication:

    * `Prieto, C. & Martínez-García, P.M. (2023) <https://doi.org/10.1093/bioinformatics/btad681>`_ scMaSigPro: a polynomial regression framework for single-cell RNA sequencing time series. Bioinformatics, 39, btad681.

* **harmony** (batch integration)

  * web page: https://github.com/immunogenomics/harmony
  * Publication:

    * `Korsunsky, I. et al. (2019) <https://www.nature.com/articles/s41592-019-0619-0>`_ Fast, sensitive and accurate integration of single-cell data with Harmony. Nat. Methods, 16, 1289–1296.

* **scCustomize**

  * web page: https://samuel-marsh.github.io/scCustomize/

Additional packages
===================

CRAN
----

* bslib: https://rstudio.github.io/bslib/
* DT: https://rstudio.github.io/DT/
* dplyr: https://dplyr.tidyverse.org/
* ggplot2: https://ggplot2.tidyverse.org/
* hdf5r: https://cran.r-project.org/web/packages/hdf5r/
* metap: https://cran.r-project.org/web/packages/metap/index.html
* patchwork: https://patchwork.data-imaginist.com/
* reactable: https://glin.github.io/reactable/
* reticulate: https://rstudio.github.io/reticulate/
* sass: https://rstudio.github.io/sass/
* shiny: https://shiny.posit.co/
* shinycssloaders: https://cran.r-project.org/web/packages/shinycssloaders/
* shinyFeedback: https://merlinoa.github.io/shinyFeedback/
* shinyWidgets: https://dreamrs.github.io/shinyWidgets/
* SeuratObject: https://cran.r-project.org/web/packages/SeuratObject/

Bioconductor
------------

* celldex: https://bioconductor.org/packages/release/data/experiment/html/celldex.html
* SingleCellExperiment: https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html
* slingshot: https://bioconductor.org/packages/release/bioc/html/slingshot.html

.. note::

   **Packages removed in v3:**
   Asc-Seurat v2.x used **Dynverse** (dynplot, dynwrap, dynfeature) for trajectory
   inference and **biomaRt** for online gene-ID mapping. Dynverse has been
   replaced with Slingshot, PAGA, and Monocle 3; bundled GO enrichment has been
   removed to keep the app lighter. See :ref:`migrating` for details.
