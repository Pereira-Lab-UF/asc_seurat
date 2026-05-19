# Trajectory inference with Slingshot, PAGA, and Monocle 3

Trajectory inference in Asc-Seurat v3 is offline-first:

- Slingshot is the default trajectory method once the R dependencies are
  installed.
- PAGA estimates cluster connectivity and supports connected-cluster
  gene discovery. Docker images include the required Python stack, and R
  package users can run
  [`ascseurat::setup_paga()`](https://pereira-lab-uf.github.io/asc_seurat/reference/setup_paga.md)
  to prepare it.
- Monocle 3 is available as an additional supported trajectory method
  when the `monocle3` package is installed, with graph-based
  trajectory-variable gene ranking.
- Trajectory differential expression is handled on Slingshot output
  through `PseudotimeDE-fast` as the recommended default, or through the
  alternate `scMaSigPro` and longer-running `tradeSeq` options.

Docker images include these trajectory dependencies. The R package
install also installs the required R trajectory packages; only the
Python stack for PAGA is prepared separately with
[`ascseurat::setup_paga()`](https://pereira-lab-uf.github.io/asc_seurat/reference/setup_paga.md).

Use a clustered Seurat object from the Single Sample or Integration
workflows as input to this tab.
