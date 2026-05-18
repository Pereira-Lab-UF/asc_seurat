test_that("doublet detection adds scDblFinder metadata when available", {
    skip_if_not_installed("scDblFinder")

    # scDblFinder expects normalized data; pre-normalize so 'data' layer is populated
    obj <- Seurat::NormalizeData(demo_fixture(), verbose = FALSE)
    dbl_obj <- ascseurat:::detect_doublets(obj)

    expect_true(all(c("scDblFinder_class", "scDblFinder_score") %in% colnames(dbl_obj@meta.data)))
    expect_true(all(!is.na(dbl_obj$scDblFinder_class)))
    expect_true(is.numeric(dbl_obj$scDblFinder_score))
})
