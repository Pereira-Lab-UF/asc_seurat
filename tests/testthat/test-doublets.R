test_that("doublet detection adds scDblFinder metadata when available", {
    skip_if_not_installed("scDblFinder")

    # scDblFinder expects normalized data; pre-normalize so 'data' layer is populated
    obj <- Seurat::NormalizeData(demo_fixture(), verbose = FALSE)
    dbl_obj <- ascseurat:::detect_doublets(obj)

    expect_true(all(c("scDblFinder_class", "scDblFinder_score") %in% colnames(dbl_obj@meta.data)))
    expect_true(all(!is.na(dbl_obj$scDblFinder_class)))
    expect_true(is.numeric(dbl_obj$scDblFinder_score))
})

test_that("doublet detection works after monocle3 is loaded", {
    skip_if_not_installed("scDblFinder")
    skip_if_not_installed("monocle3")
    loadNamespace("monocle3")
    attached_before <- search()

    obj <- Seurat::NormalizeData(demo_fixture(), verbose = FALSE)
    dbl_obj <- ascseurat:::detect_doublets(obj, sample_col = "orig.ident")

    expect_true(all(!is.na(dbl_obj$scDblFinder_class)))
    expect_identical(search(), attached_before)
})
