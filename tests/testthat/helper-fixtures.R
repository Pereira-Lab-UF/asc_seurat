demo_fixture <- local({
    cached <- NULL
    function() {
        if (is.null(cached)) {
            obj <- readRDS(system.file("extdata", "pbmc_demo.rds", package = "ascseurat"))
            # Ensure percent.mt is present so regression tests can use it
            if (!"percent.mt" %in% colnames(obj@meta.data)) {
                obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(
                    obj, pattern = "^(MT-|ATMG|ATCG)"
                )
            }
            cached <<- obj
        }
        cached
    }
})

clustered_fixture <- local({
    cached <- NULL
    function() {
        if (is.null(cached)) {
            obj <- ascseurat:::run_normalization(
                demo_fixture(),
                method = "lognorm",
                nfeatures = 500
            )
            cached <<- ascseurat:::run_clustering(
                obj,
                dims = 1:10,
                resolution = 0.4
            )
        }
        cached
    }
})

clustered_fixture_path <- local({
    cached_path <- NULL
    function() {
        if (is.null(cached_path) || !file.exists(cached_path)) {
            cached_path <<- tempfile("ascseurat_clustered_", fileext = ".rds")
            saveRDS(clustered_fixture(), cached_path)
        }
        cached_path
    }
})

markers_fixture_path <- function() {
    testthat::test_path("..", "test_markers.csv")
}
