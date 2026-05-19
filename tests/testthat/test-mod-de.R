test_that("DE module finds markers for a clustered object", {
    obj <- clustered_fixture()

    shiny::testServer(
        ascseurat::mod_de_server,
        args = list(seurat_obj = shiny::reactive(obj), allow_conserved = FALSE),
        {
            session$setInputs(
                de_type = "all",
                test_use = "wilcox",
                pval_cutoff = 0.99,
                logfc_threshold = 0,
                min_pct = 0,
                only_pos = FALSE,
                run_de = 1
            )
            res <- markers()
            expect_s3_class(res, "data.frame")
            expect_gt(nrow(res), 0)
            expect_identical(colnames(res)[1:2], c("gene", "cluster"))
            expect_identical(rownames(res), as.character(seq_len(nrow(res))))
        }
    )
})

test_that("DE module compares metadata conditions inside a selected cluster", {
    obj <- clustered_fixture()
    cluster_counts <- table(as.character(Seurat::Idents(obj)))
    cluster <- names(cluster_counts)[cluster_counts >= 6][[1]]
    cluster_cells <- names(Seurat::Idents(obj))[as.character(Seurat::Idents(obj)) == cluster]

    obj$test_condition <- "other"
    obj$test_condition[cluster_cells] <- rep(c("A", "B"), length.out = length(cluster_cells))

    shiny::testServer(
        ascseurat::mod_de_server,
        args = list(seurat_obj = shiny::reactive(obj), allow_conserved = FALSE),
        {
            session$setInputs(
                de_type = "conditions",
                cluster_1 = cluster,
                condition_column = "test_condition",
                condition_value_1 = "A",
                condition_value_2 = "B",
                test_use = "wilcox",
                pval_cutoff = 0.999999,
                logfc_threshold = 0,
                min_pct = 0,
                only_pos = FALSE,
                run_de = 1
            )
            res <- markers()
            expect_s3_class(res, "data.frame")
            expect_identical(colnames(res)[1:2], c("gene", "cluster"))
            expect_false(anyNA(res$cluster))
            expect_false(anyNA(res$condition_column))
            expect_true(all(res$cluster == cluster))
            expect_true(all(res$condition_column == "test_condition"))
            expect_identical(rownames(res), as.character(seq_len(nrow(res))))
        }
    )
})
