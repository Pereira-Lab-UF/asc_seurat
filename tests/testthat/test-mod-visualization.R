test_that("visualization helpers fetch the requested expression layer", {
    obj <- suppressWarnings(clustered_fixture())
    genes <- head(rownames(obj), 5)

    normalized <- ascseurat:::fetch_expression_layer(obj, genes, "data")
    scaled <- ascseurat:::fetch_expression_layer(obj, genes, "scale.data")

    expect_equal(colnames(normalized), colnames(scaled))
    expect_false(isTRUE(all.equal(normalized, scaled)))
})

test_that("clustered expression dot plot is generated from selected layer", {
    obj <- suppressWarnings(clustered_fixture())
    genes <- head(rownames(obj), 5)

    plot <- ascseurat:::clustered_expression_dot_plot(obj, genes, "scale.data")

    expect_s3_class(plot, "ggplot")
})

test_that("clustered expression dot plot can split clusters by sample", {
    obj <- suppressWarnings(clustered_fixture())
    obj$samples <- rep(c("sample_A", "sample_B"), length.out = ncol(obj))
    genes <- head(rownames(obj), 5)

    expect_identical(ascseurat:::visualization_sample_col(obj), "samples")
    plot <- ascseurat:::clustered_expression_dot_plot(
        obj,
        genes,
        "data",
        sample_col = "samples"
    )

    split_labels <- levels(plot$data$cluster)
    expect_true(any(grepl("sample_A", split_labels, fixed = TRUE)))
    expect_true(any(grepl("sample_B", split_labels, fixed = TRUE)))
    expect_identical(plot$labels$x, "Cluster | Sample")
})
