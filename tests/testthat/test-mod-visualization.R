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
