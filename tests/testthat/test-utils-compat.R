test_that("detect_input_format recognizes common formats", {
    csv_path <- tempfile(fileext = ".csv")
    tsv_path <- tempfile(fileext = ".tsv")
    rds_path <- tempfile(fileext = ".rds")

    write.csv(data.frame(gene = "A", cell = 1), csv_path, row.names = FALSE)
    write.table(data.frame(gene = "A", cell = 1), tsv_path,
                sep = "\t", row.names = FALSE)
    saveRDS(demo_fixture(), rds_path)

    expect_equal(ascseurat:::detect_input_format(csv_path), "csv")
    expect_equal(ascseurat:::detect_input_format(tsv_path), "tsv")
    expect_equal(ascseurat:::detect_input_format(rds_path), "rds")
    expect_equal(ascseurat:::detect_input_format(tempfile(fileext = ".foo")), "unknown")
})

test_that("format_bytes and local_future_globals_max_size return usable values", {
    expect_match(ascseurat:::format_bytes(1024), "KiB")
    info <- ascseurat:::local_future_globals_max_size()
    expect_type(info$target, "double")
    expect_true(info$target > 0)
})

test_that("QC mitochondrial default never starts at zero", {
    low_mito <- data.frame(percent.mt = c(0, 0.001, 0.002, 0.003, 0.004))
    zero_mito <- data.frame(percent.mt = rep(0, 5))

    expect_equal(ascseurat:::qc_filter_defaults(low_mito)$max_mito, 5)
    expect_equal(ascseurat:::qc_filter_defaults(zero_mito)$max_mito, 5)
})
