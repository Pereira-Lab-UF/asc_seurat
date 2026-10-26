test_that("QC preview cache key changes when a 10X file is overwritten in place", {
    sample_dir <- withr::local_tempdir()
    matrix_file <- file.path(sample_dir, "matrix.mtx")
    writeLines("old", matrix_file)

    shiny::testServer(ascseurat::mod_integration_server, {
        before <- input_version_key(sample_dir)
        writeLines("new contents", matrix_file)
        expect_false(identical(input_version_key(sample_dir), before))
    })
})
