test_that("paga availability helper returns a logical", {
    expect_type(ascseurat:::paga_available(), "logical")
    expect_length(ascseurat:::paga_available(), 1)
})

test_that("PAGA Python requirements include the supported Scanpy stack", {
    expect_true(all(c("scanpy", "anndata", "numpy", "scipy", "pandas", "leidenalg", "igraph") %in%
        ascseurat:::paga_python_packages()))
})
