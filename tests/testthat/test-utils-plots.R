test_that("save_plot writes a non-empty png", {
    plot_obj <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
    out <- tempfile(fileext = ".png")

    ascseurat:::save_plot(
        plot = plot_obj,
        filename = tools::file_path_sans_ext(basename(out)),
        format = "png",
        width = 10,
        height = 8,
        dpi = 150,
        path = dirname(out)
    )

    expect_true(file.exists(out))
    expect_gt(file.info(out)$size, 0)
})
