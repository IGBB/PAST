test_that("Data can be loaded and contains correct number of observations", {
    LD <- system.file("extdata", "LD.txt.gz",
                          package = "PAST2",
                          mustWork = TRUE
    )
    expect_equal(nrow(load_LD(LD, r_squared = "R.2")), 105714)
})
