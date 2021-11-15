test_that("Data can be loaded and contains correct number of observations", {
    single <- system.file("extdata", "single_gwas.txt.gz",
                          package = "PAST2",
                          mustWork = TRUE
    )
    expect_equal(nrow(load_GWAS_data(single,
                                     NULL,
                                     "single",
                                     "homozygous",
                                     trait = "Trait",
                                     marker = "Marker",
                                     locus = "Locus",
                                     site = "Site",
                                     p = "p",
                                     effect = "Effect")
    ), 3415)
})

test_that("All GWAS-loading functions return same results", {
    single <- system.file("extdata", "single_gwas.txt.gz",
                          package = "PAST2",
                          mustWork = TRUE
    )
    association <- system.file("extdata", "association.txt.gz",
                               package = "PAST2",
                               mustWork = TRUE
    )
    effects <- system.file("extdata", "effects.txt.gz",
                           package = "PAST2",
                           mustWork = TRUE
    )

    effects_single <- system.file("extdata", "effects.single_line.txt.gz",
                                  package = "PAST2", mustWork = TRUE
    )

    single_data <- load_GWAS_data(
        single,
        NULL,
        "single",
        "homozygous",
        trait = "Trait",
        marker = "Marker",
        locus = "Locus",
        site = "Site",
        p = "p",
        effect = "Effect"
    )
    two_data <- load_GWAS_data(
        association,
        effects_single,
        "two",
        "homozygous",
        trait = "Trait",
        marker = "Marker",
        locus = "Locus",
        site = "Site",
        p = "p",
        effect = "Effect"
    )
    TASSEL_data <- load_GWAS_data(
        association,
        effects,
        "TASSEL",
        "homozygous",
        trait = "Trait",
        marker = "Marker",
        locus = "Locus",
        site = "Site",
        p = "p",
        effect = "Effect"
    )
    expect_equal(single_data, two_data)
    expect_equal(TASSEL_data, two_data)
    expect_equal(TASSEL_data, single_data)
})
