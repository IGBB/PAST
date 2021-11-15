test_that(
    "Data can be loaded and contains correct number of observations", {
        gwas_data <- load_GWAS_data(
            system.file(
                "extdata",
                "single_gwas.txt.gz",
                package = "PAST2",
                mustWork = TRUE
            )
            ,
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
        LD <- load_LD(
            system.file(
                "extdata",
                "LD.txt.gz",
                package = "PAST2",
                mustWork = TRUE
            ),
            r_squared = "R.2"
        )
        genes <- assign_SNPs_to_genes(
            gwas_data,
            LD,
            system.file(
                "extdata",
                "genes.gff.gz",
                package = "PAST2",
                mustWork = TRUE
            ),
            "gene",
            1000,
            0.8,
            "Name"
        )

        enrichment_data <- find_pathway_significance(
            genes,
            system.file(
                "extdata",
                "pathways.txt.gz",
                package = "PAST2",
                mustWork = TRUE
            ),
            5,
            "increasing",
            1000
        )
        expect_equal(nrow(enrichment_data), 17)
        expect_lte(enrichment_data[1]$`p-value`, 0.2)
        expect_gte(enrichment_data[17]$`p-value`, 0.85)
    })
