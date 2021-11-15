test_that(
    "SNP-gene assignment contains correct number of genes in correct order", {
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
        expect_equal(nrow(genes), 183)
        expect_equal(genes[1]$marker, "S1_17462583")
        expect_equal(genes[183]$marker, "S9_82727414")
    })
