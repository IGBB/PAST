#' Load Linkage Disequilibrium
#'
#' @param LD_file the file containing linkage disequilibrium data
#' @param locus the name of the column containing the locus or chromosome
#' @param position1 the name of the column containing the first SNP's position
#' @param site1 the name of the column containing the first SNP's site
#' @param position2 the name of the column containing second SNP's position
#' @param site2 the name of the column containing the second SNP's site
#' @param distance the name of the column containing distance between the SNPs
#' @param r_squared the name of the column containing R^2
#' @return the linkage disequilibrium data in a data.table
#' @export
#' @import data.table
#' @examples
#' LD_file = system.file("extdata","LD.txt.gz", package = "PAST2",
#' mustWork = TRUE)
#' LD <- load_LD(LD_file, r_squared="R.2")
load_LD <- function(LD_file,
                    locus = "Locus1",
                    position1 = "Position1",
                    site1 = "Site1",
                    position2 = "Position2",
                    site2 = "Site2",
                    distance = NULL,
                    r_squared = "R^2") {


    arguments <- list(
        "locus" = locus,
        "position1" = position1,
        "site1" = site1,
        "position2" = position2,
        "site2" = site2,
        "distance" = distance,
        "r_squared" = r_squared
    )

    arguments <- arguments[!sapply(arguments,is.null)]

    argument_names <- names(arguments)
    LD_columns <- arguments[
        argument_names %in% argument_names[argument_names != ""]
    ]

    LD <- load_file(LD_file, LD_columns)

    # if distance wasn't provided, calculate it
    if (is.null(distance)) {
        LD[, distance := abs(position1 - position2)]
    }

    # pull the sites and remove duplicates, then calculate a new site designation
    sites = data.table::as.data.table(unique(LD[,site1]))
    data.table::setnames(sites, "V1", "site_original")
    sites[,site := as.numeric(.I)]

    # replace site1 and site2 with the newly calculated designations
    data.table::setnames(LD, "site1", "site1_original")
    LD = merge(LD, sites, by.x = "site1_original", by.y = "site_original")
    data.table::setnames(LD, "site", "site1")

    data.table::setnames(LD, "site2", "site2_original")
    LD = merge(LD, sites, by.x = "site2_original", by.y = "site_original")
    data.table::setnames(LD, "site", "site2")

    # drop the original sites
    LD[,c("site1_original", "site2_original") := NULL]

    # reorder the columns
    data.table::setcolorder(
        LD,
        c(
            "locus",
            "position1",
            "site1",
            "position2",
            "site2",
            "distance",
            "r_squared"
        )
    )

    # sort the data by locus, position1, and distance
    data.table::setorder(LD, locus, position1, distance)

    # clean up the data before returning it
    LD[, locus := as.character(locus)]
    LD <- stats::na.omit(LD, cols = c("distance", "r_squared"))
    return(LD)
}
