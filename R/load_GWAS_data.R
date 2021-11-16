#' Load GWAS data
#'
#' @param first_file a file containing combined GWAS data or containing
#' association data
#' @param second_file a file containing effects data; only use when first_file
#' contains association data
#' @param input_type one of "single", "two", or "TASSEL"
#' @param mode homozygous or heterozygous; retains either only homozygous data
#' or homozygous and heterozygous data
#' @param trait the name of the column containing the trait
#' @param marker the name of the column containing the marker
#' @param locus the name of the column containing the locus or chromosome
#' @param site the name of the column containing the site or position
#' @param p the name of the column containing the p-value
#' @param effect the name of the column containing the effect
#' @param effects_trait the name of the column containing the trait in the
#' effects file; only used when input_type is "two" or "TASSEL"
#' @param effects_marker the name of the column containing the marker in the
#' effects file; only used when input_type is "two" or "TASSEL"
#' @param effects_locus the name of the column containing the locus or
#' chromosome in the effects file; only used when input_type is "two" or
#' "TASSEL"
#' @param effects_site the name of the column containing the site or position
#' in the effects file; only used when input_type is "two" or "TASSEL"
#' @param effects_allele the name of the column containing the allele in the
#' effects file; only used when input_type is "TASSEL" and mode is
#' "increasing" or "decreasing"
#' @return GWAS data with a marker created from locus/chromosome and
#' site/position, the original marker, locus/chromosome, site/position, p-value,
#' marker R^2, and effect for every marker in the input files
#' @export
#' @import data.table
#' @import glue
#' @examples
#' association <- system.file("extdata", "association.txt.gz",
#'     package = "PAST2",
#'     mustWork = TRUE
#' )
#' effects <- system.file("extdata", "effects.txt.gz",
#'     package = "PAST2",
#'     mustWork = TRUE
#' )
#' single <- system.file("extdata", "single_gwas.txt.gz",
#'     package = "PAST2",
#'     mustWork = TRUE
#' )
#' effects_single <- system.file("extdata", "effects.single_line.txt.gz",
#'     package = "PAST2", mustWork = TRUE
#' )
#' gwas_data <- load_GWAS_data(
#'     single,
#'     NULL,
#'     "single",
#'     "homozygous",
#'     trait = "Trait",
#'     marker = "Marker",
#'     locus = "Locus",
#'     site = "Site",
#'     p = "p",
#'     effect = "Effect"
#' )
#' gwas_data <- load_GWAS_data(
#'     association,
#'     effects_single,
#'     "two",
#'     "homozygous",
#'     trait = "Trait",
#'     marker = "Marker",
#'     locus = "Locus",
#'     site = "Site",
#'     p = "p",
#'     effect = "Effect"
#' )
#' gwas_data <- load_GWAS_data(
#'     association,
#'     effects,
#'     "TASSEL",
#'     "homozygous",
#'     trait = "Trait",
#'     marker = "Marker",
#'     locus = "Locus",
#'     site = "Site",
#'     p = "p",
#'     effect = "Effect"
#' )
load_GWAS_data <- function(
    first_file,
    second_file = NULL,
    input_type,
    mode,
    trait,
    marker,
    locus,
    site,
    p,
    effect,
    effects_trait = NULL,
    effects_marker = NULL,
    effects_locus = NULL,
    effects_site = NULL,
    effects_allele = NULL
) {

    # Set some variables to NULL for BiocCheck.
    . <- marker_original <- allele <- num_occurrences <- NULL
    chromosome <- position <- NULL

    columns <- list(
        "trait" = trait,
        "marker" = marker,
        "locus" = locus,
        "site" = site,
        "p" = p,
        "effect" = effect,
        "effects_trait" = effects_trait,
        "effects_marker" = effects_marker,
        "effects_locus" = effects_locus,
        "effects_site" = effects_site,
        "effects_allele" = effects_allele
    )
    columns <- columns[!sapply(columns,is.null)]

    if (!(is.null(second_file))) {
        files <- c(first_file, second_file)
    } else {
        files <- c(first_file)
    }

    # Check for invalid settings.
    if (!(input_type) %in% c("single", "two", "TASSEL")) {
        stop(
            glue::glue(
                "Input type must be one of ",
                "\"single\", \"two\", or \"TASSEL\", ",
                "not {input_type}"
            )
        )
    }

    if (!(mode) %in% c("homozygous", "increasing", "decreasing")) {
        stop(
            glue::glue(
                "Mode must be one of ",
                "\"homozygous\", \"increasing\", or \"decreasing\", ",
                "not {mode}"
            )
        )
    }

    # Check for incompatible input_types and number of files.
    if (
        (length(files) == 1 & input_type != "single") |
        (length(files) == 2 & input_type == "single")) {
        file_word <- if (length(files) == 1) "file" else "files"
        stop(
            glue::glue(
                "Incompatible number of files and input_type: ",
                "{length(files)} {file_word} provided but ",
                "input_type is '{input_type}'"
            )
        )
    } else if (length(files) > 2) {
        stop(
            glue::glue(
                "PAST accepts 1 or 2 files as input. ",
                "You provided {length(files)} files."
            )
        )
    }

    # Check for incompatible input_types and mode.
    if (input_type != "TASSEL" & mode == "heterozygous") {
        stop(
            glue::glue(
                "Incompatible input_type: ",
                "input_type '{input_type}' is not compatible with ",
                "heterozygous mode"
            )
        )
    }

    # Set a list of required columns to check for.
    required_columns <- c(
        'trait',
        'marker',
        'locus',
        'site',
        'p',
        'effect'
    )

    if (mode %in% c("increasing", "decreasing")) {
        required_columns <- c(required_columns, "effects_allele")
    }

    # Check for missing required columns.
    missing_columns <- setdiff(required_columns, names(columns))
    if (length(missing_columns > 0)) {
        stop(
            glue::glue(
                "The following colums are required but not set: ",
                "{paste(missing_columns, collapse=', ')}"
            )
        )
    }

    if (input_type == "single") {

        # Set columns for single-file data.
        # Rename the named parameters.
        # Load data.
        gwas_columns <- columns[names(columns) %in% required_columns]
        names(gwas_columns)[[match(
            "marker",
            names(gwas_columns)
        )]] <- "marker_original"
        names(gwas_columns)[[match(
            "locus",
            names(gwas_columns)
        )]] <- "chromosome"
        names(gwas_columns)[[match(
            "site",
            names(gwas_columns)
        )]] <- "position"
        names(gwas_columns)[[match(
            "p",
            names(gwas_columns)
        )]] <- "p-value"
        gwas_data <- load_file(files[[1]], gwas_columns)

    } else if (input_type == "two" | input_type == "TASSEL") {

        # Set the columns for the association data and rename some of them.
        association_columns <- columns[names(columns) %in% required_columns]
        names(association_columns)[[match(
            "marker", names(association_columns)
        )]] <- "marker_original"
        names(association_columns)[[match(
            "locus", names(association_columns)
        )]] <- "chromosome"
        names(association_columns)[[match(
            "site", names(association_columns)
        )]] <- "position"
        names(association_columns)[[match(
            "p", names(association_columns)
        )]] <- "p-value"
        association_columns[["effect"]] <- NULL

        # Set the columns for the effects data.
        # Rename effects_allele if mode is heterozygous.
        effects_columns <- list("effect" = columns[['effect']])
        if (mode %in% c("increasing", "decreasing")) {
            effects_columns$effects_allele <- columns[['effects_allele']]
            names(effects_columns)[[match(
                "effects_allele",
                names(effects_columns)
            )]] <- "allele"
        }

        # Set a list of optional columns.
        extra_columns <- list(
            'effects_trait' = 'trait',
            'effects_marker' = 'marker',
            'effects_locus' = 'locus',
            'effects_site' = 'site'
        )

        # Determine which optional columns were provided.
        provided_extras <- names(extra_columns)[
            names(extra_columns) %in% names(columns)
        ]

        not_provided_extras <- unname(extra_columns[
            !names(extra_columns) %in% names(columns)
        ])

        # Build a list of effects columns.
        effects_columns[
            unlist(c(provided_extras, not_provided_extras))
        ] <- lapply(
            unlist(c(provided_extras, not_provided_extras)),
            function(column) {
                if (column %in% names(columns)) {
                    return(columns[[column]])
                } else {
                    return(columns[[extra_columns[[column]]]])
                }
            }
        )

        # Rename the marker column.
        if ("marker" %in% names(effects_columns)) {
            names(effects_columns)[[match(
                "marker",
                names(effects_columns)
            )]] <- "marker_original"
        } else if ("effects_marker" %in% names(effects_columns)) {
            names(effects_columns)[[match(
                "effects_marker",
                names(effects_columns)
            )]] <- "marker_original"
        }

        # Load the data.
        association_data <- load_file(files[[1]], association_columns)
        effects_data <- load_file(files[[2]], effects_columns)

        # Drop most of the columns.
        if (mode == "homozygous") {
            effects_data <- effects_data[, .(effect, marker_original)]
        } else {
            effects_data <- effects_data[, .(effect, marker_original, allele)]
        }

        # Check to make sure that shape of data matches input type.
        # Offer the option to continue the analysis as the correct mode.
        if (input_type == "two" &
            nrow(association_data) != nrow(effects_data)) {
            message("Effects file contains more than one effect per marker.")
            continue_as_TASSEL <- readline(
                prompt = glue::glue(
                    "Continue as TASSEL {mode} analysis? (Y/N) "
                )
            )
            while (!(tolower(continue_as_TASSEL) %in% c("y", "n"))) {
                message("Invalid choice.")
                continue_as_TASSEL <- readline(
                    prompt = glue::glue(
                        "Continue as TASSEL {mode} analysis? (Y/N) "
                    )
                )
            }
            if (tolower(continue_as_TASSEL) == "n") {
                stop("Analysis stopped. Check your data and try again.")
            } else {
                message("Changing input_type to \"TASSEL\" and continuing")
                input_type <- "TASSEL"
            }
        } else if (
            input_type == "TASSEL" &
            nrow(association_data) == nrow(effects_data)
        ) {
            message("Effects file contains only one effect per marker.")
            continue_as_two <- readline(
                prompt = glue::glue(
                    "Continue as two-file {mode} analysis? (Y/N) "
                )
            )
            while (!(tolower(continue_as_two) %in% c("y", "n"))) {
                message("Invalid choice.")
                continue_as_two <- readline(
                    prompt = glue::glue(
                        "Continue as TASSEL {mode} analysis? (Y/N) "
                    )
                )
            }
            if (tolower(continue_as_two) == "n") {
                stop("Analysis stopped. Check your data and try again.")
            } else {
                message("Changing input_type to \"two\" and continuing")
                input_type <- "two"
            }
        }

        if (input_type == "TASSEL" & mode == "homozygous") {
            # Check that TASSEL data is correctly formatted (every other
            #   effect == 0) and take rows with non-zero effects.
            message("homozygous mode: retain only biallelic homozygous data")
            effects_data[, row := .I]
            if (sum(effects_data[row %% 2 == 0, effect]) == 0) {
                effects_data <- effects_data[row %% 2 == 1]
                effects_data[,row := NULL]
            } else {
                stop(
                    glue::glue(
                        "PAST expected even rows in the effects file ",
                        "{files[[2]]} to have effect == 0."
                    )
                )
            }
        } else if (
            input_type == "TASSEL" &
            (mode == "increasing" | mode == "decreasing")
        ) {
            # Set heterozygous alleles.
            heterozygous <- c("AGR", "CTY", "CGS", "ATW", "GTK", "ACM")

            # Set up row choice based on mode.
            if (mode == "decreasing") {
                message("decreasing mode: ",
                        "retaining marker with largest negative effect"
                )
                one <- 1
                row_choice <- "one"
            } else if (mode == "increasing") {
                message("increasing mode: ",
                        "retaining marker with largest positive effect"
                )
                row_choice <- ".N"
            }

            # Count number of times each marker appears.
            effects_data[, num_occurrences := .(.N), by = .(marker_original)]

            # Get homozygous data.
            homozygous_data <- effects_data[num_occurrences == 2]
            homozygous_data[, row := .I]

            if (sum(homozygous_data[row %% 2 == 0, effect]) == 0) {
                homozygous_data <- homozygous_data[row %% 2 == 1]
                homozygous_data[, row := NULL]
            } else {
                stop(
                    glue::glue(
                        "PAST expected even rows in the effects file ",
                        "{files[[2]}} to have effect == 0."
                    )
                )
            }
            data.table::setcolorder(
                homozygous_data,
                c(
                    "marker_original",
                    "effect",
                    "allele",
                    "num_occurrences"
                )
            )

            # Find markers of heterozygous/triallelic data.
            number_of_markers <- data.table::uniqueN(
                effects_data[num_occurrences == 3]$marker
            )

            # Set up progress bar.
            progress <- utils::txtProgressBar(
                min = 0,
                max = number_of_markers,
                style = 3
            )

            # Set effects_data equal to homozygous data + heterozygous data.
            effects_data <- data.table::rbindlist(
                list(
                    homozygous_data,
                    effects_data[num_occurrences == 3, {
                        # Get alleles and sort them.
                        alleles <- sort(
                            paste(
                                .SD[, allele],
                                collapse = ""
                            )
                        )
                        # Mark progress.
                        utils::setTxtProgressBar(progress, .GRP)
                        # If heterozygous and not tri-allelic, take
                        #   appropriate effect.
                        if (alleles %in% heterozygous) {
                            .SD[effect != 0][get(row_choice)]
                        } else {
                            NULL
                        }
                    },
                    by = marker_original
                    ]
                )
            )
        }

        # Merge association and effects data.
        gwas_data <- association_data[effects_data, on = .(marker_original)]
    }

    # Add an internal marker and drop trait.
    gwas_data[
        ,
        marker := paste(chromosome, position, sep = "_")
    ][
        ,
        trait := NULL
    ]

    # Set the order of the columns.
    data.table::setcolorder(
        gwas_data,
        c(
            "marker",
            "chromosome",
            "position",
            "p-value",
            "effect",
            "marker_original"
        )
    )

    # Convert chromosome to character if it isn't already.
    gwas_data[, chromosome := as.character(chromosome)]

    # Key the data by chromosome/position order.
    data.table::setkeyv(gwas_data, c("chromosome", "position"))

    return(gwas_data)
}

