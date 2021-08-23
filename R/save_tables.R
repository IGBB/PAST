#' Write tables for pathways based on user-provided filtering strategy.
#'
#' @param enrichment_data the data to be plotted (returned from
#'   find_pathway_significance())
#' @param significance_measure the parameter to be used for filtering
#' @param significance_cutoff the cut-off value of the significance parameter
#' @param output_directory a directory to save results in; will be created if it
#' doesn't exist
#' @import data.table
#' @export
#' @return no return
#' @examples
#' example("find_pathway_significance")
#' write_enrichment_data(enrichment_data, "q-value", "0.05", tempdir())
#' write_enrichment_data(enrichment_data, output_directory=tempdir())
write_enrichment_data <- function(
    enrichment_data,
    significance_measure = "q-value",
    significance_cutoff = 1,
    output_directory
) {

    # Check to see if the significance measure is set correctly.
    if (!(significance_measure %in% c("q-value", "p-value"))) {
        stop("Enrichment data can only be filtered on p-value or q-value.")
    }

    # Check to see if the output directory exists.
    # Attempt to create the directory if it doesn't exist.
    if (!dir.exists(file.path(output_directory))) {
        dir.create(file.path(output_directory))
    }

    # If the significance cutoff value is 1, then don't indicate filtering
    #   information in the filename.
    # Otherwise, the output file name is set to indicate how the enrichment data
    #   was filtered.
    if (significance_cutoff == 1) {
        output_path <- file.path(output_directory, "enrichment_scores.tsv")
    } else {
        output_path <- file.path(
            output_directory,
            paste(
                "enrichment_scores",
                "filtered_by",
                significance_measure,
                "at",
                gsub("\\.", "_", significance_cutoff),
                "tsv",
                sep = "."
            )
        )
    }

    # Write the data to the file, filtering as necessary.
    # A significance of cutoff of 1 will write all data.
    data.table::fwrite(
        enrichment_data[get(significance_measure) <= significance_cutoff],
        file = output_path
    )
}
