#' Create plots for pathways based on user-provided filtering strategy.
#'
#' @param enrichment_data the data to be plotted (returned from
#'   find_pathway_significance())
#' @param significance_measure the parameter to be used for filtering
#' @param significance_cutoff the cut-off value of the significance parameter
#' @param output_directory a directory to save results in; will be created if it
#' @param theme set image to light or dark theme
#' doesn't exist
#' @import ggplot2
#' @export
#' @return no return
#' @examples
#' example("find_pathway_significance")
#' plot_enrichment_data(enrichment_data, "q-value", "0.05", tempdir())
plot_enrichment_data <- function(
    enrichment_data,
    significance_measure = "q-value",
    significance_cutoff = 1,
    output_directory,
    theme = "dark"
    analysis_mode = null,
) {

    # Assign some variables to NULL to pass checks.
    pathway_name <- enrichment_score <- gene_id <- pathway_id <- NULL

    # Check to see if the significance measure is set correctly.
    if (!(significance_measure %in% c("q-value", "p-value"))) {
        stop("Enrichment data can only be filtered on p-value or q-value.")
    }

    # Check to see if the theme is set correctly.
    if (!(theme %in% c("light", "dark"))) {
        stop("Theme must be either \"light\" or \"dark\".")
    }

    # Check to see if the output directory exists.
    # Attempt to create the directory if it doesn't exist.
    if (!dir.exists(file.path(output_directory))) {
        dir.create(file.path(output_directory))
    }

    # Plot every pathway in the filtered data.
    # A significance of cutoff of 1 plot get all pathways.
    # The plot's title is the pathway ID and the pathway name separated by a
    #   hyphen.
    # The y-intercept is the rank of the gene with the maximum enrichment score.
    # Save the plot to the output directory, naming the image using the pathway
    #   ID.
    enrichment_data[get(significance_measure) <= significance_cutoff, {
        title <- paste(.BY$pathway_id, .SD[1L, pathway_name], sep = " - ")
        intercept <- .SD[which.max(enrichment_score), rank]

        # Draw a line for running enrichment score.
        # Add marks for gene ranks.
        # Draw an intercept at the maximum running enrichment score.
        # Set the title.
        # Set the x-axis and y-axis names.
        # Set the theme to minimal.
        plot <-
            ggplot2::ggplot(.SD,
                            ggplot2::aes(
                                x = rank,
                                y = enrichment_score,
                                label = gene_id
                            )
            ) +
            ggplot2::geom_line(stat = "identity", color = "#fe4365") +
            ggplot2::geom_rug(sides = "t") +
            ggplot2::geom_vline(xintercept = intercept,
                                color = "red",
                                linetype = "longdash") +
            ggplot2::ggtitle(title) +
            ggplot2::labs(x = "Gene Rank", y = "Running Enrichment Score") +
            ggplot2::theme_minimal()

        if (!is.null(analysis_mode)) {
            plot <- plot + ggplot2::geom_point(
                size = 2,
                color = "#A7011F",
                fill = "#fe4365",
                ggplot2::aes(shape = gene_effect <= 0)
            ) +
                ggplot2::scale_shape_manual(
                    values = c(24, 25),
                    guide = "none"
                )
        }

        if (theme == "dark") {
            # Modify elements of the theme to create a dark theme.
            plot <- plot +
                ggplot2::geom_rug(sides = "t", color = "#FFFFFF") +
                ggplot2::theme(
                    axis.text = ggplot2::element_text(color = "#a1a1a1"),
                    plot.background = ggplot2::element_rect(
                        color = "#a1a1a1",
                        fill = "#222831"
                    ),
                    panel.background = ggplot2::element_rect(
                        color = "#a1a1a1",
                        fill = "#222831"
                    ),
                    plot.title = ggplot2::element_text(color = "#a1a1a1"),
                    axis.title = ggplot2::element_text(color = "#a1a1a1"),
                    panel.grid.minor = ggplot2::element_blank()
                )
        }

        ggplot2::ggsave(file.path(
            output_directory,
            paste(.BY$pathway_id, "png", sep = ".")),
            plot)
        NULL
    },
    by = pathway_id]
    return(NULL)
}
