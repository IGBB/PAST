#' Plot Rugplots for Selected Pathways
#'
#' @param rugplots_data The data to be plotted (returned from
#'   find_pathway_significance())
#' @param filter_type The parameter to be used for filtering
#' @param filter_parameter The cut-off value of the filtering parameter
#' @param mode The mode used to create the data (increasing/decreasing)
#' @param output_directory An existing directory to save results in
#' @importFrom rlang .data
#' @import dplyr
#' @import ggplot2
#' @export
#' @return Does not return a value
#' @examples
#' example("find_pathway_significance")
#' plot_pathways(rugplots_data, "pvalue", "0.03", "decreasing", tempdir())
plot_pathways <-
  function(rugplots_data,
           filter_type,
           filter_parameter,
           mode,
           output_directory) {
    
    rugplots_data <- rugplots_data %>%
      dplyr::arrange(.data$pathway_number)
    write.table(
      rugplots_data,
      file = paste0(output_directory, "/", mode, ".full.txt"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    if (filter_type == "rank") {
      rugplots_data <- rugplots_data %>%
        dplyr::filter(.data$pathway_number <= filter_parameter)
    } else if (filter_type == "pvalue") {
      rugplots_data <- rugplots_data %>%
        dplyr::filter(.data$pvalue <= filter_parameter)
    } else if (filter_type == "qvalue") {
      rugplots_data <- rugplots_data %>%
        dplyr::filter(.data$qvalue <= filter_parameter)
    } else {
      print("Incorrect filtering type. Filtering at p-value <= 0.05")
      rugplots_data <- rugplots_data %>%
        dplyr::filter(.data$pvalue <= 0.05)
    }
    
    write.table(
      rugplots_data,
      file = paste0(output_directory, "/", mode, ".filtered.txt"),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    rugplots_split <- split(rugplots_data, rugplots_data$pathway_number)
    
    for (rank in names(rugplots_split)) {
      temp_data <- rugplots_split[[rank]]
      title <-
        paste0(unique(as.character(temp_data$pathway_id)), " - ",
               unique(as.character(temp_data$pathway_name)))
      intercept <- temp_data %>%
        dplyr::arrange(desc(.data$running_enrichment_score)) %>%
        dplyr::select(.data$rank)
      intercept <- intercept[, 1][1]
      rugplot <-
        ggplot(temp_data,
               aes(x = temp_data$rank,
                   y = temp_data$running_enrichment_score)) +
        geom_line(stat = "identity") +
        geom_rug(sides = "t", position = "jitter") +
        geom_vline(xintercept = intercept,
                   color = "black",
                   linetype = "longdash") +
        ggtitle(title) +
        labs(x = "Gene Rank", y = "Running Enrichment Score") +
        scale_x_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000)) +
        theme(
          axis.text = element_text (color = "black"),
          panel.background = element_rect (color = "black", fill = "pink")
        )
      ggsave(paste0(
        output_directory,
        "/",
        mode,
        ".",
        unique(as.character(temp_data$pathway_id)),
        ".png"
      ),
      rugplot)
    }
  }
