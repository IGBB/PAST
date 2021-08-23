#' Load GWAS data
#'
#' @param file the file to read
#' @param columns a named list of the columns describing what they are and the
#' user-provided name
#' @param file_has_header a boolean indicating whether the file has a header
#' or not
#' @return file with columns specified by user named according to names of
#' columns list
#' @import data.table
#' @import glue
load_file <- function(file, columns, file_has_header=TRUE){

    # Check that the file exists.
    if (!file.exists(file)) {
        stop(glue::glue("{file} does not exist"))
    }

    # Check the headers of the data to be sure that all user-requested
    #   columns exist.
    # If they don't, throw an error and exit.
    if (file_has_header) {
        header <- unlist(data.table::fread(file, nrows = 1, header = FALSE))
        column_check <- setdiff(columns, header)
        if (length(column_check) > 0) {
            stop(
                glue::glue(
                    "The following columns were not found in {file}:\n\t",
                    "{paste(column_check, collapse=', ')}.\n",
                    "{file} has the following columns:\n\t",
                    "{paste(header, collapse=', ')}."
                )
            )
        }

        # Read the file and select its columns.
        # Rename the columns using the names of columns.
        data <- data.table::fread(
            file,
            select = unlist(unname(columns))
        )

        data.table::setnames(
            data,
            unlist(unname(columns)),
            unlist(names(columns))
        )

        return(data)
    } else {
        message("Not using header")
    }
}
