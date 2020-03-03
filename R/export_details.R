#' Exports downloaded details into seperate .csv files for each species
#'
#' \code{export_details}
#'
#' Function takes data frame with column $sci_nam and splits according to species
#' name before exporting sample details as individual .csv files.
#'
#' Main input data frame is split by species (according to $sci_nam column). For
#' each species a separate .csv file is then written out containing all the
#' information downloaded and stored about each accession.
#'
#' @param data data frame with $sci_nam column (expected to be GB_with_SeqDat,
#' the data frame result of the pipeline so far).
#' @param file_path path to newly created .csv file, defaults to working directory
#'
#' @export



export_details <- function(data, file_path = getwd()) {
    if ("sci_nam" %in% colnames(data)) {
        exporting <- split(data, data$sci_nam, drop = TRUE)

        for (n in 1:length(exporting)) {
            utils::write.csv(exporting[[n]],
                             file = paste0(file_path, "/", 
                                           exporting[[n]]$sci_nam[1], ".csv"),
                row.names = FALSE, quote = FALSE)
        }
    } else {
        stop("data doesn't have a sci_nam column")
    }
}
