#' Updates any species names user wants to alter
#'
#' \code{standardise_spp_names}
#'
#' Function takes data frame with column $sci_nam and uses a .csv file containing
#' existing names and desired replacement names to standise species nomenclature
#' across all accessions.
#'
#' Sequence data can be uploaded to GenBank with scientific names that may be
#' non-standard/of different taxonomic rank, filtering based on species name will
#' therefore not work unless names are standardised. This function reads in a .csv
#' file of old and user-defined replacement species names, these are then used to
#' update the data$sci_nam column.
#'
#' @param data data frame with `$sci_nam` column (expected to be `GB_with_SeqDat`,
#' the data frame result of the pipeline so far).
#' @param new_names_file A csv file name. csv file with two columns. First column
#' is existing scientific name, second column is name to replace existing
#' character string with. Structure built on exported "poss_synonyms.csv" file.
#'
#' @export


standardise_spp_names <- function(data, new_names_file) {
    if (inherits(new_names_file, "character")) {
        new_names <- utils::read.csv(new_names_file, header = T, stringsAsFactors = T)
        what <- ncol(new_names)

        if (what == 2) {
            for (n in 1:nrow(new_names)) {
                if (new_names$V1[n] %in% data$sci_nam) {
                    data$sci_nam <- gsub(new_names$V1[n], new_names$X[n],
                                              data$sci_nam)
                }
            }
            return(data)
        } else {
            stop("new_names_file not built correctly")
        }

    } else {
        stop("new_names_file is not character string file name")
    }
}
