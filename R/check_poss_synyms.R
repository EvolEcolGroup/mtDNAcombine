#' Checks if scientific name for sample is likely a species, or subspecies name
#'
#' \code{check_poss_synyms}
#'
#' Function takes data frame with column $sci_nam and returns a list of
#' scientific names that are longer than 2 words. Also writes out a .csv file of
#' these names.
#'
#' Sequence data can be uploaded to GenBank with scientific names that may be
#' non-standard/of different taxonomic rank, filtering based on species name will
#' therefore not work unless names are standardised. This function takes any
#' instances of names longer than 2 words, outputs the result as an object and a
#' .csv which can then be edited by the user.
#'
#' @param data data frame with `$sci_nam` column (expected to be `GB_with_SeqDat`,
#' the data frame result of the pipeline so far).
#' @param file_path path to newly created .csv file, defaults to working directory
#'
#' @export

check_poss_synyms <- function(data, file_path = getwd()) {
    poss_synyms <- NULL
    
    if ("sci_nam" %in% colnames(data)) {
        data$sci_nam <- as.character(data$sci_nam)
        for (b in 1:nrow(data)) {
            if (sapply(strsplit(data$sci_nam[b], " "), length) > 2) {
                poss_synyms <- rbind(poss_synyms, data$sci_nam[b])
            }
        }
        poss_synyms <- unique(poss_synyms)

        if (length(poss_synyms) > 0) {
            utils::write.csv(x = poss_synyms, 
                             file = paste0(file_path,"/poss_synyms.csv"), 
                             quote = FALSE, row.names = FALSE)
        }
        return(poss_synyms)
    } else {
        stop("no sci_nam column in data")
    }
}
