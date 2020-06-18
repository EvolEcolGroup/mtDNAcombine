#' Standardises the range of gene names used in GenBank
#'
#' \code{standardise_gene_names}
#'
#' This function turns an extensive but not exhaustive list of possible gene
#' name synonyms/miss-spellings/alternative abbreviations into standardised format.
#'
#' Function loads either user provided .csv file or default .csv file
#' ("gene_nomenclature_correction.csv") and uses gsub to replace all alternative
#' gene name patterns given with a single, standard, gene name.
#'
#' @param df_to_update data frame with a column name $gene_name. Any instances of
#' multiple names for what is a single gene within this column will be standardised
#' @param names_to_replace user provided .csv file or default .csv file
#' ("gene_nomenclature_correction.csv") which lists alternative names for each
#' gene. One gene for each column, header is new standardised gene name option,
#' rest of column is filled with erroneous versions of the gene name.
#'
#' @export


standardise_gene_names <- function(df_to_update, names_to_replace) {
    default <- missing(names_to_replace)
    tmp <- NULL

    if (default == TRUE) {
        # use default list
        filepath <-
            system.file("extdata", "gene_nomenclature_correction.csv",
                        package = "mtDNAcombine")
        names_to_replace <- utils::read.csv(file = filepath, header = T,
                                            stringsAsFactors = T)
    } else {
        if (inherits(names_to_replace, "character")) {
            names_to_replace <- utils::read.csv(names_to_replace, header = T,
                                                stringsAsFactors = T)
        } else {
            stop("file name for names_to_replace not a character string")
        }

    }

    for (f in 1:ncol(names_to_replace)) {
        tmp2 <- paste0(names_to_replace[, f], collapse = "", sep = "|")
        # drop | on the end of string
        while (stringr::str_sub(tmp2, start = -1) == "|") {
            tmp2 <- substr(tmp2, 1, nchar(tmp2) - 1)
        }
        tmp <- rbind(tmp, tmp2)
    }
    rows <- colnames(names_to_replace)
    rownames(tmp) <- rows


    # replacement step
    for (p in 1:length(tmp)) {
        df_to_update$gene_name <- gsub(df_to_update$gene_name, pattern = tmp[p],
                replacement = rownames(tmp)[p])
    }
    return(df_to_update)
}
