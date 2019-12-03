#' Removes multiple copies of one gene
#'
#' \code{remove_duplicates}
#'
#' Any instances of one gene being recorded multiple times for a single accession
#' number are found and reduced to a unique occourance.
#'
#' Due to the structure of NCBI xml tree, after \code{standardise_gene_names} has
#' been used, there are likely to be multiple occourances of one gene associated
#' with a single accession number, these need to be removed.
#'
#' @param df_to_update data frame with a column name `$gene_name`. Expected to be
#' `GB_data`, the result of pipeline so far.
#'
#' @export

remove_duplicates <- function(df_to_update) {
    if (inherits(df_to_update, "data.frame")) {
        df_to_update$merged <- paste(df_to_update$gene_name,
                                     df_to_update$accession_version, sep = "__")
        df_to_update <- df_to_update[match(unique(df_to_update$merged),
                                           df_to_update$merged), ]
        df_to_update <- df_to_update[, !(names(df_to_update) == "merged")]
        return(df_to_update)
    } else {
        stop("df_to_update is not a data.frame object")
    }

}
