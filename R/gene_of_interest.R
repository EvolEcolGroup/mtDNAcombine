#' Extracts data frame rows associated with a specified gene of interest
#'
#' \code{gene_of_interest}
#'
#' Function takes data frame, likely built from all available data associated
#' with the original input accession numbers, and returns only data relevant to
#' specified gene of interest.
#'
#' This function filters all data collected and curated so far from the original
#' accession list down to just information associated with a nominated gene of
#' interest.
#'
#' @param gene gene of interest written as character string. Must be in
#' "standardised" format used to simplify data in \code{standardise_gene_names}
#' function.
#' @param data data frame (expected to be `GB_data`, the data frame result of
#' the pipeline so far).
#'
#' @export


gene_of_interest <- function(gene, data) {
    if (inherits(data, "data.frame") == FALSE) {
        # check type of file being read in
        stop("df_to_update is not a data.frame object")
    }
    what <- length(colnames(data))
    if (what == 7 | 8) {
        GB_by_gene <- data[data$gene_name == gene, ]
        return(GB_by_gene)
    } else {
        stop("data not built correctly")
    }
}
