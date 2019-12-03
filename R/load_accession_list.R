#' Loads list of accession numbers
#'
#' \code{load_accession_list}
#'
#' This function reads in list of accession numbers provided by the user and
#' cleans up formatting (such as differences between "accession number" and
#' "accession version" (e.g. xxxx vs xxxx.1)).
#' Also any accession numbers starting with `[XX]_` are reviewed. These are RefSeq
#' accessions and are known to frequently be duplicates of other accessions in
#' GenBank. If the "original" accession is in the accession list already the
#' `[XX]_` accession number will be removed.
#'
#' @param accession_file_name user provided file which lists accession numbers.
#' Should be a .csv or .txt file. No header.
#'
#' @export


load_accession_list <- function(accession_file_name) {
    # check type of file being read in
    what <- stringr::str_sub(accession_file_name, start = -3)
    if (what == "csv") {
        accession_list <- utils::read.csv(accession_file_name, header = F)
    } else {
        if (what == "txt") {
            accession_list <- utils::read.table(accession_file_name, header = F)
        }
    }
    # check file format
    if (ncol(accession_list) > 1) {
        stop("This file should be a single column of accession numbers")
    }

    # be sure the accession numbers are formated the same way
    accession_list <- as.data.frame(apply(accession_list, 2, function(x)
                        gsub("\\.[[:digit:]]", "", x)))

    rm_accession_list <- NULL
    for (p in 1:nrow(accession_list)) {
        # accessions with underscores are from RefSeq, often repeates, find them
        if (stringr::str_detect(accession_list[p, 1], "_")) {
            accession_no <- as.character(accession_list[p, ])
            connection <- curl::curl(paste0(
                        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                        accession_no, "&rettype=gb&retmode=xml"))
            full_xmlTree <- readLines(connection)
            # close the connection
            close(connection)

            full_xmlTree <- XML::xmlInternalTreeParse(full_xmlTree, asText = T)
            # comments for accessions starting 'NC_' contain information on
            # duplicate accessions
            comments <- XML::getNodeSet(full_xmlTree, "//GBSeq_comment")[[1]]

            if (!is.null(comments)) {
                comments <- XML::getNodeSet(full_xmlTree, "//GBSeq_comment")[[1]]
                comments <- XML::xmlValue(comments)
                # stock phrase 'The reference sequence is identical to '
                # used when duplicate exists - find
                orig_accession <-
                    stringr::word(string = sapply(
                        strsplit(
                            comments,
                            "The reference sequence is identical to "
                        ),
                        "[",
                        2
                    ),
                    start = 1)
                # following word gives original accession number
                orig_accession <- gsub("[[:punct:]]", "", orig_accession)

                # if the orginal accession number just found is in the accession
                # list being useds then this sample has already been captured
                # and this '[XX]_' duplicate should be removed
                if (orig_accession %in% accession_list$V1) {
                    rm_accession_list <-
                        c(rm_accession_list,
                          as.character(accession_list[p, 1]))
                }
            }
        }
    }

    # drop those 'NC_' accessions numbers that are duplicates
    rm_accession_list2 <- as.data.frame(rm_accession_list)
    accession_list2 <- as.data.frame(accession_list[!accession_list$V1 %in%
                                        rm_accession_list2$rm_accession_list, ])
    accession_list <- as.matrix(accession_list2)
    colnames(accession_list) <- "accessions"

    return(accession_list)


}
