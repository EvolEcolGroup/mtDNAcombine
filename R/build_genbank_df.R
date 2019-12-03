#' Builds a reference data frame with details of all available genes associated
#' with each accession investigated
#'
#' \code{build_genbank_df}
#'
#' This function extracts and collates details of all available genes (as well
#' as the species scientific name) for each accession given.
#'
#' Function designed to construct a raw reference data frame of all sequence
#' information available for each of the accession numbers given. Uses functions
#' \code{load_accession_list}, \code{get_sample_details}, and
#' \code{get_gene_details} to extract information from NCBI website. Output will
#' need downstream processing to remove duplicated pieces of information and to
#' clarify any nomlencature issues for example.
#'
#' @param accession_file_name user provided file name which lists accession
#' numbers. Should be a .csv or .txt file without a header.
#'
#' @export

build_genbank_df <- function(accession_file_name) {
    if (inherits(accession_file_name, "character") == FALSE) {
        stop("accession_file_name is not a file name character string")
    }

    accession_list <- load_accession_list(accession_file_name)

    genbank_download <- as.data.frame(NULL)
    genbank_download2 <- as.data.frame(NULL)
    # Errors <-as.data.frame(NULL) Errors2 <- as.data.frame(NULL)


    for (i in 1:nrow(accession_list)) {
        accession_no <- as.character(accession_list[i, ])
        error_catch <- tryCatch({
            connection <- curl::curl(
                paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                       accession_no, "&rettype=gb&retmode=xml"))
            full_xmlTree <- readLines(connection)
            # close the connection
            close(connection)

            full_xmlTree <- XML::xmlInternalTreeParse(full_xmlTree, asText = T)


            sample_details <- get_sample_details(full_xmlTree = full_xmlTree)

            gene_details <-
                get_gene_details(feature_name = sample_details$feature_name,
                    freq = sample_details$freq, full_xmlTree = full_xmlTree)


            genbank_download2 <-
                as.data.frame(cbind(
                        sample_details$sci_nam,
                        as.character(gene_details$gene_nam),
                        as.numeric(as.character(gene_details$position_start)),
                        as.numeric(as.character(gene_details$position_end)),
                        sample_details$accession_version,
                        sample_details$create_date,
                        sample_details$download_date))

            genbank_download <-
                as.data.frame(rbind(genbank_download, genbank_download2))


        },
        error = function(e) {
            cat("ERROR", accession_no, ":", conditionMessage(e), "\n")
            return(NA)
        })


    }


    colnames(genbank_download) <-
        c("sci_nam", "gene_nam", "position_start", "position_end",
            "accession_version", "create_date", "download_date")

    return(genbank_download)

}
