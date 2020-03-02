#' Builds a curated data frame with sample details and sequence data for a
#' specified gene given accession numbers.
#'
#' \code{get_GB_sequence_data}
#'
#' This function extracts and collates details of the sample, as well as the
#' sequence data for a selected gene, from GenBank for each accession number given.
#'
#' Function designed to construct a comprehensive reference data frame of sample
#' and sequence data available for the gene of interest from each of the
#' accession numbers given. Uses functions \code{load_accession_list},
#' \code{get_sample_details}, and \code{get_gene_details} to extract information
#' from NCBI website. As well as \code{standardise_gene_names} and
#' \code{remove_duplicates} to clean up data frame.
#'
#' @param accessions_of_interest user provided file which lists accession numbers.
#' Should be a .csv or .txt file without a header.
#' @param gene gene of interest written as character string.
#' @param new_names_file csv file of species names that have been flagged as
#' problematic and the replacement name to be used, needed for
#' `standardise_spp_name` function
#'
#' @export




get_GB_sequence_data <-
    function(accessions_of_interest, gene, new_names_file) {
        # check type of file being read in
        if (inherits(accessions_of_interest, "data.frame")) {
            checking <- length(colnames(accessions_of_interest))
            if (checking == 7) {
                accession_list <-
                    as.data.frame(accessions_of_interest[, "accession_version"])
            } else {
                stop("data frame not structured as expected")
            }




            GB_data_details <- as.data.frame(NULL)
            GB_data2_details <- as.data.frame(NULL)


            for (i in 1:nrow(accession_list)) {
                accession_no <- as.character(accession_list[i, ])
                error_catch <- tryCatch({
                    connection <- curl::curl(
                        paste0(
                            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                                accession_no, "&rettype=gb&retmode=xml"))
                    full_xmlTree <- readLines(connection)
                    # close the connection
                    close(connection)

                    full_xmlTree <-
                        XML::xmlInternalTreeParse(full_xmlTree, asText = T)


                    sample_details <-
                        get_sample_details(full_xmlTree = full_xmlTree)

                    gene_details <-
                        get_gene_details(feature_name = sample_details$feature_name,
                            freq = sample_details$freq,
                            full_xmlTree = full_xmlTree)


                    GB_data2_details <-
                        as.data.frame(cbind(sample_details$sci_nam,
                                as.character(gene_details$gene_name),
                                as.numeric(as.character(
                                    gene_details$position_start)),
                                as.numeric(as.character(
                                    gene_details$position_end)),
                                as.character(gene_details$the_sequence),
                                sample_details$accession_version,
                                sample_details$create_date,
                                sample_details$download_date))
                    GB_data_details <-
                        as.data.frame(rbind(GB_data_details, GB_data2_details))



                },
                error = function(e) {
                    cat("ERROR", accession_no, ":", conditionMessage(e), "\n")
                    return(NA)
                })

            }

            colnames(GB_data_details) <-
                c("sci_nam", "gene_name", "position_start", "position_end",
                    "sequence", "accession_version", "create_date",
                    "download_date")

            # all filtering steps again
            cleaned <- standardise_gene_names(GB_data_details)
            cleaned <-
                standardise_spp_names(cleaned, new_names_file = new_names_file)
            cleaned <- remove_duplicates(cleaned)

            those_of_interest <-
                gene_of_interest(gene, data = cleaned)

            return(those_of_interest)
        }
    }
