#' Gets details of all available genes associated with accession in question
#'
#' \code{get_gene_details}
#'
#' This function extracts details and sequence data of genes available for each
#' accession from an XML tree
#'
#' Function designed to sit within \code{build_genbank_df}. In general, the
#' input files should be provided by earlier stages of the
#' \code{build_genbank_df} function. Extracts details about the accession being
#' investigated, such as species name.
#'
#' @param feature_name Names of genes/features found associated with the accession.
#' Function \code{get_sample_details} outputs `feature_name` value.
#' @param freq Frequency at which each of the `feature_name` genes/features were
#' found within the XML tree for this accession. Function \code{get_sample_details}
#' outputs a `freq` value associated with `feature_name` value
#' @param full_xmlTree xml tree downloaded from NCBI website and parsed to R
#' structure with \code{xmlInternalTreeParse} - should be provided by earlier
#' stages of the \code{build_genbank_df} function.
#'
#' @export


get_gene_details <- function(feature_name, freq, full_xmlTree) {
    output <- NULL
    if (inherits(full_xmlTree, "XMLInternalDocument")) {
        # check xml input is correct re-build feat_key_count
        feat_key_count <- as.data.frame(cbind(feature_name, freq))
        # using feature names from each of the nodes, go through options
        for (z in 1:nrow(feat_key_count)) {
            feature <- feat_key_count[z, 1]
            this_feat_key <-
                subset(feat_key_count,
                       feat_key_count$feature_name == feat_key_count[z, 1])
            this_feat_key$freq <-
                as.numeric(as.character(this_feat_key$freq))
            for (k in 1:this_feat_key$freq) {
                # different genes have data stored in different
                # structures/parts of xml tree
                if (feature == "gene") {
                    gene_name <- XML::getNodeSet(
                        full_xmlTree,
                        "//GBQualifier[GBQualifier_name/text()='gene']")[[k]]
                } else {
                    gene_name <-
                        XML::getNodeSet(
                            full_xmlTree,
                            paste("//GBFeature[GBFeature_key/text()=\"",
                                  feature, "\"]", sep = ""))[[k]]
                }

                if (feature %in% c("CDS", "tRNA")) {
                    gene_name <- XML::getNodeSet(
                        gene_name,
                        ".//GBQualifier[GBQualifier_name/text()='product']"
                        )[[1]]
                }

                if (feature == "D-loop") {
                    gene_name <-
                        XML::xmlValue(XML::xpathApply(
                            gene_name, ".//GBFeature_key")[[1]])
                } else {
                    gene_name <-
                        XML::xmlValue(XML::xpathApply(
                            gene_name, ".//GBQualifier_value")[[1]])
                }


                gene_info <-
                    XML::getNodeSet(
                        full_xmlTree,
                        paste("//GBFeature[GBFeature_key/text()=\"",
                            feature, "\"]", sep = ""))[[k]]
                gene_info <-
                    XML::getNodeSet(gene_info,
                                    ".//GBFeature_intervals/GBInterval")[[1]]
                # if full genome or concatinated sequence uploaded will need
                # start and stop position to crop out gene of interest from
                # full sequence - also gives data on how long sequence is
                position_start <-
                    XML::xmlValue(XML::xpathApply(
                        gene_info, ".//GBInterval_from")[[1]])
                position_end <-
                    XML::xmlValue(XML::xpathApply(
                        gene_info, ".//GBInterval_to")[[1]])

                # get full sequence as uploaded to GenBank
                raw_sequences <-
                    XML::getNodeSet(full_xmlTree, "//GBSeq/GBSeq_sequence")[[1]]
                raw_sequences <-
                    XML::xmlValue(XML::xpathApply(
                        raw_sequences, "//GBSeq_sequence")[[1]])
                raw_sequences <- as.character(raw_sequences)

                # get the section of that full sequence you are interested in
                # For loop used to deal with 'complement'
                # values where start stop numbers are reversed
                if (as.numeric(as.character(position_start))
                    > as.numeric(as.character(position_end))) {
                    the_sequence <- substr(raw_sequences,
                            start = as.numeric(as.character(position_end)),
                            stop = as.numeric(as.character(position_start)))
                } else {
                    the_sequence <- substr(raw_sequences,
                            start = as.numeric(as.character(position_start)),
                            stop = as.numeric(as.character(position_end)))
                }

                output2 <- as.data.frame(cbind(gene_name, position_start,
                        position_end, the_sequence))
                output <- as.data.frame(rbind(output, output2))

            }


        }


        names(output) <- c("gene_name", "position_start",
                           "position_end", "the_sequence")
        return(output)
    } else {
        stop("Input full_xmlTree not in a valid format")
    }
}
