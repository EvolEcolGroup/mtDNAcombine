#' Gets details of individual samples
#'
#' \code{get_sample_details}
#'
#' This function gets basic background details of the sample from the xml tree
#' input
#'
#' Function designed to sit within \code{build_GB_dataframe}. In general, the
#' full_xmlTree should be provided by earlier stages of the
#' \code{build_GB_dataframe} function. Extracts details about the accession,
#' such as species name.
#'
#' @param full_xmlTree xml tree downloaded from NCBI website and parsed to R
#' structure with \code{xmlInternalTreeParse} by earlier stages of the
#' \code{build_GB_dataframe} function.
#'
#' @export



get_sample_details <- function(full_xmlTree) {
    feat_key <- as.data.frame(NULL)
    feat_key_df <- as.data.frame(NULL)

    if (inherits(full_xmlTree, "XMLInternalDocument")) {
        # check xml input is correct how many nodes in xml
        no_nodes <-
            length(XML::getNodeSet(full_xmlTree, "//GBFeature[GBFeature_key]"))
        # loop through number of times = number of nodes
        for (v in 1:no_nodes) {
            feat_key <- XML::getNodeSet(full_xmlTree,
                                        "//GBFeature[GBFeature_key]")[[v]]
            feat_key <- XML::xmlValue(XML::xpathApply(feat_key,
                                                      ".//GBFeature_key")[[1]])
            feat_key_df <- as.data.frame(rbind(feat_key_df, feat_key))
            feat_key_df[sapply(feat_key_df, is.factor)] <-
                lapply(feat_key_df[sapply(feat_key_df, is.factor)],
                       as.character)
        }

        feat_key_count <- as.data.frame(table(feat_key_df))
        # inclusion of 'source' confuses later steps, remove
        feat_key_count <- subset(feat_key_count,
                                 feat_key_count$feat_key_df != "source")
        feat_key_count$feat_key_df <- as.character(feat_key_count$feat_key_df)

        # extract basic info about the sample itself from the first feature
        sci_nam <- XML::getNodeSet(full_xmlTree,
                    "//GBQualifier[GBQualifier_name/text()='organism']")[[1]]
        sci_nam <- XML::xmlValue(XML::xpathApply(sci_nam,
                                                 ".//GBQualifier_value")[[1]])
        create_date <- XML::getNodeSet(full_xmlTree, "//GBSet/GBSeq")[[1]]
        create_date <- XML::xmlValue(XML::xpathApply(create_date,
                                                ".//GBSeq_create-date")[[1]])
        accession_version <- XML::getNodeSet(full_xmlTree, "//GBSet/GBSeq")[[1]]
        accession_version <- XML::xmlValue(XML::xpathApply(accession_version,
                                          ".//GBSeq_accession-version")[[1]])
        # store date data accessed for later indexing/incase of updates etc
        download_date <- as.character.Date(Sys.Date())


        output <- as.list(c(feat_key_count, sci_nam, create_date,
                            accession_version, download_date))
        names(output) <- c("feature_name", "freq", "sci_nam", "create_date",
                           "accession_version", "download_date")

        return(output)
    } else {
        stop("Input full_xmlTree not in a valid format")
    }

}
