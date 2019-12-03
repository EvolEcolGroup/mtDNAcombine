#' Find instances where all accessions from one paper are unique haplotypes
#'
#' \code{hapfreq_from_paper}
#'
#' This function finds the assocatied paper in instances where every accession
#' from one study represents a novel haplotype, likely indicating data is not at
#' sampled frequency.
#'
#' Where every accession is a novel haplotype for the given species/gene this
#' function finds the associated paper.  In such instances the original paper
#' will need to be read to find values to multiple haplotypes by in order to
#' reach original sampled frequency.
#'
#' @param freq_by_study matrix created by \code{haploFreq} function (pegas) using
#' categorical variable formed of first 5 characters of accession number.
#' @param new_rownam character string of accession versions with first five
#' characters repeated and pasted with an underscore in front of full string to
#' create a two part row name ID
#' @param max_haps_found_together Threshold value for number of times a haplotype
#' can be found repeated in the uploads from one study before paper is assumed to
#' have submitted ever sample rather than unique haplotypes.
#'
#' @export



# what sequences never match other accessions from the same study
hapfreq_from_paper <- function(freq_by_study, new_rownam, max_haps_found_together) {
    check_these <- NULL
    check_these1 <- NULL
    check_hap_or_samp <- as.data.frame(NULL)
    check_hap_or_samp1 <- as.data.frame(NULL)
    haplotypes_only <- NULL
    for (i in 1:ncol(freq_by_study)) {
        # has this paper given more than one accession for this species/gene
        if (as.numeric(sum(freq_by_study[, i])) > 1) {
            # drop repeated numbers - only care about max number of times the
            # same haplotype has been found
            freqs_found <- as.numeric(max(freq_by_study[, i]))
            # if, for a single paper, haplotypes were never found more than once
            # it is likely these may not be at sampled freq
            if (freqs_found <= max_haps_found_together) {
               check_these1 <- as.data.frame(cbind(labels(freq_by_study[1, i])))
               check_these <- as.data.frame(rbind(check_these, check_these1))
            }
        }
    }

    if (!is.null(nrow(check_these)) == TRUE) {

        # get the full accession numbers for these cases haplotypes_only <- NULL
        for (p in 1:nrow(check_these)) {
            tmp <- as.data.frame(substr(new_rownam[grep(check_these[p, 1],
                        new_rownam)], 7, nchar(new_rownam[grep(
                            check_these[p, 1], new_rownam)])))
            haplotypes_only <- rbind(haplotypes_only, tmp)
        }


        # get associated paper
        for (i in 1:nrow(haplotypes_only)) {
            accession_no <- as.character(haplotypes_only[i, ])
            connection <- curl::curl(paste0(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                accession_no, "&rettype=gb&retmode=xml"))
            full_xmlTree <- readLines(connection)
            # close the connection
            close(connection)
            full_xmlTree <- XML::xmlInternalTreeParse(full_xmlTree, asText = T)
            paper <- XML::xmlValue(XML::getNodeSet(
                full_xmlTree, "//GBReference_title")[[1]])
            sci_nam <- XML::xmlValue(XML::xpathApply(XML::getNodeSet(
                full_xmlTree,
                "//GBQualifier[GBQualifier_name/text()='organism']")[[1]],
                ".//GBQualifier_value")[[1]])
            check_hap_or_samp1 <- as.data.frame(cbind(paper, sci_nam))
            check_hap_or_samp <- as.data.frame(rbind(check_hap_or_samp,
                                                     check_hap_or_samp1))
        }
        check_hap_or_samp <- as.data.frame(unique(check_hap_or_samp))
        colnames(check_hap_or_samp) <- c("paper", "species")

    }

    return(check_hap_or_samp)
}
