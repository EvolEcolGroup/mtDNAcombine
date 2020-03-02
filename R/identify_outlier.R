#' Identifies original accession number of given haplotype index number
#'
#' \code{identify_outlier}
#'
#' Uses fasta format sequence data, associated list of haplotypes, and index
#' number of interest to extract original accession number of outlying sample.
#'
#' Function uses fasta format sequence data, associated list of haplotypes
#' created by pegas \code{haplotype}, and a given index number of interest to
#' extract original accession number of outlying sample.
#'
#' @param haplotypes list of haplotypes created by pegas \code{haplotype}
#' @param fasta_data fasta format sequence data as read in by \code{read.fas}
#' @param outlier_number index number of the haplotype in network that is considered
#' to be an outlier, index numbers are attributed by pegas \code{haploNet}
#'
#' @export



identify_outlier <-
    function(haplotypes, fasta_data, outlier_number) {
        dat <- as.matrix(fasta_data)
        nam <- dimnames(dat)[[1]]
        for (i in 1:dim(haplotypes)[1])
            attr(haplotypes, "dimnames")[[1]][i] <-
            nam[attr(haplotypes, "index")[[i]][1]]

        outlier_ID <- haplotypes[outlier_number, ]
        outlier_accession <- labels(outlier_ID)

        return(outlier_accession)
    }
