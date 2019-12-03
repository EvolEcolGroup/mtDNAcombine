#' Identifies orignial accession number of given haplotype index number
#'
#' \code{outliers_dropped}
#'
#' Finds aligned sequence data for the populations with samples that are over a
#' user defined value for maximum number of mutations on a branch before being
#' considered to be an outlier.  Haplotype networks are plotted and outlier
#' sequences are dropped from fasta files.
#'
#' Function identifies populations that have samples seperated by a number of
#' mutations greater than a user defined value for being considered an outlier.
#' Outlier sequences are dropped from fasta files, haplotype networks repolotted
#' and new files written out with prefix `new_`.
#'
#' @param max_mutations user defined value for number of mutatiosn allowed before
#' haplotype/sample considered an outlier
#' @param info_df expected to be info_df or structured as info_df built from
#' earlier functions. Either read in file or use object of correct format.
#'
#' @export




# user defined value for number of mutations allowed before haplotype/sample
# considered an outlier
outliers_dropped <- function(max_mutations, info_df) {
    dropped <- NULL
    split_pop <- NULL
    if (inherits(info_df, "character")) {
        df_in <- utils::read.csv(info_df)
    } else {
        if (inherits(info_df, "data.frame")) {
            df_in <- info_df
        } else {
            stop("info_df needs to be an exisiting object or .csv file name")
        }
    }

    has_outliers <-
        df_in[as.numeric(df_in$max_step) >= max_mutations, ]
    # reconstruct aligned file names
    has_outliers$file_name <- paste0("ALIGNED_", has_outliers$spp_name, ".fas")


    for (p in 1:nrow(has_outliers)) {
        # load sequence data
        outlier1 <- ips::read.fas(has_outliers$file_name[p])
        # get unique haplotypes
        haps <- pegas::haplotype(outlier1)
        the_haplonet <- pegas::haploNet(haps)

        n <- length(the_haplonet[, "step"])
        branch_length <- max(the_haplonet[, "step"])
        if (branch_length >= max_mutations) {
            # record accession number for this oultier
            outlier_number <- the_haplonet[n, 1]
            outlier_accession <- (identify_outlier(haplotypes = haps,
                        fasta_data = outlier1, outlier_number = outlier_number))
            # check frequency of this outlier haplotype
            freq_of_outlier <- attr(the_haplonet, "freq")[the_haplonet[n, 1]]
            if (freq_of_outlier == 1) {
                split_pop <- rbind(gsub("ALIGNED_", "new_",
                        gsub(".fas", "", has_outliers$file_name[p])))
                # remove this/these accessions from the sequence file
                outlier2 <- outlier1[!labels(outlier1) %in% outlier_accession, ]

                # update to new cut
                haps <- pegas::haplotype(outlier2)
                the_haplonet <- pegas::haploNet(haps)
                branch_length <- max(the_haplonet[, "step"])




                ape::write.FASTA(
                    outlier2, file = paste0("new_", has_outliers$file_name[p]))

                # get some more details to update info_df
                spp_name <-
                    gsub("ALIGNED_",
                         "new_",
                         gsub(".fas", "", has_outliers$file_name[p]))
                dropped <- rbind(cbind(outlier_accession, spp_name))

                # split on the acession verrsion ending; normally ends [.1]
                # but not always, hence [:digit:]
                haplo_freq <-
                    pegas::haploFreq(outlier2, split = "[[:punct:][:digit:]]")
                haplo_number <- length(haplo_freq)
                sequence_length <- length(outlier2) / nrow(outlier2)
                sample_size <- nrow(outlier2)
                max_step <- max(the_haplonet[, "step"])

                # plot and store images out
                mypath <- file.path("./network_diagrams",
                              paste0("Net_", spp_name, ".jpeg"))
                grDevices::jpeg(file = mypath)
                # inclusion of 'size =' sets size of circles to hap freq values
                graphics::plot(the_haplonet, size = attr(the_haplonet, "freq"),
                    fast = FALSE)
                grDevices::dev.off()

                # for this one file - bind all info together
                info_df2 <- as.data.frame(cbind(spp_name, haplo_number,
                            sequence_length, sample_size, max_step))
                # add to larger dataframe
                updating_info_df(original_df = info_df, info_df2)
                return(dropped)

            } else {
                print(
                    "N.B. Branch lengths greater than max_mutations found,
                    however, multiple samples of these haplotypes. Manual
                    inspection of networks needed"
                )
            }
        }
    }
}
