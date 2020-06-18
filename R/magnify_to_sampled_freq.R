#' Multiplies sequence data up to sampled frequency, aligns and summarises data
#'
#' \code{magnify_to_sampled_freq}
#'
#'
#' Function multiplies data to sampled frequency, aligns sequences, writes
#' aligned fasta files, plots haplotype network and extracts summary data.
#'
#' Unaligned sequence files are read in and accessions (individual haplotypes)
#' are multiplied by given values in relevant .csv file.   Also, haplotype
#' networks are drawn, and summary data recorded.
#'
#' @param directory_path path to the directory in which to build/call files, 
#' defaults to working directory
#' @param magnify_file_list .csv file with two columns; first column lists
#' accession numbers/accession versions and the second the frequency that
#' haplotype was found at when sampled.
#'
#' @export



magnify_to_sampled_freq <- function(directory_path = getwd(), magnify_file_list) {
    mag_df <- as.data.frame(NULL)

    for (q in 1:length(magnify_file_list)) {
        what <- stringr::str_sub(magnify_file_list[q], start = -3)

        if (what == "csv") {
            spp_name <- gsub("MAGNIFY_", "",
                     gsub(".csv", "", magnify_file_list[q]))

            # read in sequences from alignment file
            alignment_file <- sub("MAGNIFY", "FOR_ALIGNMENT",
                                  sub("csv", "fasta", magnify_file_list[q]))

            my_sequence <-
                Biostrings::readDNAStringSet(paste0(directory_path, "/", 
                                                    alignment_file))
            first_alignment <-
                msa::msaClustalW(my_sequence, type = "dna")
            cleaned_matrix <-
                ape::as.matrix.DNAbin(ape::as.DNAbin(first_alignment))
            cleaned_matrix <-
                ape::as.DNAbin(gsub("[^atcg]+", "-", cleaned_matrix))
            cleaned_haplotypes <-
                ips::deleteGaps(cleaned_matrix, gap.max = 1)

            # get the new frequency values for each accession
            sample_frequency <- utils::read.csv(paste0(directory_path, "/",
                                                       magnify_file_list[q]), 
                                                stringsAsFactors = T)

            if (ncol(sample_frequency) == 2) {
                # set up 'cleaned' file to build on
                cleaned <- cleaned_haplotypes

                # for each label = sample_frequency[,1] subset cleaned_haplotypes
                for (a in 1:nrow(sample_frequency)) {
                    if (sample_frequency[a, 2] > 1) {
                        this_sequence <- cleaned_haplotypes[
                                                sample_frequency[a, 1], ]
                        for (d in 1:c(sample_frequency[a, 2] - 1)) {
                            # need to keep rownames unique for BEAST
                            rownames(this_sequence) <-
                                paste0(labels(this_sequence), ".", d)
                            cleaned <- rbind(cleaned, this_sequence)
                        }
                    }
                }

                # get a histogram of sequence length before and after clean up
                mypath <- file.path(paste0(directory_path, "/histograms"),
                                    paste0("hist_", spp_name, ".png"))
                grDevices::png(file = mypath)
                h1 <- graphics::hist(Biostrings::width(my_sequence))
                h2 <- ncol(cleaned)
                graphics::hist(x = 0, xlim = c(c(h2 - 50), c(
                    max(Biostrings::width(my_sequence)) + 50)), ylim = c(0,
                    max(h1$counts)), ylab = "Frequency", xlab = "seq lengths",
                    main = spp_name, border = "white")
                graphics::abline(v = ncol(cleaned), col = "red")
                graphics::plot(h1, add = T)
                grDevices::dev.off()

                # Write out as fasta file
                aligned_file <- sub("FOR_ALIGNMENT", "ALIGNED",
                        sub("fasta", "fas", alignment_file))

                ape::write.FASTA(cleaned, 
                                 file = paste0(directory_path, "/", aligned_file))

                # add this info to the original info_df

                # get more details about the data split on the acession verrsion
                # ending; normally ends [.1] but not always, hence [:digit:]
                haplo_freq <- pegas::haploFreq(cleaned,
                                               split = "[[:punct:][:digit:]]")
                haplo_number <- length(haplo_freq)
                sequence_length <- length(cleaned) / nrow(cleaned)
                sample_size <- nrow(cleaned)

                # mjn doesnt work with 'raw' data - can't handle repeate
                # haplotypes; find a single version of each haplotype
                haps <- pegas::haplotype(cleaned)
                the_haplonet <- pegas::haploNet(haps)
                max_step <- max(the_haplonet[, "step"])

                mypath <- file.path(paste0(directory_path, "/network_diagrams"),
                              paste0("Net_", spp_name, ".png"))
                grDevices::png(file = mypath)
                # inclusion of 'size =' sets size of circles to hap freq values
                graphics::plot(the_haplonet, size = attr(the_haplonet, "freq"),
                    fast = FALSE)
                grDevices::dev.off()


                # for this one file - bind all info together
                mag_df2 <- as.data.frame(cbind(spp_name, haplo_number,
                            sequence_length, sample_size, max_step))
                # add to larger dataframe
                mag_df <- as.data.frame(rbind(mag_df, mag_df2))
            } else {
                stop("Individual .csv files not formated as expected")
            }


        } else {
            stop("magnify_file_list not list of csv file names as expected")
        }
    }
    return(mag_df)
}
