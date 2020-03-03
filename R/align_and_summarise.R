#' Aligns sequence data, plots haplotype networks and summarises the data
#'
#' \code{align_and_summarise}
#'
#' Function aligns sequence data, plots two types of diagnostic plots (haplotype
#' networks and sequence size histograms), and summarises the data. Two directories
#' for storing these outputs created automatically,"./network_diagrams" and
#' "./histograms". Where only representative haplotypes uploaded, rather than
#' individual samples, populations recorded for later inspection.
#'
#' Unaligned sequence files are read in, aligned sequences are written out as
#' files.  Histograms of the original sequence length range and new aligned
#' sequence length is drawn and stored out.  Haplotype networks are also drawn
#' and written out.  Summary data are recorded. If only one version of each
#' haplotype ever found, assumed that accessions might represent haplotypes
#' not sample frequency, populations flagged for further inspection.
#'
#' @param directory_path path to the directory in which to build/call files, 
#' defaults to working directory
#' @param alignment_files list of .fas/.fasta files in directory that need
#' alignment, often signaled by the pattern "FOR_ALIGNMENT" in file name
#' @param max_haps_found_together input value for \code{hapfreq_from_paper} function.
#' Threshold value for number of times a haplotype can be found repeated in the
#' uploads from one study before paper is assumed to have submitted every sample
#' rather than unique haplotypes.
#' @param minbp minimum length (base pairs) for a sequence to be retained in the
#' data set
#'
#' @export

align_and_summarise <- function(directory_path = getwd(), alignment_files,
                                max_haps_found_together, minbp) {
    dir.create(paste0(directory_path, "/network_diagrams"))
    dir.create(paste0(directory_path, "/histograms"))
    info_df <- as.data.frame(NULL)
    more_info_df <- NULL
    for (a in 1:length(alignment_files)) {
        spp_name <- gsub("FOR_ALIGNMENT_",
                         "", gsub(".fasta", "", alignment_files[a]))
        # read in sequences from first file
        my_sequence_init <-
            Biostrings::readDNAStringSet(paste0(directory_path, "/", 
                                                alignment_files[a]))
        my_sequence <-
            my_sequence_init[Biostrings::width(my_sequence_init) > minbp]
        first_alignment <-
            msa::msaClustalW(my_sequence, type = "dna")
        cleaned_matrix <-
            ape::as.matrix.DNAbin(ape::as.DNAbin(first_alignment))

        # need to set any ambiguous base calls to '-' so they get swept up
        # as well convert to class()=DNAbin for
        # deleteGaps
        cleaned_matrix <-
            ape::as.DNAbin(gsub("[^atcg]+", "-", cleaned_matrix))
        # deleteGaps deletes any '-' characters.
        # default is gap.max=nrow(x) -4 => a string of four blanks.
        cleaned <- ips::deleteGaps(cleaned_matrix, gap.max = 1)

        # set up to ask if data has been loaded as haplotypes or samples
        add_IDs <- substr(rownames(cleaned), 1, 5)
        old_rownam <- rownames(cleaned)
        new_rownam <- paste(add_IDs, rownames(cleaned), sep = "_")
        rownames(cleaned) <- new_rownam
        # find a single version of each haplotype
        haps <- pegas::haplotype(cleaned)
        freq_by_study <-
            pegas::haploFreq(cleaned, what = 1, haplo = haps)


        # what accessions never come up with other accessions from
        # the same study
        check_these_papers <-
            hapfreq_from_paper(
                freq_by_study = freq_by_study,
                new_rownam = new_rownam,
                max_haps_found_together = max_haps_found_together
            )


        if (nrow(check_these_papers) > 0) {
            more_info_df <-
                as.data.frame(rbind(more_info_df, check_these_papers))

            # create a new df with the accession number and a freq. column,
            # default value of 1
            magnify_haplotypes <- as.data.frame(old_rownam)
            magnify_haplotypes$freq <- 1
            colnames(magnify_haplotypes) <-
                c("accession_version", "freq")
            magnify_file <- sub("FOR_ALIGNMENT", "MAGNIFY",
                                sub("fasta", "csv", alignment_files[a]))
            # write out a csv file to manually fill with hap freq. data
            utils::write.csv(magnify_haplotypes, 
                             file = paste0(directory_path, "/", magnify_file), 
                             quote = FALSE, row.names = FALSE)

        } else {
            # some downstream programmes need file end .FAS not .FASTA
            aligned_file <- sub("FOR_ALIGNMENT", "ALIGNED",
                                sub("fasta", "fas", alignment_files[a]))
            rownames(cleaned) <- old_rownam
            ape::write.FASTA(cleaned, 
                             file = paste0(directory_path, "/", aligned_file))


            # get some more details about the data

            # split on the acession version ending;
            # normally xxxxx[.1] but can be [.2] or other number so used
            # [:digit:]
            haplo_freq <-
                pegas::haploFreq(cleaned, split = "[[:punct:][:digit:]]")
            haplo_number <- length(haplo_freq)
            sequence_length <- length(cleaned) / nrow(cleaned)
            sample_size <- nrow(cleaned)

            the_haplonet <- pegas::haploNet(haps)
            max_step <- max(the_haplonet[, "step"])

            # get a histogram of sequence length before and after clean up
            mypath <- file.path(paste0(directory_path, "/histograms"),
                                paste0("hist_", spp_name, ".png"))
            grDevices::png(file = mypath)
            h1 <- graphics::hist(Biostrings::width(my_sequence_init))
            h2 <- min(Biostrings::width(my_sequence_init))
            graphics::hist(x = 0, xlim = c(c(h2 - 50),
                        c(max(Biostrings::width(my_sequence)) + 50)),
                        ylim = c(0, max(h1$counts)), ylab = "Frequency",
                        xlab = "seq lengths", main = spp_name, border = "white")
            graphics::abline(v = ncol(cleaned), col = "red")
            graphics::plot(h1, add = T)
            grDevices::dev.off()

            # plot network diagram
            mypath <- file.path(paste0(directory_path, "/network_diagrams"),
                          paste0("Net_", spp_name, ".png"))
            grDevices::png(file = mypath)
            # inclusion of 'size =' sets size of circles to hap freq values
            graphics::plot(the_haplonet,
                           size = attr(the_haplonet, "freq"), fast = FALSE)
            grDevices::dev.off()


            # for this one file - bind all info together
            info_df2 <-
                as.data.frame(cbind(spp_name, haplo_number, sequence_length,
                    sample_size, max_step))
            # add to larger dataframe
            info_df <- as.data.frame(rbind(info_df, info_df2))
        }
    }

    utils::write.csv(info_df, paste0(directory_path, "/", "Info_df.csv"), 
                     row.names = FALSE)
    utils::write.csv(more_info_df, paste0(directory_path, "/", "More_info_df.csv"), 
                     row.names = FALSE)

}
