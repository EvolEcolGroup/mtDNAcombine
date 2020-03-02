#' Exports downloaded details into seperate .csv files for each species
#'
#' \code{export_sequences}
#'
#' Function takes data frame with column `$sci_nam` and splits according to species
#' name before exporting sequence data from each species as seperate .fasta files.
#'
#' Input data frame is split by species (according to `$sci_nam` column). For each
#' species a separate .fasta file is then written out containing the raw sequence
#' data downloaded from each accession.
#'
#' @param data data frame with `$sci_nam`, `$sequence`, and `$accession_version`
#' column (expected to be `GB_with_SeqDat`, the data frame result of the pipeline
#' so far).
#'
#' @export


export_sequences <- function(data) {
    if (length(unique(
        c("sci_nam", "sequence", "accession_version") %in% colnames(data)
    )) == 1) {
        exporting <- split(data, data$sci_nam, drop = TRUE)
        spp <- names(exporting)


        for (k in 1:length(spp)) {
            lapply(names(exporting), function(y) {
                seqinr::write.fasta(
                    sequences = as.list(exporting[[k]]$sequence),
                    names = exporting[[k]]$accession_version,
                    file.out = paste0("FOR_ALIGNMENT_",
                    as.character(sub(" ", "_", x = spp[k])), ".fasta"),
                    open = "w", nbchar = 60)
            })
        }

    } else {
        stop("data doesn't have necessary columns")
    }

}
