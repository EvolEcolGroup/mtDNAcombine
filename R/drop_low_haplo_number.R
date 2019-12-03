#' Drops populations from main data frame when haplotype number is less than
#' threshold value
#'
#' \code{drop_low_haplo_number}
#'
#' Remove populations that, after alignment, have only haplotype count of less
#' than a threshold value.
#'
#' Find and drop from master list populations that have insufficent numbers of
#' haplotypes to be valuable for further downstream analsysis - value set by user.
#'
#' @param info_df expected to be data frame of same structure as info_df built
#' under earlier functions. Expected column names are `spp_name`, `haplo_number`,
#' `sequence_length`, `sample_size`, `max_step`
#' @param min_haps minimum number of haplotypes before user considers data will
#' be uninformative
#'
#' @export


drop_low_haplo_number <- function(info_df, min_haps) {
    if (inherits(info_df, "data.frame")) {
        expected_colnames <- c("spp_name", "haplo_number", "sequence_length",
                               "sample_size", "max_step")
        if (all(colnames(info_df) %in% expected_colnames)) {
            info_df <-
                subset(info_df, !as.numeric(info_df$haplo_number) < min_haps)
            utils::write.csv(info_df, "Info_df.csv", row.names = FALSE)
            return(info_df)
        } else {
            stop("info_df doesn't have the expected columns")
        }

    } else {
        stop("info_df not a data.frame")
    }
}
