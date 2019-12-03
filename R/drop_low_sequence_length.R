#' Drops populations from main data frame when sequence length is less than
#' threshold value
#'
#' \code{drop_low_sequence_length}
#'
#' Remove populations that, after alignment, only have sequence data of length
#' less than a threshold value.
#'
#' Find and drop from master list populations that have insufficent bp to be
#' valuable for further downstream analsysis - value set by user.
#'
#' @param info_df expected to be data frame of same structure as info_df built
#' under earlier functions. Expected column names are `spp_name`, `haplo_number`,
#' `sequence_length`, `sample_size`, `max_step`
#' @param min_length minimum length of sequence before user considers data will
#' be uninformative
#'
#' @export


drop_low_sequence_length <- function(info_df, min_length) {
    if (inherits(info_df, "data.frame")) {
        expected_colnames <- c("spp_name", "haplo_number", "sequence_length",
                               "sample_size", "max_step")
        if (all(colnames(info_df) %in% expected_colnames)) {
            info_df <- subset(info_df,
                              !as.numeric(info_df$sequence_length) < min_length)
            utils::write.csv(info_df, "Info_df.csv", row.names = FALSE)
            return(info_df)
        } else {
            stop("info_df doesn't have the expected columns")
        }

    } else {
        stop("info_df not a data.frame")
    }
}
