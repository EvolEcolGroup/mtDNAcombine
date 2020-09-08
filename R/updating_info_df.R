#' Multiplies sequence data up to sampled frequency, aligns and summarises data
#'
#' \code{updating_info_df}
#'
#'
#' Takes newly created data frame and merges with exisiting data frame of same
#' structure to ensure one master data frame is maintained.
#'
#' Function reads in or takes newly created data frame and adds information to
#' existing data frame of same structure, the combined data frame is the written
#' out as a .csv file as well as being returned by the function.  New data frame
#' is based on files not previously processed because papers uploaded haplotypes
#' not sample frequency, so data needed additional processing before inclusion.
#' This function ensures one master data frame is maintained.
#'
#' @param original_df expected to be info_df or structured as info_df from
#' earlier functions
#' @param new_df expected to be mag_df data frame or structure as mag_df
#' constructed in earlier function.
#' @param file_name defaults to 'Info_df.csv' but user can set name if desired,
#' must include file extension, e.g. '.csv' 
#'
#' @export

updating_info_df <- function(original_df, new_df, file_name = "Info_df.csv") {
    # open appropriate file
    if (inherits(original_df, "character")) {
        df_in <- utils::read.csv(original_df, stringsAsFactors = T)
    } else {
        if (inherits(original_df, "data.frame")) {
            df_in <- original_df
        } else {
            stop("original_df needs to be an existing
                 object or full .csv file name")
        }
    }

    # check structure
    if (all(colnames(new_df) == colnames(df_in))) {
        updated_df <- rbind(df_in, new_df)
        utils::write.csv(updated_df, file_name, row.names = FALSE)
        return(updated_df)
    } else {
        stop("structure of two data frames doesnt match")
    }
}
