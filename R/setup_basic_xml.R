#' Creates basic skeleton xml files
#'
#' \code{setup_basic_xml}
#'
#' Uses \code{babette} package to write out skeleton xml files for each population
#'
#'
#'
#' @param gene_name name of gene the sequence data is from
#' @param aligned_files character string of file names of aligned sequence data
#' fasta files
#'
#' @export


setup_basic_xml <- function(gene_name, aligned_files) {
    if (inherits(gene_name, "character")) {
        if (inherits(aligned_files, "character")) {
            # create a new folder for these raw files
            dir.create("./BEASTfiles")
            wd <- paste0(getwd(), "/BEASTfiles/")

            for (g in 1:length(aligned_files)) {
                # loop through all the alignment files and build the basic .xml
                # N.B. these will need more editing
                beautier::create_beast2_input_file(
                    input_filename = aligned_files[g],
                    output_filename = paste0(wd, gsub(
                        "ALIGNED", gene_name, gsub(".fas", ".xml",
                        aligned_files[g]))), mcmc = beautier::create_mcmc(
                        chain_length = 1e+08, store_every = -1),
                        tree_prior = beautier::create_cbs_tree_prior(
                        id = beautier::get_alignment_id(aligned_files[g]),
                        group_sizes_dimension = 5))
            }

        } else {
            stop("Given list of aligned_files is not a character string")
        }
    } else {
        stop("Given gene_name is not a character string")
    }
}
