#' Read in a *CopyNumberVariants.vcf file and store as an object
#'
#' @description Read in a *CopyNumberVariants.vcf file
#'
#' @param cnv_file_path a file path to a *CopyNumberVariants.vcf file
#'
#' @return A cnv.output object
#' 
#' @export
cnv <- function(cnv_file_path, local_app = FALSE) {
  new_cnv_output(cnv_file_path)
}

#' Constructor function for combined.cnv.output objects
#' Not to be called directly
#'
#' @param cnv_file_path a file path to a *CopyNumberVariants.vcf file
#' @param local_app specifies whether quality metrics are coming from local app
#'
#' @return A combined.cnv.output object
#'
#' @importFrom dplyr tibble
new_cnv_output <- function(cnv_file_path, local_app = FALSE) {

  cnv_data <- tibble(file = cnv_file_path) |>
    mutate(data = lapply(file, parse_vcf_to_df)) |>
    unnest(data) |>
    mutate(sample_id = str_replace(basename(file), "_CopyNumberVariants.vcf", "")) |>
    select(-file) |>
    relocate(sample_id)

  return(structure(cnv_data, class = "combined.cnv.output"))
}

#' Read in a batch of *CopyNumberVariants.vcf files into a list of CNV objects
#'
#' @param cnv_directory a file path to a directory containing one of more *CopyNumberVariants.vcf files
#' @param local_app specifies whether quality metrics are coming from local app
#'
#' @return A named list of combined.cnv.output objects
#'
#' @export
#'
#' @importFrom purr map set_names
#' @importFrom stringr set_remove
read_cnv_data <- function(cnv_directory, local_app = FALSE) {
  cnv_files <- list.files(
    path = cnv_directory,
    pattern = "*cnv\\.vcf$|*CopyNumberVariants\\.vcf$",
    recursive = TRUE,
    full.names = TRUE
  )
  cnv_data <- map(cnv_files, cnv, local_app)  |>
    set_names(str_remove(basename(cnv_files), "\\.vcf$"))
  cnv_data
}

#' Read in a batch of *CopyNumberVariants.vcf files into one dataframe
#'
#' @param cnv_directory a file path to a directory containing one of more
#' *tmb.json files
#'
#' @return A dataframe with the read CNV data
#'
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom stringr str_replace
#' @importFrom tidyr unnest
#' @importFrom dplyr relocate mutate
summarize_cnv_data <- function(cnv_directory) {
  cnv_files <- list.files(
    path = cnv_directory,
    pattern = "*cnv\\.vcf$|*CopyNumberVariants\\.vcf$",
    recursive = TRUE,
    full.names = TRUE
  )

  cnv_data <- tibble(file = cnv_files) |>
    mutate(data = lapply(file, parse_vcf_to_df)) |>
    unnest(data) |>
    mutate(sample_id = str_replace(basename(file), "_CopyNumberVariants\\.vcf$", "")) |>
    select(-file) |>
    relocate(sample_id)

  cnv_data
}
