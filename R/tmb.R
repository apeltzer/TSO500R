#' Read in a TMB_trace.tsv file and store as an object
#'
#' @description Read in a TMB_trace.tsv file
#'
#' @param tmb_file_path a file path to a TMB_trace.tsv file
#'
#' @return A tmb.variant.output object
#' 
#' @export
tmb <- function(tmb_file_path){
  new_tmb_variant_output(tmb_file_path)
}

#' Read in a batch of TMB_trace.tsv files into a list
#'
#' @param tmb_directory a file path to a directory containing one of more
#' TMB_trace.tsv files
#'
#' @return A named list of data frame objects
#'
#' @export
#'
#' @importFrom readr read_tsv
#' @importFrom tidyr unnest
#' @importFrom dplyr mutate relocate
#' @importFrom stringr str_replace
read_tmb_trace_data <- function(tmb_directory) {
  tmb_files <- list.files(
    path = tmb_directory,
    pattern = "*TMB_Trace\\.tsv$|tmb.trace\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )

  tmb_data = tibble(file = tmb_files) |>
    mutate(data = lapply(file, read_tsv)) |>
    unnest(data) |>
    mutate(sample_id = str_replace(basename(file), "_TMB_Trace.tsv|.tmb.trace.tsv", "")) |>
    select(-file) |>
    relocate(sample_id)

  tmb_data
}

#' Read in a batch of *tmb.json files into a list
#'
#' @param tmb_directory a file path to a directory containing one of more
#' *tmb.json files
#'
#' @return A named list of data frame objects
#'
#' @export
#'
#' @importFrom jsonlite read_json
#' @importFrom tidyr unnest_wider
#' @importFrom dplyr mutate relocate
#' @importFrom stringr str_replace
read_tmb_details_data <- function(tmb_directory) {
  tmb_files <- list.files(
    path = tmb_directory,
    pattern = "*tmb.json",
    full.names = TRUE
  )

  tmb_data <- tibble(file = tmb_files) |>
    mutate(data = lapply(tmb_files, read_json)) |>
    unnest_wider(data) |>
    unnest_wider(Settings) |>
    mutate(sample_id = str_replace(basename(file), ".tmb.json", "")) |>
    select(-file) |>
    relocate(sample_id)

  tmb_data
}

#' Read in a batch of *tmb.metrics.csv files into a list
#'
#' @param tmb_directory a file path to a directory containing one of more
#' *tmb.metrics.csv files
#'
#' @return A named list of data frame objects
#'
#' @export
#'
#' @importFrom dplyr mutate select relocate
#' @importFrom tidyr unnest_longer unnest pivot_wider
#' @importFrom tibble tibble
read_tmb_details_data_csv <- function(tmb_directory) {
  tmb_files <- list.files(
    path = tmb_directory,
    pattern = "*tmb.metrics.csv",
    full.names = TRUE
  )

  tmb_data <- tibble(file = tmb_files) |>
    mutate(data = lapply(tmb_files, read.table, header = FALSE, sep = ",")) |>
    unnest_longer(data) |>
    unnest(data) |>
    select(-c(V1, V2)) |>
    pivot_wider(names_from = V3, values_from = V4) |>
    mutate(sample_id = str_replace(basename(file), ".tmb.metrics.csv", "")) |>
    select(-file) |>
    relocate(sample_id)

  tmb_data
}
