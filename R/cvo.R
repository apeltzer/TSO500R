#' Read in a CombinedVariantOutput.tsv file and store as an object
#'
#' @description Read in a CombinedVariantOutput.tsv file
#'
#' @param cvo_file_path a file path to a CombinedVariantOutput.tsv file
#' @param local_app specifies whether data has been generated using the local app
#' @param ctdna specifies whether data has been generated using ctdna workflow
#'
#' @return A combined.variant.output object
#' @export
#'
cvo <- function(cvo_file_path, local_app = FALSE, ctdna = FALSE) {
  new_combined_variant_output(cvo_file_path, local_app, ctdna)
}

#' Constructor function for combined.variant.output objects
#' Not to be called directly
#'
#' @param cvo_file_path a file path to a CombinedVariantOutput.tsv file
#' @param local_app specifies whether data has been generated using the local app
#' @param ctdna specifies whether data has been generated using ctdna workflow
#'
#' @return A combined.variant.output object
#'
#' @importFrom readr read_file
#' @importFrom stringr str_split
#' @importFrom purrr map
#' @importFrom stats na.omit
new_combined_variant_output <- function(cvo_file_path, local_app = FALSE, ctdna = FALSE) {

  ANALYSIS_DETAILS_STRING <- "Analysis Details"
  GENE_AMP_STRING <- c("Gene Amplifications", "Copy Number Variants")
  PERCENT_MSI_STRING <- "Percent Unstable MSI Sites|SUM_JSD"

  cvo_file <- read_file(cvo_file_path)
  split_cvo_string <- str_split(string = cvo_file, pattern = "\\[") |> unlist()

  # parse analysis details, sequencing run details, TMB, and MSI part of provided CombinedVariantOutput file
  start_info <- which(grepl(ANALYSIS_DETAILS_STRING,split_cvo_string))
  end_info <- which(grepl(PERCENT_MSI_STRING, split_cvo_string))

  # handle the parts of the file that are structured as key-value pairs
  # i.e. metadata and TMB/MSI sections
  records <- map(split_cvo_string[start_info:end_info], parse_cvo_record)
  names(records) <- c("analysis_details", "sequencing_run_details", "tmb", "msi")

  # handle the parts of the file that are structured as tabular data
  # i.e. gene amplifications, splice variants, fusions, and small variants
  start_data <- na.omit(pmatch(GENE_AMP_STRING, split_cvo_string))
  end_data <- length(split_cvo_string)

  tables <- map(split_cvo_string[start_data:end_data], parse_cvo_table)

  # the number of tables is different between local app and DRAGEN analysis pipeline
  if (local_app) {
    names(tables) <- c("gene_amplifications", "splice_variants", "fusions", "small_variants")
  } else {
    if (ctdna) {
      names(tables) <- c("copy_number_variants", "dna_fusions", "small_variants")
    } else {
      names(tables) <- c("gene_amplifications", "splice_variants", "fusions", "exon_level_cnvs", "small_variants")
    }
  }

  return(structure(c(records, tables), class = "combined.variant.output"))
}

#' Validator function for combined.variant.output constructor
#' Not to be called directly
#' NOT IMPLEMENTED
#'
#' @return boolean
validate_tso500 <- function() {}

#' Read in a batch of CombinedVariantOutput.tsv files into a list
#'
#' @param cvo_directory a file path to a directory containing one of more CombinedVariantOutput.tsv files
#' @param local_app specifies whether data has been generated using the local app
#' @param ctdna specifies whether data has been generated using ctdna workflow
#'
#' @return A named list of combined.variant.output objects
#'
#' @export
read_cvo_data <- function(cvo_directory, local_app = FALSE, ctdna = FALSE) {
  cvo_files <- list.files(
    path = cvo_directory,
    pattern = "*CombinedVariantOutput\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
  )
  cvo_data <- purrr::map(cvo_files, cvo, local_app, ctdna)
  names(cvo_data) <- purrr::map(cvo_data, ~ ifelse(ctdna, .x$analysis_details$dna_sample_id,
                                                   .x$analysis_details$pair_id))
  cvo_data
}

#' Extract small variants from combined.variant.output object and
#' return in data frame format
#'
#' @param cvo_obj cvo_obj
#'
#' @return A data frame of small variants
#'
#' @export
get_small_variants <- function(cvo_obj, ...) {
  UseMethod("get_small_variants", cvo_obj)
}

#' Extract gene amplifications from combined.variant.output object and
#' return in data frame format
#'
#' @param cvo_obj cvo_obj
#'
#' @return A data frame of gene amplifications
#'
#' @export
get_gene_amplifications <- function(cvo_obj, ...) {
  UseMethod("get_gene_amplifications", cvo_obj)
}

#' Extract splice variants from combined.variant.output object and
#' return in data frame format
#'
#' @param cvo_obj cvo_obj
#'
#' @return A data frame of splice variants
#'
#' @export
get_splice_variants <- function(cvo_obj, ...) {
  UseMethod("get_splice_variants", cvo_obj)
}

#' Extract fusions from combined.variant.output object and
#' return in data frame format
#'
#' @param cvo_obj cvo_obj
#'
#' @return A data frame of fusions
#'
#' @export
get_fusions <- function(cvo_obj, ...) {
  UseMethod("get_fusions", cvo_obj)
}

#' Get small variants from combined.variant.output object
#'
#' @param cvo_obj cvo_obj
#' @return A data frame
#' @method get_small_variants combined.variant.output
#'
#' @export
#'
#' @importFrom dplyr mutate select
#' @importFrom tidyr everything
get_small_variants.combined.variant.output <- function(cvo_obj) {
  suppressWarnings(
    if(all(is.na(cvo_obj$small_variants))){
      small_variant_df <- data.frame()
    } else {
      small_variant_df <- cvo_obj$small_variants |>
        mutate(sample_id = ifelse(is.null(cvo_obj$analysis_details$pair_id),
                                  cvo_obj$analysis_details$dna_sample_id,
                                  cvo_obj$analysis_details$pair_id)) |>
        select(sample_id, everything())
  }
  )
  return(small_variant_df)
}

#' Get gene amplifications from combined.variant.output object
#'
#' @param cvo_obj cvo_obj
#' @return A data frame
#' @method get_gene_amplifications combined.variant.output
#'
#' @export
#'
#' @importFrom dplyr mutate select
#' @importFrom tidyr everything
get_gene_amplifications.combined.variant.output <- function(cvo_obj) {
  suppressWarnings(
    if (all(is.na(cvo_obj$gene_amplifications)) & all(is.na(cvo_obj$copy_number_variants))) {
      gene_amplification_df <- data.frame()
    } else {
      if (is.null(cvo_obj$gene_amplifications)) {
        gene_amplification_df <- cvo_obj$copy_number_variants
      } else {
        gene_amplification_df <- cvo_obj$gene_amplifications
      }

      gene_amplification_df <- gene_amplification_df |>
        mutate(sample_id = ifelse(is.null(cvo_obj$analysis_details$pair_id),
                                  cvo_obj$analysis_details$dna_sample_id,
                                  cvo_obj$analysis_details$pair_id)) |>
        select(sample_id, everything())
    }
  )
  return(gene_amplification_df)
}

#' Get splice variants from combined.variant.output object
#'
#' @param cvo_obj cvo_obj
#' @return A data frame
#' @method get_splice_variants combined.variant.output
#'
#' @export
#'
#' @importFrom dplyr mutate select
#' @importFrom tidyr everything
get_splice_variants.combined.variant.output <- function(cvo_obj) {
  suppressWarnings(
    if (all(is.na(cvo_obj$splice_variants))){
      splice_variant_df <- data.frame()
    } else {
      splice_variant_df <- cvo_obj$splice_variants |>
        mutate(sample_id = ifelse(is.null(cvo_obj$analysis_details$pair_id),
                                  cvo_obj$analysis_details$dna_sample_id,
                                  cvo_obj$analysis_details$pair_id)) |>
        select(sample_id, everything())
    }
  )
  return(splice_variant_df)
}

#' Get fusions from combined.variant.output object
#'
#' @param cvo_obj cvo_obj
#' @return A data frame
#' @method get_fusions combined.variant.output
#'
#' @export
get_fusions.combined.variant.output <- function(cvo_obj){
  suppressWarnings(
    if (all(is.na(cvo_obj$fusions)) & all(is.na(cvo_obj$dna_fusions))) {
      fusion_df <- data.frame()
    } else {
      if(is.null(cvo_obj$fusions)) {
        fusion_df <- cvo_obj$dna_fusions
      } else {
        fusion_df <- cvo_obj$fusions
      }

      fusion_df <- fusion_df |>
        dplyr::mutate(sample_id = ifelse(is.null(cvo_obj$analysis_details$pair_id),
                                         cvo_obj$analysis_details$dna_sample_id,
                                         cvo_obj$analysis_details$pair_id)) |>
        dplyr::select(sample_id, tidyr::everything())
    }
  )
  return(fusion_df)
}

#' Helper function to parse key-value lines in CombinedVariantOutput.tsv
#'
#' @param record_string the record content as string
#'
#' @return char vector
#'
#' @importFrom stringr str_split str_remove str_detect
#' @importFrom purrr map map_chr
#' @importFrom janitor make_clean_names
parse_cvo_record <- function(record_string) {
  intermediate <- record_string |>
    trim_cvo_header_and_footer() |>
    str_split("\n") |>
    unlist() |>
    str_remove("\\t$") |>
    str_split("\\t") |>
    rapply(., function(x) ifelse(x=="NA",NA,x), how = "replace") # replace all string NAs with NA to avoid warnings from as.numeric

  if(str_detect(record_string, "\\[TMB\\]|\\[MSI\\]")){
    record <- map(intermediate, ~ as.numeric(.x[2]))
  } else {
    record <- map(intermediate, ~ .x[2])
  }

  record_names <- map_chr(intermediate, ~ .x[1])
  names(record) <- make_clean_names(record_names)
  return(record)
}

#' Helper function to parse tabular data in CombinedVariantOutput.tsv
#'
#' @param table_string the table as string
#'
#' @return data.frame
#'
#' @importFrom stringr str_extract str_detect str_replace_all
parse_cvo_table <- function(table_string) {
  intermediate <- table_string |> trim_cvo_header_and_footer()
  header_line <- str_extract(intermediate, ".+\n")

  if (str_detect(header_line, "\\t\\n")) {
    intermediate <- str_replace_all(
                                    string = intermediate,
                                    pattern = "\\t\\n",
                                    replacement = "\n")
  }

  table_data <- intermediate |> handle_empty_cvo_table_values()
  return(table_data)
}

#' Helper function to handle empty rows
#' in CombinedVariantOutput.tsv tabular data
#'
#' @param intermediate_tbl intermediate table string
#'
#' @return data.frame
#'
#' @importFrom stringr str_detect str_replace str_length
#' @importFrom utils read.table
#' @importFrom janitor clean_names
handle_empty_cvo_table_values <- function(intermediate_tbl) {
  if (str_detect(string = intermediate_tbl, pattern = "\\nNA$")) {
    cleaned_string <- str_replace(intermediate_tbl, "\\nNA", "")
    if (str_length(cleaned_string) > 0) {
      df <- read.table(text = cleaned_string, sep = "\t", header = TRUE, fill = TRUE)
      return(df)
    } else {
      return(NA)
    }
  } else {
    to_clean <- read.table(text = intermediate_tbl, sep = "\t", header = TRUE, fill = TRUE)
    clean_name_df <- clean_names(to_clean)
    return(clean_name_df)
  }
}

#' Helper function to remove header and footer from CombinedVariantOutput.tsv
#'
#' @param string string with file content
#'
#' @return char vector
trim_cvo_header_and_footer <- function(string) {
  string |>
    stringr::str_remove(".+\\t\\t\\n") |>
    stringr::str_remove("[\\n\\t]+$")
}
