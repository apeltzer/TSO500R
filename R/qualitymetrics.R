#' Read in a MetricsOutput.tsv file and store as an object
#'
#' @description Read in a MetricsOutput.tsv file
#'
#' @param metrics_file_path a file path to a MetricsOutput.tsv file
#' @param local_app specifies whether quality metrics are coming from local app
#'
#' @return A quality.metrics.output object
#' 
#' @export
qualitymetrics <- function(metrics_file_path, local_app=FALSE, ctdna=FALSE){
  new_combined_quality_metrics_output(metrics_file_path, local_app, ctdna)
}

#' Constructor function for quality.metrics.output objects
#' Not to be called directly
#'
#' @param metrics_file_path a file path to a MetricsOutput.tsv file
#' @param local_app specifies whether quality metrics are coming from local app
#'
#' @return A quality.metrics.output object
new_combined_quality_metrics_output <- function(metrics_file_path, local_app=FALSE, ctdna=FALSE) {

  qm_file <- readr::read_file(metrics_file_path)
  split_qmo_string <- stringr::str_split(string = qm_file, pattern = "\\[") %>% unlist()
  
  notes_section_index <- ifelse(ctdna, 10, 12)

  # handle the parts of the file that are structured as key-value pairs
  # i.e. header and notes sections
  records <- purrr::map(split_qmo_string[c(2,notes_section_index)], parse_qmo_record)
  names(records) <- c("header", "notes")

  # handle the parts of the file that are structured as tabular data
  # i.e. run qc metrics, analysis status etc. 
  tables <- purrr::map(split_qmo_string[3:(notes_section_index-1)], parse_qmo_table)

  # the order of tables is different between local app and DRAGEN analysis pipeline
  if (local_app) {
    names(tables) <- c("run_qc_metrics", "analysis_status", "dna_qc_metrics", "dna_qc_metrics_snvtmb", "dna_qc_metrics_msi", "dna_qc_metrics_cnv", "dna_expanded_metrics", "rna_qc_metrics", "rna_expanded_metrics")
  }
  else {
    if (ctdna) {
      names(tables) <- c("run_qc_metrics", "analysis_status", "dna_qc_metrics", "dna_qc_metrics_snvtmb", "dna_qc_metrics_msi", "dna_qc_metrics_cnv", "dna_expanded_metrics")
    }
    else {
      names(tables) <- c("run_qc_metrics", "analysis_status", "dna_qc_metrics", "dna_qc_metrics_snvtmb", "dna_qc_metrics_msi", "dna_qc_metrics_cnv", "rna_qc_metrics", "dna_expanded_metrics", "rna_expanded_metrics")
    }
  }

  return(structure(c(records, tables), class = "combined.quality.metrics.output"))
}

#' Validator function for quality.metrics.output constructor
#' Not to be called directly
#' NOT IMPLEMENTED
#'
#' @return boolean
validate_tso500_qc <- function() {}

#' Read in a batch of MetricsOutput.tsv files into a list
#'
#' @param qmo_directory a file path to a directory containing one of more
#' MetricsOutput.tsv files
#' @param local_app specifies whether quality metrics are coming from local app (default: FALSE)
#'
#' @return A named list of combined.quality.metrics.output objects
#' @export
read_qmo_data <- function(qmo_directory, local_app=FALSE, ctdna=FALSE){
  qmo_files <- list.files(
    path = qmo_directory,
    pattern = "*MetricsOutput\\.tsv$",
    full.names = TRUE
  )
  
  qmo_data <- map(qmo_files, qualitymetrics, local_app, ctdna) %>%
    set_names(str_remove(basename(qmo_files), "\\.tsv$")) 

  qmo_data
}

#' Extract run qc metrics from combined.quality.metrics.output object and
#' return in data frame format
#'
#' @param qmo_obj qmo_obj
#'
#' @return A data frame with run qc metrics
#'
#' @export
get_run_qc_metrics <- function(qmo_obj, ...){
  UseMethod("get_run_qc_metrics", qmo_obj)
}

#' Extract the analysis status from combined.quality.metrics.output object and
#' return in data frame format
#'
#' @param qmo_obj qmo_obj
#'
#' @return A data frame with the analysis status
#'
#' @export
get_analysis_status <- function(qmo_obj, ...){
  UseMethod("get_analysis_status", qmo_obj)
}

#' Extract dna qc metrics from combined.quality.metrics.output object and
#' return in data frame format
#'
#' @param qmo_obj qmo_obj
#'
#' @return A data frame with the dna qc metrics
#'
#' @export
get_dna_qc_metrics <- function(qmo_obj, ...){
  UseMethod("get_dna_qc_metrics", qmo_obj)
}

#' Extract dna qc metrics for small variants and tmb from combined.quality.metrics.output object and
#' return in data frame format
#'
#' @param qmo_obj qmo_obj
#'
#' @return A data frame with dna qc metrics (small variants/TMB)
#'
#' @export
get_dna_qc_metrics_snvtmb <- function(qmo_obj, ...){
  UseMethod("get_dna_qc_metrics_snvtmb", qmo_obj)
}

#' Extract dna qc metrics for msi from combined.quality.metrics.output object and
#' return in data frame format
#'
#' @param qmo_obj qmo_obj
#'
#' @return A data frame with dna qc metrics (MSI)
#'
#' @export
get_dna_qc_metrics_msi <- function(qmo_obj, ...){
  UseMethod("get_dna_qc_metrics_msi", qmo_obj)
}

#' Extract dna qc metrics for CNV from combined.quality.metrics.output object and
#' return in data frame format
#'
#' @param qmo_obj qmo_obj
#'
#' @return A data frame with dna qc metrics (CNV)
#'
#' @export
get_dna_qc_metrics_cnv <- function(qmo_obj, ...){
  UseMethod("get_dna_qc_metrics_cnv", qmo_obj)
}

#' Extract expanded dna qc metrics from combined.quality.metrics.output object and
#' return in data frame format
#'
#' @param qmo_obj qmo_obj
#'
#' @return A data frame with extended dna qc metrics
#'
#' @export
get_dna_expanded_metrics <- function(qmo_obj, ...){
  UseMethod("get_dna_expanded_metrics", qmo_obj)
}

#' Extract rna qc metrics from combined.quality.metrics.output object and
#' return in data frame format
#'
#' @param qmo_obj qmo_obj
#'
#' @return A data frame with rna qc metrics
#'
#' @export
get_rna_qc_metrics <- function(qmo_obj, ...){
  UseMethod("get_rna_qc_metrics", qmo_obj)
}

#' Extract expanded rna qc metrics combined.quality.metrics.output object and
#' return in data frame format
#'
#' @param qmo_obj qmo_obj
#'
#' @return A data frame with expanded rna qc metrics
#'
#' @export
get_rna_expanded_metrics <- function(qmo_obj, ...){
  UseMethod("get_rna_expanded_metrics", qmo_obj)
}

#' Get run qc metrics from combined.quality.metrics.output object
#'
#' @param qmo_obj qmo_obj
#' @return A data frame
#' @method get_run_qc_metrics combined.quality.metrics.output
#'
#' @export
get_run_qc_metrics.combined.quality.metrics.output <- function(qmo_obj){
  suppressWarnings(
    if(all(is.na(qmo_obj$run_qc_metrics))){
      run_qc_metrics_df <- data.frame()
    } else {
      run_qc_metrics_df <- qmo_obj$run_qc_metrics %>% 
        select(metric_uom, lsl_guideline, usl_guideline, value)

    }
  )
  return(run_qc_metrics_df)
}

#' Get analysis status from combined.quality.metrics.output object
#'
#' @param qmo_obj qmo_obj
#' @return A data frame
#' @method get_analysis_status combined.quality.metrics.output
#'
#' @export
get_analysis_status.combined.quality.metrics.output <- function(qmo_obj){
  suppressWarnings(
    if(all(is.na(qmo_obj$analysis_status))){
      analysis_status_df <- data.frame()
    } else {
      analysis_status_df <- qmo_obj$analysis_status %>%
        rename(metric = x) %>%
        mutate(across(is.logical, ~as.character(.x))) %>% #otherwise pivot_longer will fail due to logical + character
        pivot_longer(!metric, names_to = "sample_id") %>%
        pivot_wider(names_from = metric, values_from = value)
    }
  )
  return(analysis_status_df)
}

#' Get dna qc metrics from combined.quality.metrics.output object
#'
#' @param qmo_obj qmo_obj
#' @return A data frame
#' @method get_dna_qc_metrics combined.quality.metrics.output
#'
#' @export
get_dna_qc_metrics.combined.quality.metrics.output <- function(qmo_obj){
    suppressWarnings(
        if(all(is.na(qmo_obj$dna_qc_metrics))){
            dna_qc_metrics_df <- data.frame()
        } else {
            dna_qc_metrics_df <- qmo_obj$dna_qc_metrics %>% 
                pivot_longer(!metric_uom, names_to = "sample_id") %>%
                pivot_wider(names_from = metric_uom)
        }
    )
    return(dna_qc_metrics_df)
}

#' Get dna qc metrics for small variants and tmb from combined.quality.metrics.output object
#'
#' @param qmo_obj qmo_obj
#' @return A data frame
#' @method get_dna_qc_metrics_snvtmb combined.quality.metrics.output
#'
#' @export
get_dna_qc_metrics_snvtmb.combined.quality.metrics.output <- function(qmo_obj){
    suppressWarnings(
        if(all(is.na(qmo_obj$dna_qc_metrics_snvtmb))){
            dna_qc_metrics_snvtmb_df <- data.frame()
        } else {
            dna_qc_metrics_snvtmb_df <- qmo_obj$dna_qc_metrics_snvtmb %>%
                pivot_longer(!metric_uom, names_to = "sample_id") %>%
                pivot_wider(names_from = metric_uom)
        }
    )
    return(dna_qc_metrics_snvtmb_df)
}

#' Get dna qc metrics for msi from combined.quality.metrics.output object
#'
#' @param qmo_obj qmo_obj
#' @return A data frame
#' @method get_dna_qc_metrics_msi combined.quality.metrics.output
#'
#' @export
get_dna_qc_metrics_msi.combined.quality.metrics.output <- function(qmo_obj){
    suppressWarnings(
        if(all(is.na(qmo_obj$dna_qc_metrics_msi))){
            dna_qc_metrics_msi_df <- data.frame()
        } else {
            dna_qc_metrics_msi_df <- qmo_obj$dna_qc_metrics_msi %>%
                pivot_longer(!metric_uom, names_to = "sample_id") %>%
                pivot_wider(names_from = metric_uom)
        }
    )
    return(dna_qc_metrics_msi_df)
}

#' Get dna qc metrics for cnv from combined.quality.metrics.output object
#'
#' @param qmo_obj qmo_obj
#' @return A data frame
#' @method get_dna_qc_metrics_cnv combined.quality.metrics.output
#'
#' @export
get_dna_qc_metrics_cnv.combined.quality.metrics.output <- function(qmo_obj){
    suppressWarnings(
        if(all(is.na(qmo_obj$dna_qc_metrics_cnv))){
            dna_qc_metrics_cnv_df <- data.frame()
        } else {
            dna_qc_metrics_cnv_df <- qmo_obj$dna_qc_metrics_cnv %>%
                pivot_longer(!metric_uom, names_to = "sample_id") %>%
                pivot_wider(names_from = metric_uom)
        }
    )
    return(dna_qc_metrics_cnv_df)
}

#' Get expanded dna qc metrics from combined.quality.metrics.output object
#'
#' @param qmo_obj qmo_obj
#' @return A data frame
#' @method get_dna_expanded_metrics combined.quality.metrics.output
#'
#' @export
get_dna_expanded_metrics.combined.quality.metrics.output <- function(qmo_obj){
    suppressWarnings(
        if(all(is.na(qmo_obj$dna_expanded_metrics))){
            dna_expanded_metrics_df <- data.frame()
        } else {
            dna_expanded_metrics_df <- qmo_obj$dna_expanded_metrics %>%
                pivot_longer(!metric_uom, names_to = "sample_id") %>%
                pivot_wider(names_from = metric_uom)
        }
    )
    return(dna_expanded_metrics_df)
}

#' Get rna qc metrics from combined.quality.metrics.output object
#'
#' @param qmo_obj qmo_obj
#' @return A data frame
#' @method get_rna_qc_metrics combined.quality.metrics.output
#'
#' @export
get_rna_qc_metrics.combined.quality.metrics.output <- function(qmo_obj){
    suppressWarnings(
        if(all(is.na(qmo_obj$rna_qc_metrics))){
            rna_qc_metrics_df <- data.frame()
        } else {
            rna_qc_metrics_df <- qmo_obj$rna_qc_metrics %>%
                pivot_longer(!metric_uom, names_to = "sample_id") %>%
                pivot_wider(names_from = metric_uom)
        }
    )
    return(rna_qc_metrics_df)
}

#' Get expanded rna qc metrics from combined.quality.metrics.output object
#'
#' @param qmo_obj qmo_obj
#' @return A data frame
#' @method get_rna_expanded_metrics combined.quality.metrics.output
#'
#' @export
get_rna_expanded_metrics.combined.quality.metrics.output <- function(qmo_obj){
    suppressWarnings(
        if(all(is.na(qmo_obj$rna_expanded_metrics))){
            rna_expanded_metrics_df <- data.frame()
        } else {
            rna_expanded_metrics_df <- qmo_obj$rna_expanded_metrics %>%
                pivot_longer(!metric_uom, names_to = "sample_id") %>%
                pivot_wider(names_from = metric_uom)
        }
    )
    return(rna_expanded_metrics_df)
}

#' Helper function to parse key-value lines in MetricsOutput.tsv
#'
#' @param record_string the record content as string
#'
#' @return char vector
parse_qmo_record <- function(record_string){

  intermediate <- record_string %>%
    trim_qmo_header_and_footer() %>%
    stringr::str_split("\n") %>%
    unlist() %>%
    stringr::str_remove("\\t$") %>%
    stringr::str_split("\\t")

  record <- purrr::map(intermediate, ~ .x[2])
  record_names <- purrr::map_chr(intermediate, ~ .x[1])
  names(record) <- janitor::make_clean_names(record_names)
  return(record)
}

#' Helper function to parse tabular data in MetricsOutput.tsv
#'
#' @param table_string the table as string
#'
#' @return data.frame
parse_qmo_table <- function(table_string){
  
  intermediate <- table_string %>% 
    trim_qmo_header_and_footer()
  
  header_line <- stringr::str_extract(intermediate, ".+\n")

  if(stringr::str_detect(header_line, "\\t\\n")){
    intermediate <- stringr::str_replace_all(
      string = intermediate,
      pattern = "\\t\\n",
      replacement = "\n")
  }

  table_data <- intermediate %>% handle_empty_qmo_table_values()
  return(table_data)
}

#' Helper function to handle empty rows in MetricsOutput.tsv tabular data
#'
#' @param intermediate_tbl intermediate table string
#'
#' @return data.frame
handle_empty_qmo_table_values <- function(intermediate_tbl){
  if(stringr::str_detect(string = intermediate_tbl, pattern = "\\nNA$")){
    return(NA)
  } else {
    to_clean <- utils::read.table(text = intermediate_tbl, sep = "\t", header = TRUE, fill = TRUE)
    clean_name_df <- janitor::clean_names(to_clean)
    return(clean_name_df)
  }
}

#' Helper function to remove header and footer from MetricsOutput.tsv
#'
#' @param string string with file content
#'
#' @return char vector
trim_qmo_header_and_footer <- function(string){
  string %>%
    stringr::str_remove(".+\\t\\t\\n") %>%
    stringr::str_remove_all("[\\t]{2,}") %>%
    stringr::str_remove("[\\n\\t]+$")
}

#' Visualize TSO500 QC results as gt-Table
#'
#' @param qc_df data frame generated with read_xxx_qc_metrics from tso500R package
#' @param id_col column containing the flowcell/run_id name when multiple qc data frames were merged.
#' @param group_name subheader for the subtable containing the actual QC metrics
#'
#' @importFrom dplyr select distinct rename bind_rows mutate case_when
#' @importFrom tidyr pivot_wider replace_na
#' @export
make_qc_table <- function(qc_df, id_col = "sample_id", group_name = "samples") {
  # field names for contamination qc
  p_value_field <- "CONTAMINATION_P_VALUE (NA)"
  contamination_field <- "CONTAMINATION_SCORE (NA)"
  
  # lower QC limit for each metric.
  lsl <- qc_df |>
      filter((!!sym(id_col)) == "lsl_guideline") |>
      mutate(n_failed = NA) |>
      distinct()
  stopifnot(nrow(lsl) == 1)

  # upper QC limit for each metric
  usl <- qc_df |>
      filter((!!sym(id_col)) == "usl_guideline") |>
      mutate(n_failed = 0) |>
      distinct()
  stopifnot(nrow(usl) == 1)

  # QC metrics for each sample
  run_metrics_df <- qc_df |>
      filter(!(!!sym(id_col)) %in% c("lsl_guideline", "usl_guideline"))

  # simple vector with all metrics
  metrics <- colnames(select(run_metrics_df, -!!id_col))

  # Count the numbers of QC failures
  qc_failure_count <- rep(0, nrow(run_metrics_df))
  names(qc_failure_count) <- run_metrics_df[[id_col]]
  for (metric in metrics) {
    # ignore contamination p value since it is only relevant in case of failed 
    # contamination score
    if (metric != p_value_field) {
      if(metric == contamination_field & (p_value_field %in% metrics)) {
        qc_failure_count <- qc_failure_count +
          replace_na(run_metrics_df[[metric]] < lsl[[metric]], 0) +
          replace_na(run_metrics_df[[metric]] > usl[[metric]] & run_metrics_df[[p_value_field]] <= usl[[p_value_field]], 0)
      }
      else { 
        qc_failure_count <- qc_failure_count +
          replace_na(run_metrics_df[[metric]] < lsl[[metric]], 0) +
          replace_na(run_metrics_df[[metric]] > usl[[metric]], 0)
      }
    }
  }

  # Merge tables into one for GT
  merged_table <- bind_rows(
      lsl |> mutate(group = "thresholds"),
      usl |> mutate(group = "thresholds"),
      run_metrics_df |> mutate(group = group_name, n_failed = qc_failure_count)
  )

  table_out <- merged_table |>
      gt::gt(rowname_col = id_col, groupname_col = "group")

  # Conditional formatting
  for (metric in c(metrics, "n_failed")) {
    table_out <- table_out |>
      gt::data_color(
        columns = metric,
        rows = group == group_name,
        fn = \(x) {
          p_value_available <- p_value_field %in% metrics
          if(p_value_available) {
            case_when(
              metric == p_value_field ~ "white",
              (p_value_available & metric == contamination_field & (x > usl[[metric]] & 
                                                                      run_metrics_df[[p_value_field]] <= usl[[p_value_field]])) ~ "#F5CDB9",
              (metric != contamination_field & x < lsl[[metric]]) ~ "#D2F2F7",
              (metric != contamination_field & x > usl[[metric]]) ~ "#F5CDB9",
              .default = "white"
            )
          }
          else {
            case_when(
              x < lsl[[metric]] ~ "#D2F2F7",
              x > usl[[metric]] ~ "#F5CDB9",
              .default = "white"
            )
          }
        }
    )
  }

  # separate the "n_failed" column visually
  table_out <- table_out |>
      gt::tab_style(
          style = gt::cell_borders(sides = "left", weight = "2px", color = "lightgrey"),
          locations = gt::cells_body(columns = "n_failed")
      )

  table_out
}
