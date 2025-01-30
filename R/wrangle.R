#' Read in a list of combined.variant.output objects
#' and return a data frame of small variants per sample
#'
#' @param cvo_list a list of combined.variant.output objects
#'
#' @return A data frame of small variants extracted from
#' each combined.variant.output object, per sample
#'
#' @export
#'
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_rows
read_small_variants <- function(cvo_list){

  small_variants <- map(cvo_list, get_small_variants) |>
    bind_rows()

  return(small_variants)
}

#' Read in a list of combined.variant.output objects
#' and return a data frame of gene amplifications per sample
#'
#' @param cvo_list a list of combined.variant.output objects
#'
#' @return A data frame of gene amplifications extracted from
#' each combined.variant.output object, per sample
#'
#' @export
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
read_gene_amplifications <- function(cvo_list) {

  gene_amplifications <- map(cvo_list, get_gene_amplifications) |>
    bind_rows()

  return(gene_amplifications)
}

#' Read in a list of combined.variant.output objects
#' and return a data frame of fusions per sample
#'
#' @param cvo_list a list of combined.variant.output objects
#'
#' @return A data frame of fusions extracted from
#' each combined.variant.output object, per sample
#'
#' @export
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
read_fusions <- function(cvo_list) {

  fusions <- map(cvo_list, get_fusions) |>
    bind_rows()

  return(fusions)
}

#' Read in a list of combined.variant.output objects
#' and return a data frame of splice variants per sample
#'
#' @param cvo_list a list of combined.variant.output objects
#'
#' @return A data frame of splice variants extracted from
#' each combined.variant.output object, per sample
#'
#' @export
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
read_splice_variants <- function(cvo_list) {

  splice_variants <- map(cvo_list, get_splice_variants) |>
    bind_rows()

  return(splice_variants)
}

#' Read in a list of combined.quality.metrics.output objects
#' and return a data frame of run qc metrics per timepoint
#'
#' @param qmo_list a list of combined.quality.metrics.output objects
#'
#' @return A data frame of small variants extracted from
#' each combined.quality.metrics.output object, per timepoint
#'
#' @export
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
read_run_qc_metrics <- function(qmo_list) {

  run_qc_metrics <- map(qmo_list, get_run_qc_metrics) |>
    bind_rows()

  return(run_qc_metrics)
}

#' Read in a list of combined.quality.metrics.output objects
#' and return a data frame of the analysis status per timepoint
#'
#' @param qmo_list a list of combined.quality.metrics.output objects
#'
#' @return A data frame of analysis status extracted from
#' each combined.quality.metrics.output object, per timepoint
#'
#' @export
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
read_analysis_status <- function(qmo_list) {

  analysis_status <- map(qmo_list, get_analysis_status) |>
    bind_rows()

  return(analysis_status)
}

#' Read in a list of combined.quality.metrics.output objects
#' and return a data frame of dna qc metrics per timepoint
#'
#' @param qmo_list a list of combined.quality.metrics.output objects
#'
#' @return A data frame of dna qc metrics extracted from
#' each combined.quality.metrics.output object, per timepoint
#'
#' @export
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
read_dna_qc_metrics <- function(qmo_list) {

  dna_qc_metrics <- map(qmo_list, get_dna_qc_metrics) |>
    bind_rows()

  return(dna_qc_metrics)
}

#' Read in a list of combined.quality.metrics.output objects
#' and return a data frame of dna qc metrics (small variants/tmb) per timepoint
#'
#' @param qmo_list a list of combined.quality.metrics.output objects
#'
#' @return A data frame of dna qc metrics (small variants/tmb) extracted from
#' each combined.quality.metrics.output object, per timepoint
#'
#' @export
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
read_dna_qc_metrics_snvtmb <- function(qmo_list) {

  dna_qc_metrics_snvtmb <- map(qmo_list, get_dna_qc_metrics_snvtmb) |>
    bind_rows()

  return(dna_qc_metrics_snvtmb)
}

#' Read in a list of combined.quality.metrics.output objects
#' and return a data frame of dna qc metrics (msi) per timepoint
#'
#' @param qmo_list a list of combined.quality.metrics.output objects
#'
#' @return A data frame of dna qc metrics (msi) extracted from
#' each combined.quality.metrics.output object, per timepoint
#'
#' @export
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
read_dna_qc_metrics_msi <- function(qmo_list) {

  dna_qc_metrics_msi <- map(qmo_list, get_dna_qc_metrics_msi) |>
    bind_rows()

  return(dna_qc_metrics_msi)
}

#' Read in a list of combined.quality.metrics.output objects
#' and return a data frame of dna qc metrics (cnv) per timepoint
#'
#' @param qmo_list a list of combined.quality.metrics.output objects
#'
#' @return A data frame of dna qc metrics (cnv) extracted from
#' each combined.quality.metrics.output object, per timepoint
#'
#' @export
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
read_dna_qc_metrics_cnv <- function(qmo_list) {

  dna_qc_metrics_cnv <- map(qmo_list, get_dna_qc_metrics_cnv) |>
    bind_rows()

  return(dna_qc_metrics_cnv)
}

#' Read in a list of combined.quality.metrics.output objects
#' and return a data frame of expanded dna qc metrics per timepoint
#'
#' @param qmo_list a list of combined.quality.metrics.output objects
#'
#' @return A data frame of expanded dna qc metrics extracted from
#' each combined.quality.metrics.output object, per timepoint
#'
#' @export
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
read_dna_expanded_metrics <- function(qmo_list) {

  dna_expanded_metrics <- map(qmo_list, get_dna_expanded_metrics) |>
    bind_rows()

  return(dna_expanded_metrics)
}

#' Read in a list of combined.quality.metrics.output objects
#' and return a data frame of rna qc metrics per timepoint
#'
#' @param qmo_list a list of combined.quality.metrics.output objects
#'
#' @return A data frame of rna qc metrics extracted from
#' each combined.quality.metrics.output object, per timepoint
#'
#' @export
#'
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
read_rna_qc_metrics <- function(qmo_list) {

  rna_qc_metrics <- map(qmo_list, get_rna_qc_metrics) |>
    bind_rows()

  return(rna_qc_metrics)
}

#' Read in a list of combined.quality.metrics.output objects
#' and return a data frame of expanded rna qc metrics per timepoint
#'
#' @param qmo_list a list of combined.quality.metrics.output objects
#'
#' @return A data frame of expanded rna qc metrics extracted from
#' each combined.quality.metrics.output object, per timepoint
#'
#' @export
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
read_rna_expanded_metrics <- function(qmo_list) {

  rna_expanded_metrics <- map(qmo_list, get_rna_expanded_metrics) |>
    bind_rows()

  return(rna_expanded_metrics)
}

#' Parse VCF files for a provided path and construct data frame.
#'
#' @param path path to VCF file in `*.vcf` or `*.vcf.gz` format
#'
#' @return {tibble} new data frame with all variants (fixed field and genotype information)
#'
#' @importFrom dplyr mutate left_join
#' @importFrom vcfR read.vcfR vcfR2tidy
#' @importFrom stringr str_split_i
parse_vcf_to_df <- function(path) {
  # parse VCF file
  vcf_content <- read.vcfR(path)
  
  # fixed field content to data frame
  fixed_df <- vcfR2tidy(vcf_content)$fix

  # GT content to data frame
  gt_df <- vcfR2tidy(vcf_content)$gt

  # create addition column with observed nucleotides in order to avoid collisions when we do the left_join
  #gt_df <- gt_df |>
  #  dplyr::mutate(ALT = str_split_i(gt_GT_alleles, "/", 2))
 
  # next use ChromKey, POS and ALT for joining vcf content data frames
  joined_vcf_df <- fixed_df |>
    left_join(gt_df, by = c("ChromKey", "POS"))

  as_tibble(joined_vcf_df)
}

#' Process and filter small variant data-frame to requirements
#'
#' @description Processes small-variant data to comply with requirements for
#' further analysis. The function:
#'
#' \enumerate{
#'   \item filters for variants that:
#'     \itemize{
#'       \item have a consequence in a pre-defined list (see details)
#'       \item are present with depth > 0
#'     }
#'   \item extracts NP ID and protein information from the P Dot-notation column
#'   \item adds columns to faciliate addition of annotation data
#' }
#'
#' @details The following variant consequences are currently included:
#' \itemize{
#'  \item frameshift_variant
#'  \item inframe_deletion
#'  \item inframe_insertion
#'  \item missense_variant
#'  \item missense_variant:splice_region_variant
#'  \item splice_acceptor_variant
#'  \item splice_donor_variant
#'  \item splice_donor_variant:intron_variant
#'  \item start_lost
#'  \item stop_gained
#'  \item stop_gained:splice_region_variant
#'  \item stop_lost
#' }
#'
#' @param small_variant_df Data-frame of small variants
#' 
#' @return processed and filtered data frame
#' 
#' @export
process_and_filter_small_variant_data <- function(small_variant_df) {

  variant_consequences <- c("3_prime_UTR_variant",
                            "5_prime_UTR_variant",
                            "downstream_gene_variant",
                            "intron_variant",
                            "splice_region_variant:intron_variant",
                            "splice_region_variant:synonymous_variant",
                            "synonymous_variant",
                            "upstream_gene_variant")

  updated_df <- small_variant_df |>
    filter_consequences(variant_consequences) |>
    filter_depth(depth_limit = 0) |>
    parse_p_dot_notation() |>
    update_annotation_join_columns()

  return(updated_df)
}

#' Filters variant data for variant consequences. Removes
#' any variant with a consequence matching the list of submitted consequences
#' or an empty consequence field
#'
#' @param variant_df Data frame with variant data
#' @param consequences List of consequences to filter out
#' @param type_column Name of column with variant types (default: consequence_s)
#'
#' @return data frame with filtered consequences
#'
#' @export
#'
#' @importFrom dplyr filter
filter_consequences <- function(variant_df, consequences, type_column = "consequence_s") {
  filtered_df <- variant_df |>
    dplyr::filter(!(get(type_column) %in% consequences | is.na(get(type_column)) | get(type_column) == ""))
  return(filtered_df)
}

#' Filters for variant data for variant consequences. Keeps
#' any variant with a consequence matching the list of submitted consequences 
#'
#' @param variant_df Data frame witgh variant data
#' @param consequences List of consequences to keep
#' @param type_column Name of column with variant types (default: consequence_s)
#'
#' @return data frame with consequences kept as specified
#'
#' @export
#' @importFrom dplyr filter
keep_consequences <- function(variant_df, consequences, type_column = "consequence_s") {
  filtered_df <- variant_df |>
    filter(get(type_column) %in% consequences )
  return(filtered_df)
}

#' Helper function to filter small variant data according
#' to specified depth
#'
#' @param small_variant_df data frame with small variants
#' @param depth_limit depth threshold
#'
#' @return data frame filtered by depth
#'
#' @export
#' 
#' @importFrom dplyr filter
filter_depth <- function(small_variant_df, depth_limit = 0) {
  filtered_df <- filter(small_variant_df, depth > depth_limit)
  return(filtered_df)
}

#' Helper function to filter small variant data according
#' to the GermlineFilterDatabase filter
#'
#' @param small_variant_df data frame with small variants
#'
#' @return data frame filtered by germlineDB annotation
#'
#' @export
#' 
#' @importFrom dplyr filter
filter_germline_db <- function(small_variant_df) {
  filtered_df <- filter(small_variant_df, !GermlineFilterDatabase)
  return(filtered_df)
}

#' Helper function to filter small variant data according
#' to the GermlineFilterProxi filter
#'
#' @param small_variant_df data frame with small variants
#'
#' @return data frame filtered by GermlineFilterProxi
#'
#' @export
#' 
#' @importFrom dplyr filter
filter_germline_proxi <- function(small_variant_df) {
  filtered_df <- filter(small_variant_df, !GermlineFilterProxi)
  return(filtered_df)
}

#' Helper function to filter small variant data and keep
#' only variants with annotated COSMIC ID(s)
#'
#' @param small_variant_df data frame with small variants
#'
#' @return data frame filtered by COSMIC ID
#'
#' @export
#'
#' @importFrom dplyr filter
filter_for_cosmic_id <- function(small_variant_df) {
  filtered_df <- filter(small_variant_df, !is.na(CosmicIDs))
  return(filtered_df)
}

#' Helper function to filter small variant data and keep
#' only variants that are included in TMB numerator
#'
#' @param small_variant_df data frame with small variants
#'
#' @return data frame filtered by TMB inclusion
#'
#' @export
#' 
#' @importFrom dplyr filter
filter_for_included_in_tmb <- function(small_variant_df) {
  filtered_df <- filter(small_variant_df, IncludedInTMBNumerator)
  return(filtered_df)
}

#' Parses P-Dot notation column, splitting it into distinct NP ID and
#' amino acid variant columns
#'
#' @param small_variant_df data frame with small variants
#'
#' @return parsed protein notation
#' 
#' @export
#'
#' @importFrom tidyr extract
parse_p_dot_notation <- function(small_variant_df) {
  extract(
    data = small_variant_df,
    col = p_dot_notation,
    into = c("np_id", "protein"),
    regex = "(NP_.+):.+\\((\\w{3}\\d+).+$",
    remove = FALSE
  )
}

#' Helper function to make small variant DF joinable to annotation data
#'
#' @param small_variant_df data frame with small variants
#'
#' @return data frame
#' 
#' @importFrom dplyr mutate
update_annotation_join_columns <- function(small_variant_df) {
  mutated_df <- small_variant_df |>
    mutate(
      protein = paste(gene, protein, sep = "_"),
      coord_id = paste(chromosome, genomic_position, sep = "_")
    )
  return(mutated_df)
}

#' Adds GEL annotation data to small variant data
#'
#' @param small_variant_df data frame with small variants
#' @param annotation_data_list annotation data list
#'
#' @return data frame with annotation data
#'
#' @export
#'
#' @importFrom purrr reduce2
add_annotation_data <- function(small_variant_df, annotation_data_list) {
  key_order <- c("protein", "coord_id", "gene", "gene", "gene", "gene", "gene")
  to_join <- c(list(small_variant_df), annotation_data_list)
  joined_data <- reduce2(.x = to_join, .y = key_order, .f = left_join)
  return(joined_data)
}

#' Adds information from TMB trace table to small variant data
#'
#' @param small_variant_df data frame with small variants
#' @param tmb_variant_df data frame with TMB trace information
#'
#' @return joined data frame
#'
#' @export
#' 
#' @importFrom dplyr left_join
add_tmb_variant_data <- function(small_variant_df, tmb_variant_df) {
  joined_data <- small_variant_df |>
    left_join(tmb_variant_df, by = c("sample_id", "chromosome" = "Chromosome", "genomic_position" = "Position", "refere
      nce_call" = "RefCall", "alternative_call" = "AltCall"))
  return(joined_data)
}

#' Adds amplifications to small variant data, 
#' variant type (consequence_s) column renamed to variant_type
#'
#' @param small_variant_df data frame with small variants
#' @param amplification_df data frame with amplifications
#'
#' @return joined data frame
#'
#' @export
#'
#' @importFrom dplyr mutate select rename bind_rows if_else
add_amplification_data <- function(small_variant_df, amplification_df) {
  prepared_amplification_df <- amplification_df |>
    mutate(variant_type = if_else(fold_change < 1.0, "DEL", "DUP")) |>
    select(-fold_change)

  joined_data <- small_variant_df |> 
    rename(variant_type = consequence_s) |>
    bind_rows(prepared_amplification_df)
  return(joined_data)
}

#' Helper function to extract TMB/MSI metrics from
#' list of combined.variant.output objects
#'
#' @param cvo_data CVO object
#' @param category category to extract, default: tmb
#'
#' @return data.frame
#'
#' @importFrom purrr map_dfr
extract_metrics <- function(cvo_data, category = "tmb") {
  df <- map_dfr(cvo_data, purrr::pluck, category)
  return(df)
}

#' Extracts all TMB/MSI metrics from list of combined.variant.output
#' objects, returning a data frame of TMB/MSI per sample
#'
#' @param cvo_data CVO object
#'
#' @return data frame with TMB/MSI metrics
#'
#' @export
#'
#' @importFrom dplyr bind_cols mutate select
#' @importForm purrr map_chr
#' @importFrom tidyr everything
get_metrics_df <- function(cvo_data) {
  tmb_df <- extract_metrics(cvo_data, category = "tmb")
  msi_df <- extract_metrics(cvo_data, category = "msi")

  metrics_df <- bind_cols(tmb_df, msi_df) |>
    mutate(sample_id = map_chr(cvo_data, ~ ifelse(is.null(.x$analysis_details$pair_id),
                                                          .x$analysis_details$dna_sample_id,
                                                          .x$analysis_details$pair_id))) |>
    select(sample_id, everything())
  return(metrics_df)
}

#' Extracts all analysis details from list of combined.variant.output
#' objects, returning a data frame of analysis details per sample
#'
#' @param cvo_data CVO object
#'
#' @return data frame with analysis details
#'
#' @export
#'
#' @importFrom dplyr mutate select
#' @importFrom purrr map_chr
#' @importFrom tidyr everything
get_analysis_details_df <- function(cvo_data) {
  analysis_details <- extract_metrics(cvo_data, category = "analysis_details")
  analysis_details_df <- analysis_details |>
    mutate(sample_id = map_chr(cvo_data, ~ ifelse(is.null(.x$analysis_details$pair_id),
                                                          .x$analysis_details$dna_sample_id,
                                                          .x$analysis_details$pair_id))) |>
    select(sample_id, everything())
  return(analysis_details_df)
}

#' Extracts all sequencing run details from list of combined.variant.output
#' objects, returning a data frame of sequencing run details per sample
#'
#' @param cvo_data CVO object
#'
#' @return data frame with sequencing run details
#'
#' @export
#'
#' @importFrom dplyr mutate select
#' @importFrom purrr map_chr
#' @importFrom tidyr everything
get_sequencing_run_details_df <- function(cvo_data) {
  sequencing_run_details <- extract_metrics(cvo_data, category = "sequencing_run_details")
  sequencing_run_details_df <- sequencing_run_details |>
    mutate(sample_id = map_chr(cvo_data, ~ ifelse(is.null(.x$analysis_details$pair_id),
                                                          .x$analysis_details$dna_sample_id,
                                                          .x$analysis_details$pair_id))) |>
    select(sample_id, everything())
  return(sequencing_run_details_df)
}

#' Helper function to extract summarized counts based on sample_id
#'
#' @param data_df data frame
#' @param column_name column name
#'
#' @return data.frame
#'
#' @importFrom dplyr group_by summarise select
#' @importFrom tibble add_column
get_summarised_statistics_df <- function(data_df, column_name) {
  if (!(nrow(data_df) == 0)) {
    summarized_data_df <- data_df |>
      group_by(sample_id) |>
      summarise({{column_name}} := n()) |>
      select(sample_id, {{column_name}})
  } else {
    summarized_data_df <- data_df |>
      add_column(sample_id = "") |>
      add_column({{column_name}} := "")
  }
  return(summarized_data_df)
}

#' Extracts the counts of the different variant types from list of combined.variant.output
#' objects, returning a data frame of counts per sample
#'
#' @param cvo_data CVO object
#'
#' @return data frame with variant type statistics
#'
#' @export
#'
#' @importFrom dplyr mutate select full_join
#' @importFrom purrr map_chr
#' @importFrom tibble add_column
get_count_df <- function(cvo_data) {
  # most of the following steps are taken to make sure that we get numbers for all samples and
  # variant types
  default_df <- extract_metrics(cvo_data, category = "sequencing_run_details") |>
    mutate(sample_id = map_chr(cvo_data, ~ ifelse(is.null(.x$analysis_details$pair_id),
                                                          .x$analysis_details$dna_sample_id,
                                                          .x$analysis_details$pair_id))) |>
    select(sample_id) |>
    add_column(number_of_amplifications = NA) |>
    add_column(number_of_small_variants = NA) |>
    add_column(number_of_fusions = NA) |>
    add_column(number_of_splice_variants = NA)

  amps <- read_gene_amplifications(cvo_data) |>
    get_summarised_statistics_df("number_of_amplifications")
  small_variants <- read_small_variants(cvo_data) |>
    get_summarised_statistics_df("number_of_small_variants")
  fusions <- read_fusions(cvo_data) |>
    get_summarised_statistics_df("number_of_fusions")
  splice_variants <- read_splice_variants(cvo_data) |>
    get_summarised_statistics_df("number_of_splice_variants")

  counts_df <- default_df |> 
    full_join(amps, by = 'sample_id') |>
    mutate(number_of_amplifications=coalesce(number_of_amplifications.x, number_of_amplifications.y)) |>
    full_join(small_variants, by = 'sample_id') |>
    mutate(number_of_small_variants=coalesce(number_of_small_variants.x,number_of_small_variants.y)) |>
    full_join(fusions, by = 'sample_id') |>
    mutate(number_of_fusions=coalesce(number_of_fusions.x,number_of_fusions.y)) |>
    full_join(splice_variants, by = 'sample_id') |>
    mutate(number_of_splice_variants=coalesce(number_of_splice_variants.x,number_of_splice_variants.y)) |>
    select(sample_id, number_of_amplifications, number_of_small_variants, number_of_fusions, number_of_splice_variants) |>
    mutate(number_of_amplifications = ifelse(is.na(number_of_amplifications), 0, number_of_amplifications),
          number_of_fusions = ifelse(is.na(number_of_fusions), 0, number_of_fusions),
          number_of_small_variants = ifelse(is.na(number_of_small_variants), 0, number_of_small_variants),
          number_of_splice_variants = ifelse(is.na(number_of_splice_variants), 0, number_of_splice_variants))

  return(counts_df)
}

#' Transforms data frame holding variant information to a matrix that can be used as OncoPrint input
#'
#' @param variant_data_frame Data frame holding variant information
#' @param id_column Column holding identifiers
#' @param gene_column Column holding genes
#' @param variant_type_column Column holding variant types
#'
#' @return data frame as needed for onocprint plot
#'
#' @export
#'
#' @importFrom dplyr pivor_widder
prepare_dataframe_for_oncoprint <- function(variant_data_frame, id_column="sample_id", gene_column="gene", variant_type_column="consequence_s") {
  oncoprint_df <- variant_data_frame |>
    pivot_wider(id_cols = id_column, names_from = gene_column, values_from = variant_type_column, values_fn = function(x) paste(x, collapse=";"))

  oncoprint_matrix <- as.matrix(oncoprint_df)
  oncoprint_matrix[is.na(oncoprint_matrix)] <- ""
  rownames(oncoprint_matrix) <- oncoprint_matrix[, 1]
  oncoprint_matrix <- oncoprint_matrix[, -1]
  oncoprint_matrix <- t(as.matrix(oncoprint_matrix))

  return(oncoprint_matrix)
}

#' Parses DRAGEN sample sheet and returns a dataframe
#'
#' @param samplesheet path to Illumina sample sheet
#' @param id_column Column holding identifiers
#' @param gene_column Column holding genes
#' @param variant_type_column Column holding variant types
#'
#' @return names list of data frames (header, settings, data) holding sample sheet content
#'
#' @export
#'
#' @importFrom readr read_file read_csv
#' @importFrom stringr str_split
#' @importFrom purrr discard
#' @importFrom tidyr pivot_wider
parse_illumina_samplesheet <- function(samplesheet) {
  # illumina sample sheet sections
  HEADER_STRING <- "[Header]"
  CHEMISTRY_STRING <- "Chemistry"
  SETTINGS_STRING <- "[Settings]"
  OVERRIDE_CYCLES_STRING <- "OverrideCycles"
  DATA_STRING <- "[Data]"

  # Read text file into string
  samplesheet_file <- readr::read_file(samplesheet)
  split_samplesheet_string <- stringr::str_split(string = samplesheet_file, pattern = "\r\n") |>
    unlist()

  # parse header part of provided sample sheet
  start_header <- pmatch(HEADER_STRING, split_samplesheet_string) + 1
  end_header <- which(grepl(CHEMISTRY_STRING, split_samplesheet_string))
  header <- read_csv(I(split_samplesheet_string[start_header:end_header]), col_names = FALSE) |>
    discard(~all(is.na(.))) |>
    pivot_wider(names_from = X1, values_from = X2)

  # parse settings part of provided sample sheet
  start_settings <- pmatch(SETTINGS_STRING, split_samplesheet_string) + 1
  end_settings <- which(grepl(OVERRIDE_CYCLES_STRING, split_samplesheet_string))
  settings <- read_csv(I(split_samplesheet_string[start_settings:end_settings]), col_names = FALSE) |>
    discard(~all(is.na(.))) |>
    pivot_wider(names_from = X1, values_from = X2)

  # prase data part of provided sample sheet
  start_data <- pmatch(DATA_STRING, split_samplesheet_string) + 1
  data <- read_csv(I(split_samplesheet_string[start_data:length(split_samplesheet_string)]))

  return(list(header = header, settings = settings, data = data))
}