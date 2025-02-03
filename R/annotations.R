#' Read Excel spreadsheet into a list of data frames, per sheet
#'
#' @param annotation_data_path Path to annotation .xlsx
#' @param sheet_names List of names of the sheets in the Excel spreadsheet
#'
#' @return list of data.frame
#' @export
read_annotation_data <- function(annotation_data_path, sheet_names){
  annotations <- map(sheet_names, ~ readxl::read_excel(annotation_data_path, sheet = .x) |> janitor::clean_names())
  names(annotations) <- sheet_names
  return(annotations)
}
