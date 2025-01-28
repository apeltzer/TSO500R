#' Plot allele-frequency for small variants per variant consequence
#'
#' @description Plot AF per-sample, grouped by variant consequence
#' @param small_variant_df Data frame containing small variants
#'
#' @return ggplot object
#' 
#' @export
plot_af_per_variant_consequence <- function(small_variant_df){
  g <- ggplot(data = small_variant_df) +
    geom_jitter(aes(x = sample_id, y = allele_frequency, colour = consequence_s),
                size= 0.5, shape = 19, alpha = 0.7) +
    facet_wrap(~consequence_s, ncol = 5) +
    theme(
      panel.grid.major = element_line(color = "grey", linetype = 3, size = 0.1),
      strip.background = element_blank(),
      strip.text = element_text(size=3),
      axis.title = element_text(size=8)
    )
  return(g)
}

#' Plot allele-frequency histogram for small variants
#'
#' @description Plot allele-frequency histogram for small variants, per-sample
#' @param small_variant_df Data frame containing small variants
#'
#' @return ggplot object
#' 
#' @export
plot_af_histogram <- function(small_variant_df){
  g <- ggplot(small_variant_df, aes(x = allele_frequency, group = sample_id,
                                    colour = sample_id, fill = sample_id)) +
    geom_histogram(alpha=0.3, binwidth = 0.025) +
    facet_wrap(~sample_id, ncol = 4) +
    theme(
      panel.grid.major = element_line(color = "grey",
                                      linetype = 3,
                                      size = 0.4))
  return(g)
}

#' Plot allele-frequency kernel density estimate (KDE) for small variants
#'
#' @description Plot allele-frequency KDE for small variants, per-sample
#' @param small_variant_df Data frame containing small variants
#'
#' @return ggplot object
#' @export
plot_af_density <- function(small_variant_df){
  g <- ggplot(small_variant_df, aes(x = allele_frequency, group = sample_id,
                                    colour = sample_id, fill = sample_id)) +
    geom_density(alpha=0.3) +
    facet_wrap(~sample_id, ncol = 4) +
    theme(panel.grid.major = element_line(color = "grey",
                                          linetype = 3,
                                          size = 0.4))
  return(g)
}

#' Plot OncoPrint plot (heatmap) for variants
#'
#' @description Plot OncoPrint plot for specified list of variants
#' @param variant_matrix matrix with variants for samples/entities
#' @param column_title Title of OncoPrint plot
#' @param alter_list List with alter_graphics functions for variant types
#' @param variant_colors List of defined colors for variant types
#' @param heatmap_legend List containing legend specs (title, at, labels)
#' @param bottom_annotation the annotation object for the bottom, default: empty
#' @param top_annotation the annotation object for the top, default: empty
#' @param right_annotation the annotation object for the right, default: empty
#' @param left_annotation the annotation object for the left, default: empty
#' @param top Plot top x genes for given variant matrix (default: all)
#'
#' @return ComplexHeatmap object
#'
#' @export
plot_onco_print <- function(variant_matrix, column_title, alter_list, variant_colors, heatmap_legend, bottom_annotation=ComplexHeatmap::anno_empty(border = FALSE), top_annotation=anno_empty(border = FALSE), right_annotation=anno_empty(border = FALSE), left_annotation=anno_empty(border = FALSE), top=nrow(variant_matrix)){
  # get top x genes
  top_index = order(apply(variant_matrix, 1, function(x) sum(x != "")), decreasing = TRUE)[1:top]
  variant_matrix_top <- variant_matrix[top_index, ]

  onco_print <- ComplexHeatmap::oncoPrint(variant_matrix_top,
          alter_fun = alter_list, col = variant_colors, 
          column_title = column_title,
          show_column_names = TRUE,
          show_row_names = TRUE,
          show_pct = TRUE,
          remove_empty_rows = TRUE,
          heatmap_legend_param = heatmap_legend_param,

          top_annotation = ComplexHeatmap::HeatmapAnnotation(ta = top_annotation, 
            cbar = anno_oncoprint_barplot(height = unit(4, "cm")), show_annotation_name = FALSE),
          bottom_annotation = HeatmapAnnotation(ba = bottom_annotation, show_annotation_name = FALSE),
          right_annotation = rowAnnotation(ra = right_annotation, show_annotation_name = FALSE),
          left_annotation = rowAnnotation(la = left_annotation, show_annotation_name = FALSE)
        )
  return(onco_print)
}

#' Add theme elements for visual customisation
#'
#' @param ggplot_object ggplot2 plot object
#'
#' @return ggplot object
#' 
#' @export
add_common_theme_elements <- function(ggplot_object){
  common_theme_elements <- theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0, size = 5),
    axis.text.y = element_text(size = 5),
    axis.text = element_text(size=5),
    legend.title = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = 'white'),
    panel.border = element_rect(fill= NA, color = "black")
  )
  g2 <- ggplot_object + common_theme_elements
  return(g2)
}
