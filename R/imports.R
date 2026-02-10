# Imports of frequently used functions from dependencies

#' @importFrom rlang `.data`
#' @importFrom rlang `!!`
#' @importFrom rlang `:=`
#' @importFrom rlang set_names
#' @importFrom rlang list2
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 scale_size_area
#' @importFrom ggplot2 scale_radius
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_color_gradient2
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 geom_errorbarh
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 position_jitter
#' @importFrom ggplot2 position_dodge
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 guide_axis
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 margin
#' @importFrom ggplot2 geom_label
#'
#' @importFrom ggwordcloud geom_text_wordcloud
#' @importFrom ggwordcloud geom_text_wordcloud
#'
#' @importFrom cowplot plot_grid
#' @importFrom cowplot get_legend
#' @importFrom cowplot get_title
#' @importFrom cowplot get_subtitle
#'
#' @importFrom purrr map
#' @importFrom purrr map_lgl
#' @importFrom purrr map2
#' @importFrom purrr transpose
#' @importFrom purrr pmap
#' @importFrom purrr map_dbl
#' @importFrom purrr map_dfr
#' @importFrom purrr map2_dfr
#' @importFrom purrr map_dfc
#' @importFrom purrr compact
#' @importFrom purrr map_chr
#' @importFrom purrr map2_chr
#' @importFrom purrr reduce
#' @importFrom purrr pluck
#'
#' @importFrom dplyr tibble
#' @importFrom dplyr as_tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr full_join
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr transmute
#' @importFrom dplyr count
#' @importFrom dplyr top_n
#' @importFrom dplyr summarise
#' @importFrom dplyr relocate
#' @importFrom dplyr arrange
#' @importFrom dplyr all_of
#' @importFrom dplyr pick
#' @importFrom dplyr desc
#'
#' @importFrom tidyr pivot_longer
#'
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames
#' @importFrom tibble is_tibble
#'
#' @importFrom stats median
#' @importFrom stats complete.cases
#' @importFrom stats p.adjust
#' @importFrom stats var
#' @importFrom stats as.formula
#' @importFrom stats dist
#' @importFrom stats as.dist
#' @importFrom stats confint
#' @importFrom stats anova
#' @importFrom stats nobs
#' @importFrom stats reorder
#' @importFrom stats quantile
#' @importFrom stats t.test
#' @importFrom stats aov
#' @importFrom stats summary.lm
#' @importFrom stats rnorm
#' @importFrom stats cmdscale
#' @importFrom stats hclust
#' @importFrom stats na.omit
#' @importFrom stats pt
#' @importFrom stats qt
#'
#' @importFrom umap umap
#'
#' @importFrom utils combn
#'
#' @importFrom car recode
#'
#' @importFrom stringi stri_detect
#' @importFrom stringi stri_replace
#' @importFrom stringi stri_split
#'
#' @importFrom future plan
#'
#' @importFrom furrr future_map2
#' @importFrom furrr furrr_options
#' @importFrom furrr future_map
#' @importFrom furrr future_map_dfc
#'
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggrepel geom_label_repel
#'
#' @importFrom ggforce geom_ellipse
#'
#' @importFrom limma goana
#'
#' @importFrom GOSemSim geneSim
#' @importFrom GOSemSim goSim
#'
#' @importFrom eulerr euler
#'
#' @importFrom ggdendro ggdendrogram
#'
#' @importFrom generics components
#'
#' @importFrom Rcpp evalCpp
#'
#' @useDynLib microViz

  NULL
