# Oncoplot for binary data

# tools -------

  library(tidyverse)
  library(rlang)
  library(microViz)
  library(trafo)

  library(perich)

  select <- dplyr::select
  reduce <- purrr::reduce

# data -------

  ## 0/ coded mutation data in three molecular clusters of
  ## bladder cancers

  data("tcga_mutations")

  ## the most common genes present in at least 40 samples

  mut_frequency <- tcga_mutations %>%
    select(-sample_id, -clust_id) %>%
    colSums

  genes <- names(mut_frequency)[mut_frequency >= 40]

  ## classification of the genes

  gene_classification <-
    tibble(gene_symbol = genes,
           subset = 'other') %>%
    mutate(subset = ifelse(gene_symbol %in% c('ATM', 'BRCA2', 'TP53'),
                           'repair',
                           ifelse(gene_symbol %in% c('AHNAK', 'AHNAK2',
                                                     'DNAH11', 'DNAH3', 'DNAH5',
                                                     'DNAH8', 'FAT1', 'FAT3', 'FAT4',
                                                     'FLG', 'GPR98', 'MACF1',
                                                     'MUC16', 'MUC17', 'NEB',
                                                     'OBSCN', 'PCLO', 'PDE4DIP',
                                                     'SPTA1', 'SPTAN1', 'STAG2',
                                                     'SYNE1', 'SYNE2', 'TTN', 'USH2A',
                                                     'XIRP2'),
                                  'adhesion/cytoskeleton',
                                  ifelse(gene_symbol %in% c('AKAP9', 'ELF3',
                                                            'ERBB2', 'ERBB3',
                                                            'FGFR3', 'HMCN1',
                                                            'PIK3CA'),
                                         'signaling',
                                         ifelse(gene_symbol %in% c('ARID1A', 'EP300',
                                                                   'KDM6A', 'KMT2A', 'KMT2C',
                                                                   'KMT2D'),
                                                'epigenetics',
                                                ifelse(gene_symbol %in% c('BIRC6', 'RB1'),
                                                       'cell cycle', subset))))),
           subset = factor(subset))


# plots --------

  ## onco plot panel

  onco_plot <-
    plot_bionco(data = tcga_mutations,
                variables = genes,
                split_fct = 'clust_id',
                #variable_classification = gene_classification,
                plot_title = 'Mutation landscape of BLCA molecular clusters',
                hide_x_axis_text = TRUE,
                x_lab = 'TCGA BLCA cancer sample',
                y_lab = NULL,
                y_text_face = 'italic',
                one_plot = TRUE,
                cust_theme = theme_micro(),
                legend_position = 'right',
                y_rug_scale = 'percentage')

  ## stack plot

  stack_plot <-
    plot_bistack(data = tcga_mutations,
                 variables = genes[1:10],
                 split_fct = 'clust_id',
                 variable_classification = gene_classification,
                 scale = 'percentage',
                 show_frequencies = TRUE)


# END -------
