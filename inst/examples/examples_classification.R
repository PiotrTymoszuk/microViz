# Classification of variables by specificity for a strata of a splitting factor

# tools -------

  library(tidyverse)
  library(rlang)
  library(microViz)
  library(trafo)

  library(pROC)

  select <- dplyr::select
  reduce <- purrr::reduce

# data -------

  data("brca")

  counts <- brca %>%
    mutate(histology = droplevels(histology))

  ## listing all available genes

  genes <- names(brca)

  genes <- genes[!genes %in% c('sample_id',
                               'patient_id',
                               'timepoint',
                               'metastasis',
                               'histology',
                               'er_status')]

  ## a data frame with log2-transformed variables
  ## increasing by 1 to avoid potential log2(zero)

  log_expression <- counts

  log_expression[genes] <- log_expression[genes] %>%
    map_dfc(~log2(.x + 1))

# Classification example --------

  roc_results <-
    classify(log_expression,
             variables = genes,
             split_fct = 'histology')$classification

  ## top markers of the clusters

  top_markers <- roc_results %>%
    group_by(histology) %>%
    filter(auc > 0.8) %>%
    ungroup %>%
    .$variable

  ## heat map representation of Z-scores of expression
  ## levels

  top_heat_map <-
    heat_map(log_expression,
             variables = top_markers,
             split_fct = 'histology',
             limits = c(-5, 5),
             midpoint = 0,
             oob = scales::squish,
             cust_theme = theme_micro() +
               theme(axis.title.y = element_blank(),
                     axis.text.y = element_text(face = 'italic'),
                     axis.text.x = element_blank()))


# END -------

