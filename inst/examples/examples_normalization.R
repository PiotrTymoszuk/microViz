# Normalization of numeric variables

# tools -------

  library(tidyverse)
  library(rlang)
  library(microViz)
  library(trafo)

  select <- dplyr::select
  reduce <- purrr::reduce

# data -------

  data("brca")

  counts <- brca

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

  ## example genes

  test_genes <- genes[1:100]

# Min-max scaling ---------

  minMax(log_expression$BRCA2)

  minMax(as.matrix(log_expression[test_genes]),
         feature_type = 'columns')

  minMax(as.matrix(log_expression[test_genes]),
         feature_type = 'rows')

  minMax(log_expression, variables = test_genes)

  minMax(log_expression[, test_genes])

# Z-scores ---------

  zScores(as.matrix(log_expression[test_genes]),
          feature_type = 'columns',
          center = 'median')

  zScores(as.matrix(log_expression[test_genes]),
          feature_type = 'rows')

  zScores(log_expression$TSPAN6)

  zScores(log_expression[, test_genes])

  zScores(log_expression,
          variables = test_genes)

  zScores(log_expression,
          variables = test_genes,
          center = 'median')

# END ------
