# Example: integration of gene expression by arithmetic mean

# tools and data -------

  library(GEOquery)
  library(tidyverse)

  library(trafo)

  library(org.Hs.eg.db)
  library(AnnotationDbi)

  select <- dplyr::select

# Fetching raw GEO data sets --------

  geo_raw <- getGEO('GSE70769')

# Extraction of the annotation and expression data sets ------

  gse70769_annotation <- fData(geo_raw[[1]]) %>%
    transmute(probe = ID,
              gene = ILMN_Gene) %>%
    filter(complete.cases(.)) %>%
    as_tibble

  gse70769_annotation

  ## roughly 1/4 of probes are duplicated

  gse70769_annotation %>%
    filter(duplicated(gene))

# Extraction of expression and integration -------

  gse70769_expression <- Biobase::exprs(geo_raw[[1]])

  gse70769_expression[1:5, 1:5]

  gse70769_agg_expression <-
    integrate_expression(data = gse70769_expression,
                         annotation = gse70769_annotation,
                         fun = colMeans,
                         transpose = FALSE,
                         trans_fun = identity)

  gse70769_agg_df <-
    integrate_expression(data = gse70769_expression,
                         annotation = gse70769_annotation,
                         fun = colMeans,
                         transpose = TRUE,
                         trans_fun = identity)

# END ------
