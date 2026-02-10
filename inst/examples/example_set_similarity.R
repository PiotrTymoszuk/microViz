# Examples of set similarity calculation

  library(tidyverse)
  library(trafo)
  library(microViz)

# analysis data ------

  set.seed(1234)

  ## character sets of differing lengths

  test_sets <-
    list(set_number = 1:40,
         set_size = sample(1:20, 40, replace = TRUE)) %>%
    pmap(function(set_number, set_size) LETTERS %>%
           sample(size = set_size,
                  replace = FALSE)) %>%
    set_names(paste0('obs_', 1:40))

  ## Jaccard

  length(intersect(test_sets[[1]], test_sets[[2]]))/
    length(union(test_sets[[1]], test_sets[[2]]))

  vector_similarity(test_sets[[1]], test_sets[[2]])

  set_similarity(test_sets)

  ## overlap coefficient

  length(intersect(test_sets[[1]], test_sets[[2]]))/
    min(length(test_sets[[1]]), length(test_sets[[2]]))

  vector_similarity(test_sets[[1]], test_sets[[2]], 'overlap')

  set_similarity(test_sets, 'overlap')

  ## Dice coefficient

  2 * length(intersect(test_sets[[1]], test_sets[[2]]))/
    (length(test_sets[[1]]) + length(test_sets[[2]]))

  vector_similarity(test_sets[[1]], test_sets[[2]], 'dice')

  set_similarity(test_sets, 'dice')

  ## Tversky coefficient: for a != b, asymmetric distance matrices are produced

  vector_similarity(test_sets[[1]], test_sets[[2]], 'tversky', a = 0.25, b = 1)

  set_similarity(test_sets, 'tversky', a = 0.25, b = 1)

# Dimensionality reduction and plotting of the set distances -------

  set_dists <- set_similarity(test_sets)

  dist_red <- set_dists %>%
    components(type = 'mds')

  dist_clust <- set_dists %>%
    components(type = 'hcl')

  plot(set_dists,
       type = 'heat_map',
       order_by_similarity = 'ascending')

  plot(set_dists,
       type = 'mds',
       point_wjitter = 0.002,
       point_hjitter = 0.002,
       show_txt = FALSE)

  plot(set_dists,
       type = 'umap',
       point_wjitter = 0.002,
       point_hjitter = 0.002,
       show_txt = FALSE)

  plot(set_dists,
       type = 'dendrogram',
       method = 'ward.D2')

# END -------
