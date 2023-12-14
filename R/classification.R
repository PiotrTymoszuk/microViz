# Classification of selected features by their specificity for data subsets
# and similarity between data subsets.

# Classification ---------

#' Classify variables by their specificity for data set subsets.
#'
#' @description
#' Classifies selected features, e.g. differentially expressed genes, by their
#' specificity for the data subsets defined by a splitting factor.
#'
#' @details
#' Specificity for a given data subset is determined by calculation of area
#' under the receiver-operating characteristic (ROC) in a comparison of the
#' subset with the pooled remaining subsets (i.e. the rest of the data set).
#' Of note, these remaining subsets are set as the baseline or control in the
#' ROC analysis. AUC ROC is computed with \code{\link[pROC]{roc}}.
#'
#' @param data a data frame.
#' @param variables a vector of names of the variables of interest.
#' @param split_fct splitting factor defining the data set subsets.
#' @param direction direction of the comparison in the ROC analysis.
#' See: \code{\link[pROC]{roc}} for details.
#'
#' @return
#' A list with two elements, `auc` and `classification`.
#' The `auc` element is a data frame that stores the complete analysis results.
#' `classification` lists variables along with their most specific subsets of
#' the data set. Both elements are data frames with the following columns:
#'
#' * `variable`: variable name
#' * `split_fct`: a column named after the splitting factor
#' * `auc`: are under the ROC curve for detection of the particular data subset
#' * `delta_auc`: difference in AUC ROC for the particular subset and mean AUC
#' ROC computed across all subsets.
#'
#' @export

  classify <- function(data,
                       variables,
                       split_fct,
                       direction = c('<', '>', 'auto')) {

    ## classifies variables for their ROC AUC at detection of the subsets

    ## input control --------

    if(!is.data.frame(data)) {

      stop("'data' hs to be a data frame.", call. = FALSE)

    }

    if(any(!variables %in% names(data))) {

      stop('Some variables missing from the input data frame.', call. = FALSE)

    }

    if(!split_fct %in% names(data)) {

      stop('The splitting factor is missing from the input data frame.',
           call. = FALSE)

    }

    if(!is.factor(data[[split_fct]])) {

      stop("'split_fct' argument has to be a factor.",
           call. = FALSE)

    }

    direction <-
      match.arg(direction[1],
                c('<', '>', 'auto'))

    ## analysis data ------

    ## subsets are compared in a pairwise manner

    data <- data[c(split_fct, variables)]

    split_levs <- levels(data[[split_fct]])

    split_levs <- set_names(split_levs, split_levs)

    sub_class <- NULL

    split_data <-
      map(split_levs,
          ~mutate(data,
                  sub_class = ifelse(.data[[split_fct]] == .x,
                                     'yes', 'no'),
                  sub_class = factor(sub_class, c('no', 'yes'))))

    ## ROC analysis, extraction of the AUC ------

    ## ROC AUC for detection of the particular subset versus the pooled rest

    roc_objects <-
      map(split_data,
          function(df) map(variables,
                           ~pROC::roc(df[['sub_class']],
                                      df[[.x]],
                                      levels = c('no', 'yes'),
                                      direction = direction)))

    roc_objects <- transpose(roc_objects)

    roc_objects <- set_names(roc_objects, variables)

    auc_tbl <- map(roc_objects, map_dbl, ~.x$auc)

    auc <- NULL
    delta_auc <- NULL
    variable <- NULL

    auc_tbl <-
      map(auc_tbl,
          ~tibble(!!split_fct := names(.x),
                  auc = .x))

    auc_tbl <-
      map2_dfr(auc_tbl, names(auc_tbl),
               ~mutate(.x, variable = .y))

    ## computing deviance from the mean AUC

    auc_tbl <- group_by(auc_tbl, variable)

    auc_tbl <- mutate(auc_tbl, delta_auc = auc - mean(auc))

    auc_tbl <- ungroup(auc_tbl)

    auc_tbl <-
      mutate(auc_tbl,
             !!split_fct := factor(.data[[split_fct]],
                                   split_levs))

    ## selecting the cluster with the best detection accuracy ------

    best_clusters <- group_by(auc_tbl, variable)

    best_clusters <- filter(best_clusters, auc == max(auc))

    best_clusters <- ungroup(best_clusters)

    best_clusters <- arrange(best_clusters, .data[[split_fct]], delta_auc)

    list(auc = auc_tbl[c('variable',
                         split_fct,
                         'auc',
                         'delta_auc')],
         classification = best_clusters[c('variable',
                                          split_fct,
                                          'auc',
                                          'delta_auc')])

  }

# Similarity and distance --------

#' Similarity and distance between data subsets.
#'
#' @description
#' The function `subset_distance()` computes distances such as Euclidean,
#' Manhattan or cosine distance, or similarity measures between subsets of
#' a data set in respect to variables of interest.
#'
#' @details
#' Since there are multiple ways to calculate distances in R, the function
#' allows the user to work with their favorite one - provided it returns a
#' numeric matrix or a `dist` class instance. We can recommend
#' \code{\link[stats]{dist}} and \code{\link[clustTools]{calculate_dist}}.
#' The function returns a list of numeric distance matrices of
#' the \code{\link[clustTools]{cross_dist}} class, with multiple analytic and
#' visualization tools provided by the
#' [clustTools](https://github.com/PiotrTymoszuk/clustTools) package.
#'
#' @param data a data frame.
#' @param variables a vector of names of the variables of interest.
#' @param split_fct name of the splitting factor, the variable defining the
#' data set subsets.
#' @param normalize logical, should the data frame variables be normalized prior
#' to plotting?
#' @param norm_center defines centering of the variable during scaling: mean
#' (default) or median. Ignored if `normalize = FALSE`.
#' @param dist_FUN a distance-computing function.
#' @param ... additional arguments passed to `dist_FUN`.
#'
#' @export

  subset_distance <- function(data,
                              variables,
                              split_fct,
                              normalize = TRUE,
                              norm_center = c('mean', 'median'),
                              dist_FUN = stats::dist, ...) {

    ## input control -------

    stopifnot(is.data.frame(data))

    if(!is.data.frame(data)) {

      stop("'data' hs to be a data frame.", call. = FALSE)

    }

    if(any(!variables %in% names(data))) {

      stop('Some variables missing from the input data frame.', call. = FALSE)

    }

    if(!split_fct %in% names(data)) {

      stop('The splitting factor is missing from the input data frame.',
           call. = FALSE)

    }

    if(!is.factor(data[[split_fct]])) {

      stop("'split_fct' argument has to be a factor.",
           call. = FALSE)

    }

    if(!is.function(dist_FUN)) {

      stop("'dist_FUN' has to be a function.", call. = FALSE)

    }

    stopifnot(is.logical(normalize))

    norm_center <- match.arg(norm_center,
                             c('mean', 'median'))

    ## data pre-processing --------

    data <- data[c(split_fct, variables)]

    if(normalize) {

      center_fun <- switch(norm_center,
                           mean = function(x) mean(x, na.rm = TRUE),
                           median = function(x) median(x, na.rm = TRUE))

      data[variables] <-
        map_dfc(data[variables],
                ~scale(.x, center = center_fun(.x))[, 1])

    }

    data <- filter(data, complete.cases(data))

    data_mtx <- as.matrix(data[variables])

    if(is.null(rownames(data_mtx))) {

      observations <- paste0('obs_', 1:nrow(data_mtx))

    } else {

      observations <- rownames(data_mtx)

    }

    observation <- NULL

    assignment <- tibble(observation = observations,
                         !!split_fct := data[[split_fct]])

    ## distance calculation --------

    dist_mtx <- dist_FUN(data_mtx, ...)

    err_txt <- "'dist_FUN' has to return a 'dist' or a numeric matrix object."

    if(!is.matrix(dist_mtx)) {

      if(!inherits(dist_mtx, 'dist')) stop(err_txt, call. = FALSE)

    }

    dist_mtx <- as.matrix(dist_mtx)

    if(!is.numeric(dist_mtx)) stop(err_txt, call. = FALSE)

    if(nrow(dist_mtx) != length(observations) |
       ncol(dist_mtx) != length(observations)) {

      stop('Improper format of the distance matrix.',
           call. = FALSE)

    }

    rownames(dist_mtx) <- observations
    colnames(dist_mtx) <- observations

    ## pairwise distances between the clusters ------

    split_levs <- levels(data[[split_fct]])

    observation_map <- split(assignment$observation,
                             assignment[[split_fct]])

    subset_identities <- map(split_levs, ~c(.x, .x))

    subset_pairs <- utils::combn(split_levs, m = 2, simplify = FALSE)

    subset_pairs <- c(subset_pairs, subset_identities)

    subset_pair_names <-
      map(subset_pairs,
          paste,
          collapse = ' vs ')

    subset_pairs <- set_names(subset_pairs, subset_pair_names)

    subset_dists <-
      map(subset_pairs,
          ~dist_mtx[observation_map[[.x[1]]],
                    observation_map[[.x[2]]]])

    ## cross-distance object -------

    dots <- rlang::list2(...)

    if(!is.null(dots$method)) {

      distance_method <- dots$method

    } else {

      distance_method <- 'euclidean'

    }

    subset_dists <- structure(subset_dists,
                              class = c('cross_dist', class(subset_dists)))

    attr(subset_dists, 'type') <- 'homologous'
    attr(subset_dists, 'dist_method') <- distance_method
    attr(subset_dists, 'x_levels') <- split_levs
    attr(subset_dists, 'y_levels') <- split_levs

    subset_dists

  }

# END -----
