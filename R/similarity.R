# Expression and semantic similarity between genes or GO terms

#' @include imports.R

  NULL

# Semantic similarity ---------

#' Semantic similarity between GO terms or genes.
#'
#' @description
#' The `gen_sem()` and `GO_sem()` functions compute semantic similarities
#' between genes and gene ontology (GO) terms using the interface provided
#' by the `GOSemSim` package. Given a vector of gene Entrez ID or GO IDs,
#' a similarity or distance matrix is returned.
#'
#' @details
#' Please refer to the
#' GOSemSim package manual, and to the downstream functions
#' \code{\link[GOSemSim]{geneSim}} and \code{\link[GOSemSim]{goSim}}.
#' The parallel backend is provided by the `furrr`
#' package tools and may not work with some platforms. Distances are computed
#' with the `distance = 1 - similarity` formula.
#'
#' @references
#' Yu G, Li F, Qin Y, Bo X, Wu Y, Wang S. GOSemSim: an R package for measuring
#' semantic similarity among GO terms and gene products. Bioinformatics (2010)
#' 26:976–978. doi:10.1093/BIOINFORMATICS/BTQ064
#'
#' @return a numeric matrix or an object of the `dist` class containing pairwise
#' distances or similarities between the genes or GO terms-
#'
#' @param x a vector of Entrez IDs, gene symbols or genes or GO IDs, specific
#' for the `semData` argument.
#' @param semData a GOSemSimDATA object, created by e.g.
#' \code{\link[GOSemSim]{godata}}.
#' @param as_dist logical, if `TRUE` an object of the `dist` class is returned.
#' @param pad_na `NULL` or a numeric value used to replace any NAs in the output
#' matrix. If `NULL`, NAs are left as they are.
#' @param .parallel logical, should the analysis be run in parallel?
#' @param measure distance measure, please consult
#' \code{\link[GOSemSim]{geneSim}} or \code{\link[GOSemSim]{goSim}}.
#' @param drop evidence rules for dropping annotations, please consult
#' \code{\link[GOSemSim]{geneSim}} or \code{\link[GOSemSim]{goSim}}.
#' @param combine a method for combining semantic similarity scores, please
#' consult \code{\link[GOSemSim]{geneSim}} or \code{\link[GOSemSim]{goSim}}
#'
#' @export

  gene_sem <- function(x,
                       semData,
                       measure = 'Wang',
                       drop = 'IEA',
                       combine = 'BMA',
                       as_dist = TRUE,
                       pad_na = NULL,
                       .parallel = FALSE) {

    stopifnot(is.character(x))
    stopifnot(is.logical(as_dist))
    stopifnot(is.logical(.parallel))

    if(!is.null(pad_na)) stopifnot(is.numeric(pad_na))

    pad_na <- pad_na[1]

    start_time <- Sys.time()

    if(length(x) < 2) stop('At least two features are required.', call. = FALSE)

    message(paste('Measuring similarity for', length(x), 'features'))

    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))
    on.exit(plan('sequential'), add = TRUE)

    if(.parallel) plan('multisession')

    x <- unname(x)

    x <- x[!is.na(x)]

    ind_combi <- combn(seq_along(x), m = 2, simplify = FALSE)
    x_combi <- combn(x, m = 2, simplify = FALSE)

    sim_lst <- future_map(x_combi,
                          ~geneSim(.x[1],
                                   .x[2],
                                   semData = semData,
                                   measure = measure,
                                   drop = drop,
                                   combine = combine),
                          .options = furrr_options(seed = TRUE))

    sim_lst <- map_dbl(sim_lst, ~.x[[1]])

    ind_mat <- map(ind_combi, reduce, cbind)

    ind_mat <- do.call('rbind', ind_mat)

    res <- matrix(1, nrow = length(x), ncol = length(x))

    for(i in seq_along(sim_lst)) {

      res[ind_mat[i, 1], ind_mat[i, 2]] <- sim_lst[i]
      res[ind_mat[i, 2], ind_mat[i, 1]] <- sim_lst[i]

    }

    rownames(res) <- x
    colnames(res) <- x

    if(!is.null(pad_na)) {

      res <- ifelse(is.na(res), pad_na, res)

    }

    if(as_dist) res <- as.dist(1 - res)

    res

  }

#' @rdname gene_sem
#' @export

  go_sem <- function(x, semData,
                     measure = 'Wang',
                     as_dist = TRUE,
                     pad_na = NULL,
                     .parallel = FALSE) {

    stopifnot(is.character(x))
    stopifnot(is.logical(as_dist))
    stopifnot(is.logical(.parallel))

    if(!is.null(pad_na)) stopifnot(is.numeric(pad_na))

    pad_na <- pad_na[1]

    start_time <- Sys.time()

    if(length(x) < 2) stop('At least two features are required.', call. = FALSE)

    message(paste('Measuring similarity for', length(x), 'features'))

    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))
    on.exit(plan('sequential'), add = TRUE)

    if(.parallel) plan('multisession')

    x <- unname(x)

    x <- x[!is.na(x)]

    ind_combi <- combn(seq_along(x), m = 2, simplify = FALSE)
    x_combi <- combn(x, m = 2, simplify = FALSE)

    sim_lst <- future_map(x_combi,
                          ~goSim(.x[1],
                                 .x[2],
                                 semData = semData,
                                 measure = measure),
                          .options = furrr_options(seed = TRUE))

    sim_lst <- map_dbl(sim_lst, ~.x[[1]])

    ind_mat <- map(ind_combi, reduce, cbind)

    ind_mat <- do.call('rbind', ind_mat)

    res <- matrix(1, nrow = length(x), ncol = length(x))

    for(i in seq_along(sim_lst)) {

      res[ind_mat[i, 1], ind_mat[i, 2]] <- sim_lst[i]
      res[ind_mat[i, 2], ind_mat[i, 1]] <- sim_lst[i]

    }

    rownames(res) <- x
    colnames(res) <- x

    if(!is.null(pad_na)) {

      res <- ifelse(is.na(res), pad_na, res)

    }

    if(as_dist) res <- as.dist(1 - res)

    return(res)

  }

# END  -------
