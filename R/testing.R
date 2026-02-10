# Functions for testing gene expression differences between two groups.
# Retired and removed

#' @include imports.R

  NULL

# Testers for two groups ------

#' Test gene expression differences between two groups.
#'
#' @description
#' The testing functions are not active anymore, because we offer a faster and
#' more versatile tools for high-throughput statistical hypothesis testing.
#' Please use functions from `fastTest` package instead
#' (https://github.com/PiotrTymoszuk/fastTest/)
#'
#' @details
#' These functions are not active anymore and return `NULL` with a warning.
#'
#' @param data a data frame with the expression values.
#' @param split_fct name of the group-defining variable, needs to be a factor.
#' @param variables a vector with dependent variable names.
#' @param type test type, currently 't' for T test (default) or 'wilcox' for
#' the Wilcoxon/Mann-Whitney test.
#' @param adj_method p value adjustment method,
#' see: \code{\link[stats]{p.adjust}}.
#' @param .parallel should the testing be run in parallel?
#' @param ... additional arguments passed to \code{\link[stats]{t.test}} or
#' \code{\link[stats]{wilcox.test}}.
#'
#' @return `NULL` with a warning.
#'
#' @export

  test_two_groups <- function(data,
                              split_fct,
                              variables,
                              type = c('t', 'wilcox'),
                              adj_method = 'BH',
                              .parallel = FALSE, ...) {

    warning(paste("This function is not active anymore.",
                  "see: https://github.com/PiotrTymoszuk/fastTest/",
                  "for faster and more robust alternatives."),,
            call. = FALSE)

    return(NULL)

  }

# One-way ANOVA with LM -----

#' Test gene expression differences between multiple groups with one-way ANOVA.
#'
#' @description
#' The testing functions are not active anymore, because we offer a faster and
#' more versatile tools for high-throughput statistical hypothesis testing.
#' Please use functions from `fastTest` package instead
#' (https://github.com/PiotrTymoszuk/fastTest/)
#'
#' @details
#' These functions are not active anymore and return `NULL` with a warning
#'
#' @param data a data frame with the expression values.
#' @param split_fct name of the group-defining variable, needs to be a factor.
#' @param variables a vector with dependent variable names.
#' @param confounder name of the confounding variable. Defaults to NULL, which
#' means that no confounder is included in the model.
#' @param adj_method p value adjustment method,
#' see: \code{\link[stats]{p.adjust}}.
#' @param .parallel should the testing be run in parallel?
#' @param ... additional arguments passed to \code{\link[stats]{aov}}.
#'
#' @return `NULL` with a warning.
#'
#' @export

  test_anova <- function(data,
                         split_fct,
                         variables,
                         confounder = NULL,
                         adj_method = 'BH',
                         .parallel = FALSE, ...) {

    warning(paste("This function is not active anymore.",
                  "see: https://github.com/PiotrTymoszuk/fastTest/",
                  "for faster and more robust alternatives."),,
            call. = FALSE)

    return(NULL)

  }

# Negative binomial modeling ------

#' Test gene expression differences with negative binomial modeling.
#'
#' @description
#' The testing functions are not active anymore, because we offer a faster and
#' more versatile tools for high-throughput statistical hypothesis testing.
#' Please use functions from `fastTest` package instead
#' (https://github.com/PiotrTymoszuk/fastTest/).
#' When looking specifically for negative binomial modeling of gene expression
#' data, please use EdgeR (https://bioconductor.org/packages/release/bioc/html/edgeR.html)
#' or DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
#' instead.
#'
#' @details
#' These functions are not active anymore and return `NULL` with a warning
#'
#' @param data a data frame with the expression values.
#' @param split_fct name of the group-defining variable, needs to be a factor.
#' @param variables a vector with dependent variable names.
#' @param confounder name of the confounding variable. Defaults to NULL, which
#' means that no confounder is include in the model.
#' @param dev_test name of the test used for analysis of deviance. One of:
#' 'Rao', 'LRT', 'Chisq' (default), 'F', 'Cp'.
#' @param adj_method p value adjustment method,
#' see: \code{\link[stats]{p.adjust}}.
#' @param .parallel should the testing be run in parallel?
#' @param ... additional arguments passed to \code{\link[MASS]{glm.nb}}.
#'
#' @return `NULL` with a warning.
#'
#' @export

  test_nb <- function(data,
                      split_fct,
                      variables,
                      confounder = NULL,
                      adj_method = 'BH',
                      dev_test = 'Chisq',
                      .parallel = FALSE, ...) {

    warning(paste("This function is not active anymore.",
                  "see: https://github.com/PiotrTymoszuk/fastTest/, ",
                  "https://bioconductor.org/packages/release/bioc/html/edgeR.html, and ",
                  "https://bioconductor.org/packages/release/bioc/html/DESeq2.html",
                  "for faster and more robust alternatives."),,
            call. = FALSE)

    return(NULL)

  }

# END -----
