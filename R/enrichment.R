# GO enrichment analysis

#' @include imports.R

  NULL

# Numeric analyses -------

#' Gene Ontology enrichment analysis.
#'
#' @description
#' The function performs Gene Ontology (GO) enrichment analysis with the
#' `limma` function \code{\link[limma]{goana}} and enhances its results by
#' computing odds ratios (OR) for enrichment of the differentially regulated
#' genes associated with the particular GO term over the genome or universe
#' frequency.
#'
#' @details
#' The OR is computed with the following formula:
#'
#' \deqn{OR = \frac{N_{GO DGE} \times N_{universe}}{N_{GO} \times N_{DGE}}}
#'
#' where \eqn{N_{GO DGE}} stands for the number of differentially regulated
#' genes associated with the particular GO term, \eqn{N_{universe}} is the
#' number of all investigated genes, \eqn{N_{GO}} is the number of all genes
#' associated with the particular GO, and \eqn{N_{DGE}} represents the total
#' number of differentially regulated genes.
#' As such, the OR metric may be used as a measure of effect size of the
#' enrichment. See: \code{\link[effectsize]{interpret_oddsratio}} for details.
#'
#' @param de a vector with Entrez IDs of differentially regulated genes.
#' @param universe an obligatory vector of Entrez IDs for all investigated genes.
#' @param ontology specifies ontology of the GO term to be returned. If `all`
#' (default), all possible GO terms are returned.
#' @param adj_method method of adjustment for multiple testing, `BH`
#' (Benjamini-Hochberg or FDR) by default.
#' @param ... extra arguments passed to \code{\link[limma]{goana}}.
#'
#' @return
#' A data frame with the following columns:
#'
#' * `go_id`: GO term identifier
#' * `term`: name of the GO term
#' * `n_go_dge`: number of the differentially regulated genes associated with
#' the particular GO term
#' * `n_go_total`: total number of genes associated with the particular GO term
#' * `n_dge`: number of the differentially regulated genes.
#' * `n_total`: number of all analyzed genes
#' * `or`: OR for the GO term enrichment over the frequency in the universe
#' * `p_value`: raw, uncorrected p value of the enrichment, as returned
#' by \code{\link[limma]{goana}}.
#' * `p_adjusted`: p value corrected for multiple testing
#'
#' @export

  GOana <- function(de,
                    universe,
                    ontology = c('all', 'BP', 'CC', 'MF'),
                    adj_method = 'BH', ...) {

    ## computes GO term enrichment for the differentially regulated gene set
    ## over the genome/universe with goana()

    ## input control and enrichment analysis

    ontology <- match.arg(ontology[1],
                          c('all', 'BP', 'CC', 'MF'))

    if(is.null(universe)) {

      stop('Please provide a vector of Entrez IDs for the universe.',
           call. = FALSE)

    }

    Term <- NULL
    Ont <- NULL
    DE <- NULL
    N <- NULL
    P.DE <- NULL

    res <- goana(de = de, universe = universe, ...)

    ## formatting of the results ------

    if(ontology != 'all') {

      res <- filter(res, .data[['Ont']] == ontology)

    }

    res <- rownames_to_column(res, 'go_id')

    go_id <- NULL
    term <- NULL
    n_go_dge <- NULL
    n_go_total <- NULL
    n_dge <- NULL
    n_total <- NULL
    or <- NULL
    p_value <- NULL
    p_adjusted <- NULL

    res <-
      transmute(res,
                go_id = go_id,
                term = Term,
                ontology = Ont,
                n_go_dge = DE,
                n_go_total = N,
                n_dge = length(de),
                n_total = length(universe),
                or = (n_go_dge/n_go_total)/(n_dge/n_total),
                p_value = P.DE,
                p_adjusted = p.adjust(P.DE, method = adj_method))

    res <- arrange(res, -or)

    as_tibble(res)

  }

# END -----
