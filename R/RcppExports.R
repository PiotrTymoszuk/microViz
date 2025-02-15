# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

colMed <- function(x, na_rm = TRUE) {
    .Call('_microViz_colMed', PACKAGE = 'microViz', x, na_rm)
}

colGeoMean <- function(x, na_rm = TRUE) {
    .Call('_microViz_colGeoMean', PACKAGE = 'microViz', x, na_rm)
}

colHarmMean <- function(x, na_rm = TRUE) {
    .Call('_microViz_colHarmMean', PACKAGE = 'microViz', x, na_rm)
}

colVariance <- function(x, na_rm = TRUE) {
    .Call('_microViz_colVariance', PACKAGE = 'microViz', x, na_rm)
}

colSD <- function(x, na_rm = TRUE) {
    .Call('_microViz_colSD', PACKAGE = 'microViz', x, na_rm)
}

colGi <- function(x, unbiased = TRUE, na_rm = TRUE) {
    .Call('_microViz_colGi', PACKAGE = 'microViz', x, unbiased, na_rm)
}

colQuant <- function(x, probs, na_rm = TRUE) {
    .Call('_microViz_colQuant', PACKAGE = 'microViz', x, probs, na_rm)
}

colPerCi <- function(x, conf_level, na_rm = TRUE) {
    .Call('_microViz_colPerCi', PACKAGE = 'microViz', x, conf_level, na_rm)
}

colBcaCi <- function(x, conf_level, na_rm = TRUE) {
    .Call('_microViz_colBcaCi', PACKAGE = 'microViz', x, conf_level, na_rm)
}

colFreqRatio <- function(x) {
    .Call('_microViz_colFreqRatio', PACKAGE = 'microViz', x)
}

colPercUnique <- function(x) {
    .Call('_microViz_colPercUnique', PACKAGE = 'microViz', x)
}

Delta <- function(x, mu) {
    .Call('_microViz_Delta', PACKAGE = 'microViz', x, mu)
}

colMi <- function(x, na_rm = TRUE) {
    .Call('_microViz_colMi', PACKAGE = 'microViz', x, na_rm)
}

colMa <- function(x, na_rm = TRUE) {
    .Call('_microViz_colMa', PACKAGE = 'microViz', x, na_rm)
}

colMis <- function(x) {
    .Call('_microViz_colMis', PACKAGE = 'microViz', x)
}

minMaxVec <- function(x) {
    .Call('_microViz_minMaxVec', PACKAGE = 'microViz', x)
}

minMaxCols <- function(x) {
    .Call('_microViz_minMaxCols', PACKAGE = 'microViz', x)
}

minMaxRows <- function(x) {
    .Call('_microViz_minMaxRows', PACKAGE = 'microViz', x)
}

zScoreVec <- function(x, center = "mean", dispersion = "sd") {
    .Call('_microViz_zScoreVec', PACKAGE = 'microViz', x, center, dispersion)
}

zScoreCols <- function(x, center = "mean", dispersion = "sd") {
    .Call('_microViz_zScoreCols', PACKAGE = 'microViz', x, center, dispersion)
}

zScoreRows <- function(x, center = "mean", dispersion = "sd") {
    .Call('_microViz_zScoreRows', PACKAGE = 'microViz', x, center, dispersion)
}

offsetVector <- function(x) {
    .Call('_microViz_offsetVector', PACKAGE = 'microViz', x)
}

offsetMatrix <- function(x) {
    .Call('_microViz_offsetMatrix', PACKAGE = 'microViz', x)
}

poolVarPaired <- function(x, y) {
    .Call('_microViz_poolVarPaired', PACKAGE = 'microViz', x, y)
}

poolVarWelch <- function(x, y) {
    .Call('_microViz_poolVarWelch', PACKAGE = 'microViz', x, y)
}

poolVarStandard <- function(x, y) {
    .Call('_microViz_poolVarStandard', PACKAGE = 'microViz', x, y)
}

poolVarLM <- function(x) {
    .Call('_microViz_poolVarLM', PACKAGE = 'microViz', x)
}

rocThreshold <- function(outcome, marker, threshold, direction = "<", skip_control = FALSE) {
    .Call('_microViz_rocThreshold', PACKAGE = 'microViz', outcome, marker, threshold, direction, skip_control)
}

rocMarker <- function(outcome, marker, direction = "auto", skip_control = FALSE) {
    .Call('_microViz_rocMarker', PACKAGE = 'microViz', outcome, marker, direction, skip_control)
}

aucVec <- function(outcome, marker, direction = "auto", skip_control = FALSE) {
    .Call('_microViz_aucVec', PACKAGE = 'microViz', outcome, marker, direction, skip_control)
}

aucMtx <- function(outcome, markers, direction = "auto") {
    .Call('_microViz_aucMtx', PACKAGE = 'microViz', outcome, markers, direction)
}

rowMed <- function(x, na_rm = TRUE) {
    .Call('_microViz_rowMed', PACKAGE = 'microViz', x, na_rm)
}

rowGeoMean <- function(x, na_rm = TRUE) {
    .Call('_microViz_rowGeoMean', PACKAGE = 'microViz', x, na_rm)
}

rowHarmMean <- function(x, na_rm = TRUE) {
    .Call('_microViz_rowHarmMean', PACKAGE = 'microViz', x, na_rm)
}

rowVariance <- function(x, na_rm = TRUE) {
    .Call('_microViz_rowVariance', PACKAGE = 'microViz', x, na_rm)
}

rowSD <- function(x, na_rm = TRUE) {
    .Call('_microViz_rowSD', PACKAGE = 'microViz', x, na_rm)
}

rowGi <- function(x, unbiased = TRUE, na_rm = TRUE) {
    .Call('_microViz_rowGi', PACKAGE = 'microViz', x, unbiased, na_rm)
}

rowQuant <- function(x, probs, na_rm = TRUE) {
    .Call('_microViz_rowQuant', PACKAGE = 'microViz', x, probs, na_rm)
}

rowPerCi <- function(x, conf_level, na_rm = TRUE) {
    .Call('_microViz_rowPerCi', PACKAGE = 'microViz', x, conf_level, na_rm)
}

rowBcaCi <- function(x, conf_level, na_rm = TRUE) {
    .Call('_microViz_rowBcaCi', PACKAGE = 'microViz', x, conf_level, na_rm)
}

rowFreqRatio <- function(x) {
    .Call('_microViz_rowFreqRatio', PACKAGE = 'microViz', x)
}

rowPercUnique <- function(x) {
    .Call('_microViz_rowPercUnique', PACKAGE = 'microViz', x)
}

rowMi <- function(x, na_rm = TRUE) {
    .Call('_microViz_rowMi', PACKAGE = 'microViz', x, na_rm)
}

rowMa <- function(x, na_rm = TRUE) {
    .Call('_microViz_rowMa', PACKAGE = 'microViz', x, na_rm)
}

rowMis <- function(x) {
    .Call('_microViz_rowMis', PACKAGE = 'microViz', x)
}

vecSimilarity <- function(x, y, method = "jaccard") {
    .Call('_microViz_vecSimilarity', PACKAGE = 'microViz', x, y, method)
}

vecTversky <- function(x, y, a = 1, b = 1) {
    .Call('_microViz_vecTversky', PACKAGE = 'microViz', x, y, a, b)
}

setSimilarity <- function(x, method = "jaccard") {
    .Call('_microViz_setSimilarity', PACKAGE = 'microViz', x, method)
}

setTversky <- function(x, a, b) {
    .Call('_microViz_setTversky', PACKAGE = 'microViz', x, a, b)
}

Median <- function(x) {
    .Call('_microViz_Median', PACKAGE = 'microViz', x)
}

geoMean <- function(x) {
    .Call('_microViz_geoMean', PACKAGE = 'microViz', x)
}

harmMean <- function(x) {
    .Call('_microViz_harmMean', PACKAGE = 'microViz', x)
}

Var <- function(x) {
    .Call('_microViz_Var', PACKAGE = 'microViz', x)
}

SD <- function(x) {
    .Call('_microViz_SD', PACKAGE = 'microViz', x)
}

SEM <- function(x) {
    .Call('_microViz_SEM', PACKAGE = 'microViz', x)
}

GiniCpp <- function(x, unbiased = TRUE) {
    .Call('_microViz_GiniCpp', PACKAGE = 'microViz', x, unbiased)
}

Quantile <- function(x, probs) {
    .Call('_microViz_Quantile', PACKAGE = 'microViz', x, probs)
}

perci <- function(theta, conf_level = 0.95) {
    .Call('_microViz_perci', PACKAGE = 'microViz', theta, conf_level)
}

bca <- function(theta, conf_level = 0.95) {
    .Call('_microViz_bca', PACKAGE = 'microViz', theta, conf_level)
}

Table <- function(x) {
    .Call('_microViz_Table', PACKAGE = 'microViz', x)
}

freqRatioCpp <- function(x) {
    .Call('_microViz_freqRatioCpp', PACKAGE = 'microViz', x)
}

percUniqueCpp <- function(x) {
    .Call('_microViz_percUniqueCpp', PACKAGE = 'microViz', x)
}

fillMat <- function(x, ind, dim) {
    .Call('_microViz_fillMat', PACKAGE = 'microViz', x, ind, dim)
}

missingCpp <- function(x) {
    .Call('_microViz_missingCpp', PACKAGE = 'microViz', x)
}

sortByFirstColumn <- function(mat, decreasing = FALSE) {
    .Call('_microViz_sortByFirstColumn', PACKAGE = 'microViz', mat, decreasing)
}

sortIndexes <- function(x, decreasing = FALSE) {
    .Call('_microViz_sortIndexes', PACKAGE = 'microViz', x, decreasing)
}

trapezoidal_integration <- function(x, y) {
    .Call('_microViz_trapezoidal_integration', PACKAGE = 'microViz', x, y)
}

