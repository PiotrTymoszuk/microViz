% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{brca}
\alias{brca}
\title{Breast cancer transcriptome data.}
\format{
A data frame with 157 rows and 38957 variables:
\itemize{
\item \strong{sample_id}: sample identifier
\item \strong{patient_id}: patient identifier
\item \strong{timepoint}: time point of the sample collection
\item \strong{metastasis}: presence of distant metastases
\item \strong{histology}: tumor histology
\item \strong{er_status}: estrogen receptor (ER) positivity of the cancer
\item gene expression values as RSEM in columns named after the official HGNC
gene symbol
}
}
\source{
cBioportal, the Metastatic Breast Cancer Project,
https://www.cbioportal.org/study/summary?id=brca_mbcproject_2022
}
\usage{
data(brca)
}
\description{
Whole-genome expression data of breast cancers collected by the
Metastatic Breast Cancer Project.
}
