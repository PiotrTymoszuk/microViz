% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiple_plotting.R
\name{plot_common_change}
\alias{plot_common_change}
\title{Before-after or change plots for multiple data frames.}
\usage{
plot_common_change(
  data,
  variables,
  split_fct,
  normalize = TRUE,
  norm_center = c("mean", "median"),
  average_fun = colMeans,
  show_whiskers = FALSE,
  whisker_size = 0.5,
  probs = c(0.025, 0.975),
  dodge = 0,
  point_size = 2,
  point_alpha = 1,
  linetype = "solid",
  line_alpha = 0.75,
  line_color = "gray70",
  show_txt = TRUE,
  txt_size = 2.75,
  txt_position = c("left", "right"),
  txt_color = NULL,
  txt_face = "plain",
  labeller = identity,
  cust_theme = microViz::theme_micro(),
  plot_titles = variables,
  plot_subtitles = NULL,
  x_lab = split_fct,
  y_lab = "average Z-score",
  palette = NULL
)
}
\arguments{
\item{data}{a list of data frames.}

\item{variables}{a vector with names of the variables of interest.}

\item{split_fct}{name of the splitting factor.}

\item{normalize}{logical, should the data frame variables be normalized prior
to plotting?}

\item{norm_center}{defines centering of the variable during scaling: mean
(default) or median. Ignored if \code{normalize = FALSE}.}

\item{average_fun}{a function used to calculate average values of the
variables within levels of the splitting factor.}

\item{show_whiskers}{logical, should whiskers representing a quantile range
be displayed in the plot? Defaults to false.}

\item{whisker_size}{size of the whiskers.}

\item{probs}{a numeric vector of two quantiles defining the range for the
whiskers. Ignored if \code{show_whiskers = FALSE}.}

\item{dodge}{a numeric specifying dodging of the points, whiskers, labels,
and lines.}

\item{point_size}{size of the data points.}

\item{point_alpha}{alpha of the points and, optionally, whiskers.}

\item{linetype}{type of the line connecting the points.}

\item{line_alpha}{alpha of the connecting lines.}

\item{line_color}{color of the connecting lines.}

\item{show_txt}{logical, should labels with names of the data sets be
displayed?}

\item{txt_size}{size of the label text.}

\item{txt_position}{placement of the text labels.}

\item{txt_color}{color of the text. If \code{NULL}, the color is specified by
levels of the splitting factor.}

\item{txt_face}{fornt face of the text labels.}

\item{labeller}{a function transforming the data set names to the text labels
in the plot. Defaults to \code{identity}, i.e. the names are displayed as the are.}

\item{cust_theme}{a custom \code{ggplot} theme.}

\item{plot_titles}{a character vector with plot titles.}

\item{plot_subtitles}{a character vector with lot subtitles.}

\item{x_lab}{X axis title.}

\item{y_lab}{Y axis title.}

\item{palette}{colors for levels of the splitting factor. If \code{NULL}, default
R color and fill scales will be used.}
}
\value{
a list of \code{ggplot} objects.
}
\description{
Draws a list of plots of average statistics (i.e. mean or median) of numeric
variables split by a factor splitting variable. Averages for the same data
set are connected by lines.
}
