% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{box_from_stats}
\alias{box_from_stats}
\title{Plot a box plot given median, interquartile and 95\% percentile range.}
\usage{
box_from_stats(
  data,
  x_variable,
  fill_variable = NULL,
  fill_color = "steelblue",
  width = 0.8,
  alpha = 1,
  plot_title = NULL,
  plot_subtitle = NULL,
  x_lab = NULL,
  y_lab = NULL,
  cust_theme = microViz::theme_micro()
)
}
\arguments{
\item{data}{a data frame with the following columns: \code{median},
\code{perc25}, \code{perc75}, \code{perc025}, and \code{perc975} storing the stats.}

\item{x_variable}{name of the variable presented in the X axis.}

\item{fill_variable}{name of the variable defining color of the box plots.}

\item{fill_color}{fill color for the boxes, ignored if \code{fill_variable = NULL}.}

\item{width}{width of the boxes.}

\item{alpha}{alpha of the box plots.}

\item{plot_title}{plot title,}

\item{plot_subtitle}{plot subtitle.}

\item{x_lab}{x axis title.}

\item{y_lab}{y axis title.}

\item{cust_theme}{custom \code{ggplot} theme}
}
\value{
a \code{ggplot} object.
}
\description{
Generates a \code{ggplot} box plot from distribution statistics: median,
25\% and 75\% percentiles, and 2.5\% and 95\% percentiles.
The function may be useful at plotting distribution statistic for very large
data sets (e.g. whole transcriptome). In this case, a box plot made with
pre-calculated distribution statistic is going to spare a lot of memory.
}
