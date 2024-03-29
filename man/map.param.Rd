% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_map.R
\name{map.param}
\alias{map.param}
\title{Generate MAP using the information on data}
\usage{
map.param(data, skel, ...)
}
\arguments{
\item{data}{A dataframe}

\item{skel}{An instance of skeleton of MAP}

\item{...}{Others}
}
\value{
An instance of MAP
}
\description{
Generate MAP randomly and adjust parameters to fit its first moment to
the first moment of data.
}
