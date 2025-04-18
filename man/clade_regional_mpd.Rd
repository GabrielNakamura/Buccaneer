% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clade_regional_mpd.R
\name{clade_regional_mpd}
\alias{clade_regional_mpd}
\title{Compute regional clade Mean Pairwise Distance in a time series}
\usage{
clade_regional_mpd(
  df.TS.TE,
  time.slice,
  trait,
  round.digits = 1,
  species = "species",
  TS = "TS",
  TE = "TE",
  compute.ses.mpd = TRUE,
  null.model = "taxa.labels",
  runs = 999,
  group = NULL,
  group.focal.compare = NULL,
  type.comparison = NULL
)
}
\arguments{
\item{type.comparison}{}
}
\description{
Compute regional clade Mean Pairwise Distance in a time series
}
