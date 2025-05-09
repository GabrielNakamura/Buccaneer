% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_occ_fossil.R
\name{clean_occ_fossil}
\alias{clean_occ_fossil}
\title{Cleaning and flagging fossil occurrence records according different criteria}
\usage{
clean_occ_fossil(
  df.occ.fossil,
  method.ages = c("midpoint", "upper", "lower"),
  thresh.age.range = 10,
  species = "species",
  Max.age = "Maximum_Age",
  Min.age = "Minimum_Age",
  remove.sub.species = TRUE,
  comp.TS.TE = TRUE,
  lat = NULL,
  lng = NULL,
  site = NULL,
  group = NULL,
  trait = NULL
)
}
\arguments{
\item{...}{}
}
\description{
Cleaning and flagging fossil occurrence records according different criteria
}
