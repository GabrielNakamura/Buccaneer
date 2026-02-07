#' Longevities for 133 Canidae species
#'
#' Data frame containing Time of Speciation (TS) and Time of Extinction (TE) for
#'     133 Canidae species
#'
#'
#' @format
#' A data frame object with 133 rows and 3 columns:
#' \describe{
#'     \item{species}{Species names}
#'     \item{TS}{Time of Speciation}
#'     \item{TE}{Time of Extinction}
#' }
#'
#' @source <https://academic.oup.com/evolut/article/79/9/1835/8140865#535027902>
"df_longevities_canidae"

#' Trait data for 133 Canidae species
#'
#' Data frame containing five traits measured for 133 Canidae species
#'
#' @format
#' A data frame object with 133 rows and 5 columns, each row
#'     correspond to a species name their respective traits:
#' \describe{
#'     \item{species}{Species name}
#'     \item{diet}{Numeric indicating the diet category of species}
#'     \item{LD1}{Numeric indicating the level of hypercarnivore of species}
#'     \item{log_mass}{Numeric indicating the mass of species in log scale}
#'     \item{diet_cat}{Character indicating the diet category of species.
#'         hyper indicate hypercarnivores, meso indicates mesocarnivores and
#'         hypo hypocarnivores}
#' }
#'
#' @source <https://academic.oup.com/evolut/article/79/9/1835/8140865#535027902>
"traits_canidae"


#' Occurrence records
#'
#' Data frame containing spatial occurrence records of 133 species of canidae
#'     in North America
#'
#' @format
#' A data frame object with 133 rows and 14 columns:
#' \describe{
#'     \item{species}{Species name}
#'     \item{MaxT}{Maximum (oldest) age estimate for each occurrence record}
#'     \item{MinT}{the minimum (youngest) age estimate for each occurrence record}
#'     \item{lng}{The longitude coordinate of the occurrence record}
#'     \item{lat}{The latitude coordinate of the occurrence record}
#'     \item{Site}{Name of the paleontological site location identifier}
#'     \item{max_low_res}{Same as MaxT but with lower temporal precision}
#'     \item{min_lo_res}{Same as MinT but with lower temporal precision}
#'     \item{site.char}{Name of the paleontological site identifier}
#'     \item{midpoint}{Age midpoint of occurrence record}
#' }
#'
#' @source <https://academic.oup.com/evolut/article/79/9/1835/8140865#535027902>
"df_occurrence_canidae"
