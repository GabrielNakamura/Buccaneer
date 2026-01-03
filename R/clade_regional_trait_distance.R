#' Calculate Regional Mean Trait Distances Across Time Slices
#'
#' This function computes mean pairwise trait distances between species at the
#' regional scale across different time slices. For each time slice, it calculates
#' the mean distance among all species whose temporal ranges overlap with that
#' interval, supporting both mean nearest neighbor distance (MNND) and mean
#' pairwise distance (MPD) metrics. The function can perform comparisons between
#' groups (e.g., clades, families) or within a single group.
#'
#' @param df.TS.TE A data frame containing species temporal and trait data with
#'     at least three columns: species names, origination times (TS), extinction
#'     times (TE). Trait values are optional in this object.
#'     Additional columns may include group assignments.
#' @param time.slice Numeric. The time interval (in the same units as TS and TE)
#'     between consecutive time slices for temporal binning.
#' @param dist.trait A distance matrix object (class \code{dist} or \code{matrix})
#'     containing pairwise trait distances between species. Row and column names
#'     must match species names in \code{df.TS.TE}. If NULL, distances will be
#'     computed from the trait column using Euclidean distance.
#' @param nearest.taxon Numeric or character. The number of nearest neighbors to
#'     consider when calculating mean distances. Use \code{1} for mean nearest
#'     neighbor distance (MNND), or \code{"all"} for mean pairwise distance (MPD).
#' @param trait Character. The name of the column in \code{df.TS.TE} containing
#'     trait values. If NULL (default), \code{dist.trait} must be provided.
#' @param round.digits Integer. The number of decimal places to round time slice
#'     values. Default is 1. This affects temporal binning precision.
#' @param species Character. The name of the column in \code{df.TS.TE} containing
#'     species identifiers. Default is "species".
#' @param TS Character. The name of the column in \code{df.TS.TE} containing
#'     origination (first appearance) times for each species. Default is "TS".
#' @param TE Character. The name of the column in \code{df.TS.TE} containing
#'     extinction (last appearance) times for each species. Default is "TE".
#' @param group Character. The name of the column in \code{df.TS.TE} containing
#'     group assignments for species (e.g., clade, family). Required if using
#'     \code{group.focal.compare}. Default is NULL.
#' @param group.focal.compare Character vector of length 2. The first element
#'     specifies the focal group and the second specifies the comparison group.
#'     If NULL (default), distances are calculated across all species regardless
#'     of group membership.
#' @param type.comparison Character. Specifies the type of distance comparison:
#'     \itemize{
#'       \item \code{"between"}: Calculate distances only between species from
#'             the focal and comparison groups.
#'       \item \code{"within"}: Calculate distances only among species within
#'             the focal group.
#'       \item NULL (default): Calculate distances among all species together.
#'     }
#'
#' @return A data frame with three columns:
#'   \item{mean.distance}{Numeric. The mean trait distance among species present
#'       in each time slice. Returns NA when insufficient data is available
#'       (e.g., only one species present, or no representatives from required groups).}
#'   \item{var.distance}{Numeric. The variance of trait distances among species
#'       in each time slice. Returns NA when insufficient data is available.}
#'   \item{time.slice}{Numeric. The time point representing each slice, typically
#'       the upper (older) boundary of the time bin.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates a sequence of time slices from maximum TS to minimum TE
#'   \item Generates regional co-occurrence matrices using \code{aux_matrix_regional_coex()}
#'   \item For each time slice, identifies species whose temporal ranges overlap
#'   \item Computes pairwise trait distances among overlapping species
#'   \item Optionally filters comparisons by group membership
#'   \item Calculates mean and variance of distances based on \code{nearest.taxon}
#' }
#'
#' Distance calculation options:
#' \itemize{
#'   \item \strong{MPD (nearest.taxon = "all")}: Calculates mean pairwise distance
#'         considering all pairwise comparisons among species
#'   \item \strong{MNND (nearest.taxon = 1)}: Calculates mean nearest neighbor
#'         distance using only the closest species for each focal species
#'   \item \strong{Threshold (nearest.taxon = n)}: Uses the \code{n} nearest
#'         neighbors for each focal species
#' }
#'
#' Missing values (NA) are returned for time slices where:
#' \itemize{
#'   \item Only one species is present
#'   \item No species from required groups are present (when using group comparisons)
#'   \item Insufficient data for distance calculations
#' }
#'
#' This function calculates distances at the regional (pool) scale, considering
#' all species present during each time slice regardless of their geographic
#' distribution. For site-level distance calculations, see \code{IndivSpec_site_distance()}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example fossil data with traits
#' df_temporal <- data.frame(
#'   species = c("sp1", "sp2", "sp3", "sp4"),
#'   TS = c(100, 95, 90, 85),
#'   TE = c(50, 45, 40, 35),
#'   trait = c(1.2, 2.5, 3.1, 4.0),
#'   group = c("A", "A", "B", "B")
#' )
#'
#' # Calculate regional MPD through time
#' result_mpd <- clade_regional_distance(
#'   df.TS.TE = df_temporal,
#'   time.slice = 10,
#'   dist.trait = NULL,
#'   nearest.taxon = "all",
#'   trait = "trait"
#' )
#'
#' # View results
#' head(result_mpd)
#'
#' # Plot mean distance through time
#' plot(result_mpd$time.slice,
#'      result_mpd$mean.distance,
#'      type = "l",
#'      xlab = "Time (Ma)",
#'      ylab = "Mean Trait Distance",
#'      main = "Regional Trait Distance Through Time")
#'
#' # Calculate MNND between groups
#' result_between <- clade_regional_distance(
#'   df.TS.TE = df_temporal,
#'   time.slice = 10,
#'   dist.trait = NULL,
#'   nearest.taxon = 1,
#'   trait = "trait",
#'   group = "group",
#'   group.focal.compare = c("A", "B"),
#'   type.comparison = "between"
#' )
#'
#' # Calculate distances within a single group
#' result_within <- clade_regional_distance(
#'   df.TS.TE = df_temporal,
#'   time.slice = 10,
#'   dist.trait = NULL,
#'   nearest.taxon = "all",
#'   trait = "trait",
#'   group = "group",
#'   group.focal.compare = c("A", "B"),
#'   type.comparison = "within"
#' )
#'
#' # Using a pre-computed distance matrix
#' dist_matrix <- dist(df_temporal$trait)
#' result_custom_dist <- clade_regional_distance(
#'   df.TS.TE = df_temporal,
#'   time.slice = 10,
#'   dist.trait = dist_matrix,
#'   nearest.taxon = "all"
#' )
#' }
clade_regional_distance <-
  function(df.TS.TE,
           time.slice,
           dist.trait,
           nearest.taxon,
           trait = NULL,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE",
           group = NULL,
           group.focal.compare = NULL,
           type.comparison = NULL
  ){

    # subseting columns
    if(!is.null(group) == TRUE){
      df.TS.TE <- df.TS.TE[, c(species, trait, TS, TE, group)]
      if(is.null(trait) == TRUE){
        colnames(df.TS.TE) <- c("species", "TS", "TE", "group")
      } else{
        colnames(df.TS.TE) <- c("species", "trait", "TS", "TE", "group")
      }
    } else{
      df.TS.TE <- df.TS.TE[, c(species, trait, TS, TE)]
      if(is.null(trait) == TRUE){
        colnames(df.TS.TE) <- c("species", "TS", "TE")
      } else{
        colnames(df.TS.TE) <- c("species", "trait", "TS", "TE")
      }
    }

    # creating time slices
    seq_interval <- seq(from = ceiling(max(df.TS.TE[, "TS"])),
                        to = ceiling(min(df.TS.TE[, "TE"])),
                        by = -time.slice)

    # co-occurrence matrix
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE,
                               time.slice,
                               round.digits = 1,
                               species = "species",
                               TS = "TS",
                               TE = "TE")

    # community composition matrix

    spp_slice <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })


    # trait distance
    if(!is.null(dist.trait) == TRUE){
      matrix_dist_trait <- as.matrix(dist.trait)
    } else{
      matrix_dist_trait <-
        as.matrix(dist(x = df.TS.TE[, "trait"], method = "euclidean", upper = T, diag = T))
    }
    rownames(matrix_dist_trait) <- df.TS.TE$species
    colnames(matrix_dist_trait) <- df.TS.TE$species


    # modified trait matrix containing group comparison
    if(!is.null(group.focal.compare) == TRUE){
      focal <- group.focal.compare[1]
      compare <- group.focal.compare[2]
      spp_focal <- df.TS.TE[which(df.TS.TE$group == focal), "species"]
      spp_compare <- df.TS.TE[which(df.TS.TE$group == compare), "species"]

      if(type.comparison == "between"){# comparison between species of two groups
        matrix_dist_trait_comp <- matrix_dist_trait[spp_focal, spp_compare] # focal speices in lines and comparison in columns
      }
      if(type.comparison == "within"){ # comparison only within the focal group
        matrix_dist_trait_comp <- matrix_dist_trait[spp_focal, spp_focal]
      }
    } else{
      matrix_dist_trait_comp <- matrix_dist_trait
    }

    # filtering by timeslices

    # spp_slice <- spp_slice[which(seq_interval == 10.00)] # test for slice 10
    mean_dist_timeslice <- vector(length = length(spp_slice))
    var_dist_timeslice <- vector(length = length(spp_slice))
    for(i in 1:length(spp_slice)){

      if(length(spp_slice[[i]]) == 1){
        mean_dist_timeslice[i] <- NA
        var_dist_timeslice[i] <- NA
      } else{
        rows <- match(spp_slice[i][[1]], rownames(matrix_dist_trait_comp))
        cols <- match(spp_slice[i][[1]], colnames(matrix_dist_trait_comp))
        # checking the presence of representants of focal and comparison groups
        if(all(is.na(rows)) | all(is.na(cols))){
          mean_dist_timeslice[i] <- NA
          var_dist_timeslice[i] <- NA
        } else{
          # if there is only one representant of focal and comparison group
          if(length(rows) == 1 | length(cols) == 1){
            mean_dist_timeslice[i] <- NA
            var_dist_timeslice[i] <- NA
          } else{
            matrix_dist_comp2 <- matrix_dist_trait_comp[na.omit(rows),
                                                        na.omit(cols)]
            if(!is.matrix(matrix_dist_comp2) == TRUE){
              matrix_dist_comp2 <- as.matrix(matrix_dist_comp2)

            }
            matrix_dist_comp3 <- apply(matrix_dist_comp2, 1, function(x) sort(x))
            if(is.null(type.comparison) == TRUE){
              matrix_dist_comp3 <- matrix_dist_comp3[-1, ]
            } else{
              if(type.comparison == "within"){
                if(is.vector(matrix_dist_comp3) == TRUE){
                  matrix_dist_comp3 <- matrix_dist_comp3
                } else{
                  matrix_dist_comp3 <- matrix_dist_comp3[-1,]
                }
              }
            }

            # filtering by the threshold and keeping only the n nearest species
            # if the matrix has only one closest distance
            if(is.vector(matrix_dist_comp3) == TRUE){
              mean_dist_timeslice[i] <- mean(matrix_dist_comp3)
              var_dist_timeslice[i] <- var(as.vector(matrix_dist_comp3))
            } else{ # if the matrix has less close taxon than the threshold get the total number of comparison of the matrix
              if(nearest.taxon == "all"){ # used to compute mpd - considering all distances
                mean_dist_timeslice[i] <- mean(as.vector(matrix_dist_comp3[1:nrow(matrix_dist_comp3), ]))
                var_dist_timeslice[i] <- var(as.vector(matrix_dist_comp3[1:nrow(matrix_dist_comp3), ]))
              } else{
                if(is.numeric(nearest.taxon) == TRUE){ # used to compute mean distances considering thresholds
                  if(nearest.taxon <= dim(matrix_dist_comp3)[1]){
                    mean_dist_timeslice[i] <- mean(as.vector(matrix_dist_comp3[1:nearest.taxon, ]))
                    var_dist_timeslice[i] <- var(as.vector(matrix_dist_comp3[1:nearest.taxon, ]))
                  } else{
                    tmp_nearest_taxon <- dim(matrix_dist_comp3)[1] # n neighbours
                    mean_dist_timeslice[i] <- mean(as.vector(matrix_dist_comp3[1:tmp_nearest_taxon, ]))
                    var_dist_timeslice[i] <- var(as.vector(matrix_dist_comp3[1:tmp_nearest_taxon, ]))
                  } # if the nearest neighbor value are higher than the dimension of the matrix
                }
              }
            }
          }
        }
      }
    }

    # data frame with the results
    df_res <-
      data.frame(mean.distance = mean_dist_timeslice,
                 var.distance = var_dist_timeslice,
                 time.slice = seq_interval)


    # data frame with results
    return(df_res)

  }
