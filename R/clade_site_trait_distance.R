#' Calculate Mean Site-Level Trait Distances Across Time Slices
#'
#' This function computes mean trait distances between co-occurring species at
#' individual sites across different time slices. For each time slice, it calculates
#' the mean distance for all species across all sites where they co-occur, then
#' aggregates these individual species distances to obtain a time slice-level mean
#' and variance.
#'
#' @param df.TS.TE A data frame containing species temporal and trait data with
#'     at least four columns: species names, origination times (TS), extinction
#'     times (TE), and trait values. Additional columns may include group assignments.
#' @param df.occ A data frame containing fossil occurrence records with at least
#'     four columns: species names, minimum age, maximum age, and site location ID.
#'     Each row represents a single occurrence record at a specific site.
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
#' @param species Character. The name of the column in \code{df.TS.TE} and
#'     \code{df.occ} containing species identifiers. Default is "species".
#' @param TS Character. The name of the column in \code{df.TS.TE} containing
#'     origination (first appearance) times for each species. Default is "TS".
#' @param TE Character. The name of the column in \code{df.TS.TE} containing
#'     extinction (last appearance) times for each species. Default is "TE".
#' @param Max.age Character. The name of the column in \code{df.occ} containing
#'     the maximum (oldest) age estimate for each occurrence record. Default is "Max.age".
#' @param Min.age Character. The name of the column in \code{df.occ} containing
#'     the minimum (youngest) age estimate for each occurrence record. Default is "Min.age".
#' @param site Character. The name of the column in \code{df.occ} containing
#'     site location identifiers. Default is "site".
#' @param remove.singletons Logical. Should singleton species (species occurring
#'     alone at a site with no co-occurring species) be excluded from mean and
#'     variance calculations? Default is TRUE. When TRUE, singletons are assigned
#'     a distance of 0 but may be excluded from aggregation depending on implementation.
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
#'   \item{mean.distance}{Numeric. The mean trait distance across all species-site
#'       combinations within each time slice. This is calculated by first computing
#'       mean distances for each species to its co-occurring species at sites,
#'       then averaging these values across all species in the time slice.
#'       Returns NA when there is no species occurrence in the time slice and
#'       NA_singleton when there is only one species occurring }
#'   \item{var.distance}{Numeric. The variance of trait distances across all
#'       species-site combinations within each time slice.}
#'   \item{time.slice}{Numeric. The time point representing each slice, typically
#'       the upper (older) boundary of the time bin.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates time slices from maximum TS to minimum TE
#'   \item Generates regional co-occurrence matrices using \code{aux_matrix_regional_coex()}
#'   \item Determines which species occur at which sites in each time slice
#'   \item Creates site-based co-occurrence matrices using \code{comp_site_cooccurr()}
#'   \item For each species at each site, calculates mean distances to co-occurring species
#'   \item Aggregates individual species distances to obtain time slice-level mean and variance
#'   \item Optionally filters by group membership
#' }
#'
#' Distance calculation hierarchy:
#' \enumerate{
#'   \item \strong{Species-site level}: For each species at each site, calculate
#'         mean distance to co-occurring species (based on \code{nearest.taxon})
#'   \item \strong{Time slice level}: Average all species-site distances within
#'         each time slice to obtain overall mean and variance
#' }
#'
#' Special cases:
#' \itemize{
#'   \item \strong{Singleton species}: Species with no co-occurring taxa at a site
#'         are assigned distance = 0 and flagged (treatment depends on \code{remove.singletons})
#'   \item \strong{Missing data}: Time slices with no occurrences return NA
#'   \item \strong{Group comparisons}: When using \code{group.focal.compare},
#'         only distances between/within specified groups are computed
#' }
#'
#' This function differs from \code{IndivSpec_site_distance()} by aggregating
#' to the time slice level rather than returning individual species-level results.
#' For regional-scale (non-site-based) distances, see \code{clade_regional_distance()}.
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
#' df_occurrences <- data.frame(
#'   species = c("sp1", "sp1", "sp2", "sp3", "sp4", "sp4"),
#'   Max.age = c(100, 95, 95, 90, 85, 85),
#'   Min.age = c(90, 85, 85, 80, 75, 75),
#'   site = c("site1", "site2", "site1", "site1", "site2", "site3")
#' )
#'
#' # Calculate mean site-level MPD through time
#' result_mpd <- clade_site_distance(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
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
#'      ylab = "Mean Site-Level Trait Distance",
#'      main = "Mean Trait Distance at Sites Through Time")
#'
#' # Add variance as error bands
#' polygon(c(result_mpd$time.slice, rev(result_mpd$time.slice)),
#'         c(result_mpd$mean.distance - sqrt(result_mpd$var.distance),
#'           rev(result_mpd$mean.distance + sqrt(result_mpd$var.distance))),
#'         col = rgb(0, 0, 1, 0.2), border = NA)
#'
#' # Calculate MNND at sites
#' result_mnnd <- clade_site_distance(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10,
#'   dist.trait = NULL,
#'   nearest.taxon = 1,
#'   trait = "trait"
#' )
#'
#' # Calculate distances between groups at sites
#' result_between <- clade_site_distance(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10,
#'   dist.trait = NULL,
#'   nearest.taxon = "all",
#'   trait = "trait",
#'   group = "group",
#'   group.focal.compare = c("A", "B"),
#'   type.comparison = "between"
#' )
#'
#' # Calculate distances within a single group at sites
#' result_within <- clade_site_distance(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10,
#'   dist.trait = NULL,
#'   nearest.taxon = "all",
#'   trait = "trait",
#'   group = "group",
#'   group.focal.compare = c("A", "B"),
#'   type.comparison = "within",
#'   remove.singletons = TRUE
#' )
#' }
clade_site_distance <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           dist.trait,
           nearest.taxon,
           group = NULL,
           group.focal.compare = NULL,
           type.comparison = NULL,
           trait = NULL,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE",
           Max.age = "Max.age",
           Min.age = "Min.age",
           site = "site"){
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

    df_occ <-
      df.occ[, c(species, Max.age, Min.age, site)]
    vars <- list(species, Max.age, Min.age, site)
    name_vars <- c("species", "Max.age", "Min.age", "site")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names
    df_occ$site <- as.factor(df_occ$site)

    # Generating time intervals used to compute temporal coexistence
    seq_interval <- seq(from = ceiling(max(df.TS.TE[, "TS"])),
                        to = ceiling(min(df.TS.TE[, "TE"])),
                        by = -time.slice)
    # coexistence matrix
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = 1,
                               species = "species",
                               TS = "TS",
                               TE = "TE")

    # species composition at each timeslice
    spp_slice <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    names(spp_slice) <- format(seq_interval, trim = TRUE, scientific = FALSE)

    # calculating trait distances for all clades

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


    # calculating matrix of species cooccurrence for site

    list_matrix_cooccur_site <-
      comp_site_cooccurr(spp_slice = spp_slice, df.occ = df_occ)


    # Ensure same species and same order in both cooccurrence matrix and distance matrix
    list_dist_spp <-
      lapply(list_matrix_cooccur_site, function(x){
        species_row <- intersect(rownames(x), rownames(matrix_dist_trait_comp))
        species_col <- intersect(colnames(x), colnames(matrix_dist_trait_comp))
        cooccur_matrix <- x[species_row, species_col, drop = FALSE]
        dist_matrix <- matrix_dist_trait_comp[species_row, species_col, drop = FALSE]

        # Create matrix to store the results of mean pairwise distances
        mean_distances <- matrix(NA, nrow = length(species_row), ncol = 1,
                                 dimnames = list(species_row, paste("mean.dist.to.cooccur", nearest.taxon, sep = ".")))

        # calculating distances for all species
        for (sp in species_row) {

          #checking if there are no species in the slice or
          if(length(species_row) == 0 | length(species_col) == 0){
            mean_distances[sp, 1] <- NA
            if(length(species_row) != 0){
              mean_distances[sp, 1] <- "NA_singleton"
            }
          } else{
            if(length(species_row) == 1 & length(species_col) == 1){
              mean_distances[sp, 1] <- dist_matrix
            } else{
              cooccur_species <- names(which(cooccur_matrix[sp, ] > 0 & names(cooccur_matrix[sp, ]) != sp))

              # If there are co-occurring species, compute mean distance
              if (length(cooccur_species) > 0) {
                distance_sorted <- sort(dist_matrix[sp, cooccur_species], decreasing = FALSE)
                if(nearest.taxon == "all"){ # calculating for all taxon
                  mean_distances[sp, 1] <- mean(as.numeric(dist_matrix[sp, cooccur_species]), na.rm = TRUE)
                } else{ # using the threshold distance set by the user
                  mean_distances[sp, 1] <- mean(distance_sorted[1:nearest.taxon], na.rm = TRUE)
                }
              }
            }
          }
        }
        return(mean_distances)
      })

    mean_dist_timeslice <- lapply(list_dist_spp, function(x) mean(x, na.rm = TRUE))
    var_dist_timeslice <- lapply(list_dist_spp, function(x) var(x, na.rm = TRUE))

    df_res <-
      data.frame(mean.distance = unlist(mean_dist_timeslice),
                 var.distance = unlist(var_dist_timeslice),
                 time.slice = seq_interval)

    return(df_res)

  }
