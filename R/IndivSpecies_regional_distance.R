#' Calculate Mean Trait Distance for Individual Species Co-occurring at Regional Scale
#'
#' This function computes mean pairwise trait distances between species that
#' co-occur in regional scale across different time slices. For each
#' species present at a time, it calculates the mean distance to its co-occurring
#' species, supporting both mean nearest neighbor distance (MNND) and mean
#' pairwise distance (MPD) metrics. The function can perform comparisons between
#' groups or within a single group.
#'
#' @param df.TS.TE A data frame containing species temporal and, optionally,
#'     trait data with at least three columns: species names, origination times (TS), extinction
#'     times (TE), and, optionally, trait values.
#'     Additional columns may include group assignments.
#' @param time.slice Numeric. The time interval (in the same units as TS and TE)
#'     between consecutive time slices for binning occurrences.
#' @param dist.trait A distance matrix object (class \code{dist} or \code{matrix})
#'     containing pairwise trait distances between species. Row and column names
#'     must match species names in \code{df.TS.TE}. If NULL, distances will be
#'     computed from the trait column in \code{df.TS.TE} using Euclidean distance.
#' @param nearest.taxon Numeric or character. The number of nearest neighbors to
#'     consider when calculating mean distances. Use \code{1} for mean nearest
#'     neighbor distance (MNND), or \code{"all"} for mean pairwise distance (MPD).
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
#' @param trait Character. The name of the column in \code{df.TS.TE} containing
#'     trait values. If NULL (default), \code{dist.trait} must be provided.
#' @param round.digits Integer. The number of decimal places to round time slice
#'     values. Default is 1. This affects temporal binning precision.
#' @param species Character. The name of the column in \code{df.TS.TE} and
#'     \code{df.occ} containing species identifiers. Default is "species".
#' @param TS Character. The name of the column in \code{df.TS.TE} containing
#'     origination times for each species. Default is "TS".
#' @param TE Character. The name of the column in \code{df.TS.TE} containing
#'     extinction times for each species. Default is "TE".
#'
#' @return A data frame with four columns:
#'   \item{species}{Character. The name of each species.}
#'   \item{time.slice}{Numeric. The time slice identifier.}
#'   \item{mean_dist_to_cooccur}{Numeric. The mean trait distance from each
#'       species to its co-occurring species at each site and time slice.
#'       Returns 0 for singleton species (species with no co-occurring taxa),
#'       and NA when data is insufficient.}
#'   \item{is_singleton}{Logical. TRUE if the species occurred alone at the
#'       site during that time slice, FALSE otherwise.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Divides the temporal range into discrete time slices based on \code{time.slice}
#'   \item Determines which species co-occur at each site within each time slice
#'   \item For each species at each site, calculates mean trait distances to co-occurring species
#'   \item Optionally filters comparisons by group membership (focal vs. comparison groups)
#' }
#'
#' Species-level results include:
#' \itemize{
#'   \item \strong{Singleton species}: Species occurring alone at a site receive a
#'         \code{mean_dist_to_cooccur} value of 0 and \code{is_singleton = TRUE}
#'   \item \strong{Missing data}: Time slices or sites with insufficient data return NA
#'   \item \strong{Group filtering}: When using \code{group.focal.compare}, only
#'         distances between specified groups are computed
#' }
#'
#' @export
#'
#' @examples
#' # Create example fossil data
#' df_longevities <- data.frame(
#'   species = c("sp1", "sp2", "sp3", "sp4"),
#'   TS = c(100, 98, 98, 85),
#'   TE = c(50, 45, 40, 35),
#'   trait = c(1.2, 2.5, 3.1, 4.0),
#'   group = c("A", "A", "B", "B")
#' )
#'
#' # Calculate MPD for all species at each site
#' result <- IndivSpec_regional_distance(
#'   df.TS.TE = df_longevities,
#'   trait = "trait",
#'   time.slice = 5,
#'   dist.trait = NULL,
#'   nearest.taxon = "all"
#' )
#'
#' # Calculate MNND between groups at each site
#' result_between <- IndivSpec_regional_distance(
#'   df.TS.TE = df_longevities,
#'   time.slice = 5,
#'   trait = "trait",
#'   dist.trait = NULL,
#'   nearest.taxon = "all",
#'   group = "group",
#'   group.focal.compare = c("A", "B"),
#'   type.comparison = "between"
#' )
#'
#'
#' result_within <- IndivSpec_regional_distance(
#'   df.TS.TE = df_longevities,
#'   time.slice = 5,
#'   trait = "trait",
#'   dist.trait = NULL,
#'   nearest.taxon = "all",
#'   group = "group",
#'   group.focal.compare = c("A", "B"),
#'   type.comparison = "within"
#' )
IndivSpec_regional_distance <-
  function(df.TS.TE,
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
           TE = "TE"){
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

    # Generating time intervals used to compute temporal coexistence
    seq_interval <- seq(from = ceiling(max(df.TS.TE[, "TS"])),
                        to = ceiling(min(df.TS.TE[, "TE"])),
                        by = -time.slice)

    # regional coexistence matrix based on longevities
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE,
                               time.slice,
                               round.digits = round.digits,
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
    # if trait is not informed in longevities data frame use trait distance
    if(!is.null(dist.trait) == TRUE){
      matrix_dist_trait <- as.matrix(dist.trait)
    } else{ # otherwise use the trait provided to create the distance matrix
      matrix_dist_trait <-
        as.matrix(dist(x = df.TS.TE[, "trait"], method = "euclidean", upper = T, diag = T))
    }
    # matching name orders in longevities and trait matrix
    rownames(matrix_dist_trait) <- df.TS.TE$species
    colnames(matrix_dist_trait) <- df.TS.TE$species
    matrix_dist_trait <-
      matrix_dist_trait[match(rownames(matrix_dist_trait), df.TS.TE$species),
                        match(colnames(matrix_dist_trait), df.TS.TE$species)]



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

    # Ensure same species and same order in both cooccurrence matrix and distance matrix
    list_dist_spp <-
      lapply(matrix_coex, function(x){
        species_row <- intersect(rownames(x), rownames(matrix_dist_trait_comp))
        species_col <- intersect(colnames(x), colnames(matrix_dist_trait_comp))
        cooccur_matrix <- x[species_row, species_col, drop = FALSE]
        dist_matrix <- matrix_dist_trait_comp[species_row, species_col, drop = FALSE]

        # Create matrix to store the results of mean pairwise distances
        mean_distances <-
          matrix(NA, nrow = length(species_row), ncol = 1,
                 dimnames = list(species_row,
                                 paste("mean.dist.to.cooccur",
                                       nearest.taxon, sep = ".")
                 )
          )

        # calculating distances for all species
        for (sp in species_row) {
          #checking if there are no species in the slice or
          if(length(species_row) == 0 | length(species_col) == 0){
            mean_distances[sp, 1] <- NA # no cooccurrence in sites
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

    names(list_dist_spp) <- format(seq_interval, trim = TRUE, scientific = FALSE)

    df_dist_spp <-
      do.call(rbind, lapply(names(list_dist_spp), function(age) {
        element <- list_dist_spp[[age]]

        if ((is.matrix(element) || is.data.frame(element)) && nrow(element) > 0) {
          data.frame(
            species = rownames(element),
            time.slice = as.numeric(age),
            mean_dist_to_cooccur = element[, 1],
            row.names = NULL
          )
        } else {
          # Return NA row if element is not a matrix/data.frame or has 0 rows
          data.frame(
            species = NA,
            time.slice = as.numeric(age),
            mean_dist_to_cooccur = NA
          )
        }
      }))

    # Ensure the column is character first
    df_dist_spp$mean_dist_to_cooccur <- as.character(df_dist_spp$mean_dist_to_cooccur)

    # Create flag column for singleton entries
    df_dist_spp$is_singleton <- df_dist_spp$mean_dist_to_cooccur == "NA_singleton"

    # Replace "NA_singleton" with "0"
    df_dist_spp$mean_dist_to_cooccur[df_dist_spp$mean_dist_to_cooccur == "NA_singleton"] <- "0"

    # Replace "<NA>" strings with real NA
    df_dist_spp$mean_dist_to_cooccur[df_dist_spp$mean_dist_to_cooccur == "<NA>"] <- NA

    # Convert mean_dist_to_cooccur to numeric
    df_dist_spp$mean_dist_to_cooccur <- as.numeric(df_dist_spp$mean_dist_to_cooccur)


    return(df_dist_spp)

  }
