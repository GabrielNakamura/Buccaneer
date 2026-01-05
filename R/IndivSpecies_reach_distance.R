#' Calculate Mean Individual Species Trait Distances Across Time Slices With Reach Criteria
#'
#' @param df.TS.TE A data frame containing species temporal data with at least
#'     three columns: species names, origination times (TS), and extinction
#'     times (TE). Additional columns may include group assignments.
#' @param df.occ  data frame containing fossil occurrence records with at least
#'     five columns: species names, minimum age, maximum age, and latitude and
#'     longitude.
#'     Each row represents a single occurrence record at a specific site.
#' @param time.slice Numeric. The time interval (in the same units as TS and TE)
#'     between consecutive time slices for temporal binning.
#' @param round.digits Integer. The number of decimal places to round time slice
#'     values. Default is 1. This affects temporal binning precision.
#' @param species Character. The name of the column in \code{df.TS.TE} and
#'     \code{df.occ} containing species identifiers. Default is "species".
#' @param TS Character. The name of the column in \code{df.TS.TE} containing
#'     origination (first appearance) times for each species. Default is "TS".
#' @param TE Character. The name of the column in \code{df.TS.TE} containing
#'     extinction (last appearance) times for each species. Default is "TE".
#' @param lat Numeric. The latitude coordinate of the occurrence record in
#'     \code{df.occ}.
#' @param lon Numeric. The longitude coordinate of the occurrence records in
#'     \code{df.occ}.
#' @param Max.age Character. The name of the column in \code{df.occ} containing
#'     the maximum (oldest) age estimate for each occurrence record. Default is "Max.age".
#' @param Min.age Character. The name of the column in \code{df.occ} containing
#'     the minimum (youngest) age estimate for each occurrence record. Default is "Min.age".
#' @param crs Numeric. The code indicating the coordinate reference system
#'     to be used for latitude and longitude of occurrence records in \code{df.occ}
#' @param remove.singletons Logical. Should singleton species (species occurring
#'     alone at a site with no co-occurring species) be excluded from mean and
#'     variance calculations? Default is TRUE. When TRUE, singletons are treated
#'     as NA; when FALSE, they contribute 0 to the mean.
#' @param group Character. The name of the column in \code{df.TS.TE} containing
#'     group assignments for species (e.g., clade, family). Required if using
#'     \code{group.focal.compare}. Default is NULL.
#' @param group.focal.compare Character vector of length 2. The first element
#'     specifies the focal group and the second specifies the comparison group.
#'     If NULL (default), coexistence is calculated across all species regardless
#'     of group membership.
#' @param type.comparison Character. Specifies the type of coexistence comparison:
#'     \itemize{
#'       \item \code{"between"}: Count only co-occurrences between species from
#'             the focal and comparison groups.
#'       \item \code{"within"}: Count only co-occurrences among species within
#'             the focal group.
#'       \item NULL (default): Count all co-occurrences regardless of group.
#'     }
#'
#' @returns A data frame with three columns:
#'   \item{mean.distance}{Numeric. The mean distance for each individual species in
#'       each time slice for all other cooccurring species according with
#'       reach criteria. This represents average distance of each species in the
#'       morphospace relative to all other cooccurring species
#'       (excluding or including singletons based on \code{remove.singletons}).}
#'   \item{var.distance}{Numeric. The variance in the mean distance for each
#'       individual species in each time slice.}
#'   \item{time.slice}{Numeric. The time point representing each slice, typically
#'       the upper (older) boundary of the time bin.}
#'
#' @export
#'
#' @examples
IndivSpec_reach_distance <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           dist.trait,
           nearest.taxon,
           round.digits,
           species,
           TS,
           TE,
           lat,
           lon,
           Max.age = "Max.age",
           Min.age = "Min.age",
           crs = 4326,
           group = NULL,
           group.focal.compare = NULL,
           type.comparison = NULL){
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
      df.occ[, c(species, lat, lon, Max.age, Min.age)]
    vars <- list(species, lat, lon, Max.age, Min.age)
    name_vars <- c("species", "lat", "lon", "Max.age", "Min.age")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names

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

    # calculating matrix of species cooccurrence using reach criteria
    # Removing non numeric values in lat long in occurrence data frame
    df_long_slice_coord_spp2 <-
      df_occ |>
      filter(if_any(c(lat, lon), ~ !is.na(as.numeric(.)))) |>
      mutate(lat = as.numeric(lat), lon = as.numeric(lon))

    # transforming to a coordinate system  in occurrence data frame
    df_long_slice_coord_spp3 <-
      sf::st_as_sf(x = df_long_slice_coord_spp2,
                   coords = c("lon", "lat"),
                   crs = crs)

    df_occ_occurr <-
      comp_slice_occ_cooccurr(spp_slice = spp_slice,
                              df.occ = df_long_slice_coord_spp3,
                              species = "species",
                              Max.age = "Max.age",
                              Min.age = "Min.age")

    # list with matrices containing geographical distances for all occurrence records in all time slices
    list_matrix_dist_occ <-
      lapply(df_occ_occurr, function(x){
        dist_occ <- sf::st_distance(x = x)
        colnames(dist_occ) <- x$species
        rownames(dist_occ) <- x$species
        return(dist_occ)
      })

    # getting pairwise combination between all species cooccurring in all time slices
    combination <-
      lapply(df_occ_occurr, function(x){
        spp_names <- unique(x$species)
        if(length(spp_names) == 1 | length(spp_names) == 0){
          df <- NA
        } else{
          df <- data.frame(combn(spp_names, m = 2, simplify = TRUE))
        }
        return(df)
      })

    # removing NAs - no co-occurrence, only one species occurrence in the timeslice
    combination2 <- combination[-which(is.na(combination) == TRUE)]
    matrix_all_dist2 <- list_matrix_dist_occ[-which(is.na(combination) == TRUE)]
    matrix_coex2 <- matrix_coex[-which(is.na(combination) == TRUE)]
    # filtering the combinations and filling the co-occurrence matrix with zeroes and 1 accordingly with reach criteria
    list_res <- vector(mode = "list", length = length(matrix_coex2))
    for(i in 1:length(matrix_all_dist2)){

      res <- matrix_coex2[[i]] # regional cooccurrence matrix
      for(j in 1:ncol(combination2[[i]])){
        # j = 5

        min_between <- matrix_all_dist2[[i]][which(combination2[[i]][1, j] == rownames(matrix_all_dist2[[i]])),
                                             which(combination2[[i]][2, j] == colnames(matrix_all_dist2[[i]])),
                                             drop = FALSE]
        max_1 <- matrix_all_dist2[[i]][which(combination2[[i]][1, j] == rownames(matrix_all_dist2[[i]])),
                                       which(combination2[[i]][1, j] == colnames(matrix_all_dist2[[i]])),
                                       drop = FALSE]
        max_2 <- matrix_all_dist2[[i]][which(combination2[[i]][2, j] == rownames(matrix_all_dist2[[i]])),
                                       which(combination2[[i]][2, j] == colnames(matrix_all_dist2[[i]])),
                                       drop = FALSE]
        res[combination2[[i]][1, j], combination2[[i]][2, j]] <- ifelse(min(min_between) <= sum(max(max_1), max(max_2)), 1, 0)

      }
      list_res[[i]] <- res
    }

    # subsetting focal and compare groups
    # modified coexistence matrix containing group comparison, in rows are focal species, in colums comparison
    # the "within" argument place focal species in both columns and rows
    if(!is.null(group.focal.compare) == TRUE){
      focal <- group.focal.compare[1]
      compare <- group.focal.compare[2]
      spp_focal <- df.TS.TE[which(df.TS.TE$group == focal), "species"]
      spp_compare <- df.TS.TE[which(df.TS.TE$group == compare), "species"]

      if(type.comparison == "between"){# comparison between species of two groups
        matrix_coex2 <- lapply(list_res, function(x){
          if(is.matrix(x) != TRUE){
            NA
          } else{
            x[spp_focal, spp_compare]
          }
        }) # focal species in lines and comparison in columns
      }
      if(type.comparison == "within"){ # comparison only within the focal group
        matrix_coex2 <- lapply(list_res, function(x){
          if(is.matrix(x) != TRUE){
            NA
          } else{
            x[spp_focal, spp_focal]
          }
        })
      }
    } else{
      matrix_coex2 <- matrix_coex
    }

    # naming coexistence matrices
    names(matrix_coex2) <- format(seq_interval, trim = TRUE, scientific = FALSE)

    # Ensuring same species and same order in both cooccurrence matrix and distance matrix
    list_dist_spp <-
      lapply(matrix_coex2, function(x){
        if(is.matrix(x) != TRUE){
          mean_distances <-  NA
        } else{
          species_row <- intersect(rownames(x), rownames(matrix_dist_trait_comp))
          species_col <- intersect(colnames(x), colnames(matrix_dist_trait_comp))
          cooccur_matrix <- x[species_row, species_col, drop = FALSE]
          dist_matrix <- matrix_dist_trait_comp[species_row, species_col, drop = FALSE]

          # Create matrix to store the results of mean pairwise distances
          mean_distances <-
            matrix(NA, nrow = length(species_row), ncol = 1,
                   dimnames = list(species_row, paste("mean.dist.to.cooccur",
                                                      nearest.taxon, sep = ".")))

          # calculating distances for all species accordingly to the criteria chosen
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
                cooccur_species <-
                  names(which(cooccur_matrix[sp, ] > 0 &
                                names(cooccur_matrix[sp, ]) != sp))

                # If there are co-occurring species, compute mean distance
                if (length(cooccur_species) > 0) {
                  distance_sorted <- sort(dist_matrix[sp, cooccur_species],
                                          decreasing = FALSE)
                  if(nearest.taxon == "all"){ # calculating for all taxon
                    mean_distances[sp, 1] <- mean(as.numeric(dist_matrix[sp, cooccur_species]),
                                                  na.rm = TRUE)
                  } else{ # using the threshold distance set by the user
                    mean_distances[sp, 1] <- mean(distance_sorted[1:nearest.taxon],
                                                  na.rm = TRUE)
                  }
                }
              }
            }
          }
        }
        return(mean_distances)
      })

    # naming list elements with time slice
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
