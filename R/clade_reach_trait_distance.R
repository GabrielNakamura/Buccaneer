#' Calculate Mean Clade-Level Trait Distances Across Time Slices With Reach Criteria
#'
#' @param df.TS.TE
#' @param df.occ
#' @param time.slice
#' @param round.digits
#' @param species
#' @param TS
#' @param TE
#' @param lat
#' @param lon
#' @param Max.age
#' @param Min.age
#' @param crs
#' @param remove.singletons
#' @param group
#' @param group.focal.compare
#' @param type.comparison
#'
#' @returns
#' @export
#'
#' @examples
clade_reach_distance <-
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
           remove.singletons = TRUE,
           group = NULL,
           group.focal.compare = NULL,
           type.comparison = NULL){

    # renaming longevity data frame
    if(!is.null(group) == TRUE){
      df.TS.TE <- df.TS.TE[, c(species, TS, TE, group)]
      colnames(df.TS.TE) <- c("species", "TS", "TE", "group")
    } else{
      df.TS.TE <- df.TS.TE[, c(species, TS, TE)]
      colnames(df.TS.TE) <- c("species", "TS", "TE")
    }

    # renaming occurrence data frame
    df_occ <-
      df.occ[, c(species, lat, lon, Max.age, Min.age)]
    vars <- list(species, lat, lon, Max.age, Min.age)
    name_vars <- c("species", "lat", "lon", "Max.age", "Min.age")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names

    # time slices
    names_slice <- seq(from = ceiling(max(df.TS.TE[, "TS"])),
                       to = ceiling(min(df.TS.TE[, "TE"])),
                       by = -time.slice)

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


    # Time coexistence matrix for all species
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE,
                               time.slice = time.slice,
                               round.digits = round.digits,
                               species = "species",
                               TS = "TS",
                               TE = "TE")

    # modified coexistence matrix containing group comparison, in rows are focal species, in colums comparison
    # the "within" argument place focal species in both columns and rows
    if(!is.null(group.focal.compare) == TRUE){
      focal <- group.focal.compare[1]
      compare <- group.focal.compare[2]
      spp_focal <- df.TS.TE[which(df.TS.TE$group == focal), "species"]
      spp_compare <- df.TS.TE[which(df.TS.TE$group == compare), "species"]

      if(type.comparison == "between"){# comparison between species of two groups
        matrix_coex <- lapply(matrix_coex, function(x) x[spp_focal, spp_compare]) # focal species in lines and comparison in columns
      }
      if(type.comparison == "within"){ # comparison only within the focal group
        matrix_coex <- lapply(matrix_coex, function(x) x[spp_focal, spp_focal])
      }
    } else{
      matrix_coex <- matrix_coex
    }

    # species composition at each time slice
    spp_slice <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    # naming list with time slices
    names(spp_slice) <- format(names_slice, trim = TRUE, scientific = FALSE)

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

    # obtaining occurrence records for each time slice based on species composition in time slices
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

    matrix_all_dist2 <- list_matrix_dist_occ

    # filtering the combinations and filling the co-occurrence matrix with zeroes and 1s accordingly with reach criteria
    list_res <- vector(mode = "list", length = length(matrix_coex))
    for(i in 1:length(matrix_all_dist2)){
      # i = 100
      res <- matrix_coex[[i]] # regional cooccurrence matrix
      if(is.data.frame(combination[[i]]) != TRUE){
        res <- NA
      } else{
        for(j in 1:ncol(combination[[i]])){

          min_between <-
            matrix_all_dist2[[i]][which(combination[[i]][1, j] ==
                                          rownames(matrix_all_dist2[[i]])),
                                  which(combination[[i]][2, j] ==
                                          colnames(matrix_all_dist2[[i]])),
                                  drop = FALSE]
          max_1 <-
            matrix_all_dist2[[i]][which(combination[[i]][1, j] ==
                                          rownames(matrix_all_dist2[[i]])),
                                  which(combination[[i]][1, j] ==
                                          colnames(matrix_all_dist2[[i]])),
                                  drop = FALSE]
          max_2 <-
            matrix_all_dist2[[i]][which(combination[[i]][2, j] ==
                                          rownames(matrix_all_dist2[[i]])),
                                  which(combination[[i]][2, j] ==
                                          colnames(matrix_all_dist2[[i]])),
                                  drop = FALSE]
          res[combination[[i]][1, j], combination[[i]][2, j]] <-
            ifelse(min(min_between) <=
                     sum(max(max_1),
                         max(max_2)), 1, 0)
        }
      }
      list_res[[i]] <- res
    }

    # naming coexistence matrices
    names(list_res) <- format(names_slice, trim = TRUE, scientific = FALSE)

    # Ensure same species and same order in both cooccurrence matrix and distance matrix
    list_dist_spp <-
      lapply(list_res, function(x){
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

    # calculating mean distances and variance
    mean_dist_timeslice <- lapply(list_dist_spp, function(x) mean(x, na.rm = TRUE))
    var_dist_timeslice <- lapply(list_dist_spp, function(x) var(x, na.rm = TRUE))

    # wraping the results in a matrix
    df_res <-
      data.frame(mean.distance = unlist(mean_dist_timeslice),
                 var.distance = unlist(var_dist_timeslice),
                 time.slice = names_slice)

    return(df_res)

  }

