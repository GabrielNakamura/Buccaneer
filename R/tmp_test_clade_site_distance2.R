
df.TS.TE = df_TS_TE_mass
df.occ = df_occ_faurby
time.slice = 0.1
round.digits = 1
species = "species"
TS = "TS"
TE = "TE"
dist.trait = dist_body_mass
nearest.taxon = 1
trait = NULL
group = "group"
group.focal.compare = c("Caniformia", "Feliformia")
type.comparison = "between"

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
    seq_interval <- seq(from = max(df.TS.TE[, "TS"]), to = min(df.TS.TE[, "TE"]), by = -time.slice)
    seq_interval <- c(round(seq_interval, digits = round.digits))

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

    names(spp_slice) <- seq_interval

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
      spp_focal <- df.TS.TE[which(df.TS.TE$group == focal), "species"]$species
      spp_compare <- df.TS.TE[which(df.TS.TE$group == compare), "species"]$species

      if(type.comparison == "between"){# comparison between species of two groups
        matrix_dist_trait_comp <- matrix_dist_trait[spp_focal, spp_compare] # focal speices in lines and comparison in columns
      }
      if(type.comparison == "within"){ # comparison only within the focal group
        matrix_dist_trait_comp <- matrix_dist_trait[spp_focal, spp_focal]
      }
    }


    # calculating the number of species per site in each time slice

    list_site_interval <-
      lapply(1:length(spp_slice), function(x){
        names_slices <- names(spp_slice)[[x]]
        names_species <- spp_slice[[x]]
        filtered_df <-
          df.occ %>% filter(species %in% names_species)

        filtered_df_site <-
          filtered_df |>
          filter(Max.age >= names_slices & Min.age <= names_slices)

        spp_coex_site <-
          filtered_df_site |>
          distinct(species, .keep_all = TRUE) |>
          select(species, TS, TE)

        filtered_df_site2 <-
          filtered_df_site |>
          group_by(site) |>
          distinct(species, .keep_all = T) |>
          add_count(site, name = "species.per.site") |>
          select(species, site, TS, TE, species.per.site)

        return(filtered_df_site2)
      })

    # list of cooccurrences filtered by site occurrences
    list_site_cooccur <-
      lapply(list_site_interval, function(x){
        calc_cooccur_site_matrix(x)
      })

    # trait matrix filtered by site
    list_trait_site_cooccur <-
      lapply(list_site_cooccur, function(x){
        matrix_dist_trait_comp[rownames(matrix_dist_trait_comp) %in% rownames(x),
                               colnames(matrix_dist_trait_comp) %in% colnames(x)]
      })

    # getting all distances for all species from focal group coexisting with species from comparison group
    lapply(list_site_cooccur, function(x){
      matrix_dist_trait_comp[1,
                             colnames(matrix_dist_trait_comp) %in% names(which(x[1, ] == 1))]
    })
    vec_dist_trait_site <-
      matrix_dist_trait_comp[rownames(list_site_cooccur[[550]])[1],
                             colnames(matrix_dist_trait_comp) %in% names(which(list_site_cooccur[[550]][1, ] == 1))] # take only the comparison group species


    # calculating mnnd for all timeslices based on site cooccurrence
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
            if(type.comparison == "within"){
              if(is.vector(matrix_dist_comp3) == TRUE){
                matrix_dist_comp3 <- matrix_dist_comp3
              } else{
                matrix_dist_comp3 <- matrix_dist_comp3[-1,]
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

    return(df_mean_rich_site)

  }


# matrix co-occurrence in timeslices (W)

matrix_coex <-
  aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = 1,
                           species = "species",
                           TS = "TS",
                           TE = "TE")


# matrix trait distance species (T)

if(!is.null(dist.trait) == TRUE){
  matrix_dist_trait <- as.matrix(dist.trait)
} else{
  matrix_dist_trait <-
    as.matrix(dist(x = df.TS.TE[, "trait"], method = "euclidean", upper = T, diag = T))
}
rownames(matrix_dist_trait) <- df.TS.TE$species
colnames(matrix_dist_trait) <- df.TS.TE$species

# matrix co-occurrence in sites (S)

list_site_cooccur <-
  lapply(list_site_interval, function(x){
    calc_cooccur_site_matrix(x, spp_slice = )
  })

# joint matrix trait (T) x timeslice (W) (TW)


# joint matrix TW x S
