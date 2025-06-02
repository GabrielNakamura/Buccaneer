
dados_caninae
trait <- dados_caninae[, c("mass", "rbl")]
dist_body_mass <- dist(trait)
df.TS.TE = dados_caninae
time.slice = 0.1
dist.trait = dist_body_mass
nearest.taxon = "all"
trait = NULL
round.digits = 1
species = "spp"
TS = "TS"
TE = "TE"
group = "group"
group.focal.compare = c("focal_group", "focal_group")
type.comparison = "within"

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
    seq_interval <- seq(from = max(df.TS.TE[, "TS"]), to = min(df.TS.TE[, "TE"]), by = -time.slice)
    seq_interval <- round(seq_interval, digits = round.digits)

    # co-occurrence matrix
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = 1,
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
      spp_focal <- df.TS.TE[which(df.TS.TE$group == focal), "species"]$species
      spp_compare <- df.TS.TE[which(df.TS.TE$group == compare), "species"]$species

      if(type.comparison == "between"){# comparison between species of two groups
        matrix_dist_trait_comp <- matrix_dist_trait[spp_focal, spp_compare] # focal speices in lines and comparison in columns
      }
      if(type.comparison == "within"){ # comparison only within the focal group
        matrix_dist_trait_comp <- matrix_dist_trait[spp_focal, spp_focal]
      }
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
              }
              if(is.numeric(nearest.taxon) == TRUE){ # used to compute mean distances considering thresholds
                if(nearest.taxon <= dim(matrix_dist_comp3)[1]){
                  mean_dist_timeslice[i] <- mean(as.vector(matrix_dist_comp3[1:nearest.taxon, ]))
                  var_dist_timeslice[i] <- var(as.vector(matrix_dist_comp3[1:nearest.taxon, ]))
                } else{
                  nearest.taxon <- dim(matrix_dist_comp3)[1] # n neighbours
                  mean_dist_timeslice[i] <- mean(as.vector(matrix_dist_comp3[1:nearest.taxon, ]))
                  var_dist_timeslice[i] <- var(as.vector(matrix_dist_comp3[1:nearest.taxon, ]))
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
                 time.slice = paste("slice", seq_interval, sep = "_")
      )


    # data frame with results
    return(df_res)

  }
