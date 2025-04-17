#' Compute regional clade Mean Trait Distance in a time series
#'
#' @param df.TS.TE Data frame object containing at least four columns. Species names,
#'     origination time, extinction time and a trait value for each species.
#' @param time.slice Scalar indicating the time interval between consecutive time slices.
#' @param trait Character indicating the name of the column containing values of the traits for each species. If NULL,
#'     the default, the user must provide a distance matrix.
#' @param round.digits Scalar indicating the precision of time slices.
#' @param species Character indicating the name of the column of the data frame
#'     containing the species name information.
#' @param TS Character indicating the name of the columns of the data frame
#'     containing the information on origination time.
#' @param TE Character indicating the name of the column of the data frame
#'     containing the information on extinction time.
#' @param group Character indicating the name of the column that contain the groups that will be used
#'     in comparison.
#' @param group.focal.compare Character vector indicating the focal (first element) and comparison (second element).
#' @param type.comparison Character. It can be "between" to compute distances only between species/genus of two groups
#'     or "within" to calculate distance only inside the focal group.
#' @param dist.trait A dist object containing the pairwise distance matrices
#' @param nearest.taxon A scalar indicating the number of nearest species/genus that will be used.
#'     1 computes mnnd metric and the option "all" computes mpd.
#'
#' @return
#' @export
#'
#' @examples
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

    comm_all <- matrix(0, nrow = length(spp_slice), ncol = length(df.TS.TE$species),
                       dimnames = list(paste("slice", seq_interval, sep = "_"), df.TS.TE$species))

    # community matrix
    for(i in 1:length(spp_slice)){
      # i = 552
      pos_spp <- match(spp_slice[[i]], colnames(comm_all))
      pos_comm <- i
      comm_all[i, pos_spp] <- 1
    }

    # trait distance
    if(!is.null(dist.trait) == TRUE){
      matrix_dist_trait <- as.matrix(dist.trait)
    } else{
      matrix_dist_trait <-
        as.matrix(dist(x = df.TS.TE[, "trait"], method = "euclidean", upper = T, diag = T))
    }
    rownames(matrix_dist_trait) <- df.TS.TE$species
    colnames(matrix_dist_trait) <- df.TS.TE$species

    # calculating the distances by each timeslice according to at threshold
    matrix_dist_trait_times <-
      aux_mpd_general(df.TS.TE = df.TS.TE,
                      dist.trait = dist.trait,
                      nearest.taxon = nearest.taxon)


    # modified trait matrix containing group comparison
    if(!is.null(group.focal.compare) == TRUE){
      focal <- group.focal.compare[1]
      compare <- group.focal.compare[2]
      spp_focal <- df.TS.TE[which(df.TS.TE$group == focal), "species"]$species
      spp_compare <- df.TS.TE[which(df.TS.TE$group == compare), "species"]$species

      if(type.comparison == "between"){# comparison between species of two groups
        matrix_dist_trait[spp_focal, spp_focal] <- 0 # removing within groups distances
        matrix_dist_trait[spp_compare, spp_compare] <- 0 # removing within groups distances
      }
      if(type.comparison == "within"){ # comparison only within the focal group
        spp_focal <- df.TS.TE[which(df.TS.TE$group == focal), "species"]$species
        spp_remove <- df.TS.TE[which(focal != df.TS.TE$group), "species"]$species
        matrix_dist_trait[spp_remove, spp_remove] <- 0
        matrix_dist_trait[spp_remove, spp_focal] <- 0
        matrix_dist_trait[spp_focal, spp_remove] <- 0
      }
    }



    # calculating mpd for all species considering group comparisons and otu co-occurrences in each timeslice
    res_vec_mean_dist <- numeric(length = length(matrix_dist_trait_times))
    res_vec_sd_dist <- numeric(length = length(matrix_dist_trait_times))
    for(i in 1:length(matrix_dist_trait_times)){
      # i = 211
      if(all(is.na(matrix_dist_trait_times[[i]])) == TRUE){
        res_vec_mean_dist[i] <- NA
      } else{
        matrix_dist_trait_time <- matrix_dist_trait[match(names(matrix_dist_trait_times[[i]]), colnames(matrix_dist_trait)),
                                                    match(names(matrix_dist_trait_times[[i]]), colnames(matrix_dist_trait))]
        # matrix with only distances from desired group comparison and otu that occur in a given time slice
        spp_rows_keep <- spp_focal[match(rownames(matrix_dist_trait_time), spp_focal)]
        spp_columns_keep <- spp_compare[match(rownames(matrix_dist_trait_time), spp_compare)]
        if(all(is.na(spp_rows_keep)) == TRUE){
          res_vec_mean_dist[i] <- NA
        } else{
          mat_comp_group_time <-
            matrix_dist_trait_time[na.omit(spp_rows_keep),
                                   na.omit(spp_columns_keep)]
          res_vec_mean_dist[i] <- mean(as.vector(mat_comp_group_time))
          res_vec_sd_dist[i] <- sd(as.vector(mat_comp_group_time))
        }
      }
    }

    df_res <-
      data.frame(mean.distance = res_vec_mean_dist,
                 sd.distance = res_vec_sd_dist,
                 time.slice = paste("slice", seq_interval, sep = "_")
      )


    # data frame with results
    return(df_res)

  }
