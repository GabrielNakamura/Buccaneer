#' Compute regional clade Mean Pairwise Distance in a time series
#'
#' @param df.TS.TE
#' @param time.slice
#' @param trait
#' @param round.digits
#' @param species
#' @param TS
#' @param TE
#' @param compute.ses.mpd
#' @param null.model
#' @param runs
#' @param group
#' @param group.focal.compare
#' @param type.comparison
#'
#' @return
#' @export
#'
#' @examples
clade_regional_mpd <-
  function(df.TS.TE,
           time.slice,
           trait,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE",
           compute.ses.mpd = TRUE,
           null.model = "taxa.labels",
           runs = 999,
           group = NULL,
           group.focal.compare = NULL,
           type.comparison = NULL
  ){

    # subseting columns

    if(!is.null(group) == TRUE){
      df.TS.TE <- df.TS.TE[, c(species, trait, TS, TE, group)]
      colnames(df.TS.TE) <- c("species", "trait", "TS", "TE", "group")
    } else{
      df.TS.TE <- df.TS.TE[, c(species, trait, TS, TE)]
      colnames(df.TS.TE) <- c("species", "trait", "TS", "TE")
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

    comm_all <- matrix(0, nrow = length(seq_interval), ncol = length(df.TS.TE$species),
                       dimnames = list(paste("slice", seq_interval, sep = "_"), df.TS.TE$species))

    spp_slice <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    # community matrix
    for(i in 1:length(spp_slice)){
      pos_spp <- match(spp_slice[[i]], colnames(comm_all))
      pos_comm <- i
      comm_all[i, pos_spp] <- 1
    }

    # trait distance
    matrix_dist_trait <-
      as.matrix(dist(x = df.TS.TE[, "trait"], method = "euclidean", upper = T, diag = T))

    rownames(matrix_dist_trait) <- df.TS.TE$species
    colnames(matrix_dist_trait) <- df.TS.TE$species


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



    # calculating mpd for all species

    # raw value
    res_mpd_regional <- picante::ses.mpd(samp = comm_all, dis = matrix_dist_trait, abundance.weighted = FALSE, null.model = null.model, runs = 1)
    res_mpd_regional2 <-
      res_mpd_regional |>
      dplyr::select(ntaxa, mpd.obs) |>
      dplyr::mutate(richness = ntaxa, mpd = mpd.obs, time.slice = seq_interval) |>
      dplyr::select(richness, mpd, time.slice)

    # effect size
    if(compute.ses.mpd == TRUE){
      res_mpd_regional <- picante::ses.mpd(samp = comm_all,
                                           dis = matrix_dist_trait,
                                           null.model = null.model,
                                           runs = runs)
      res_mpd_regional2 <-
        res_mpd_regional |>
        dplyr::select(ntaxa, mpd.obs, mpd.obs.z, mpd.obs.p) |>
        dplyr::mutate(richness = ntaxa,
                      ses.mpd = mpd.obs.z,
                      p.value = mpd.obs.p,
                      time.slice = seq_interval) |>
        dplyr::select(richness, ses.mpd, p.value, time.slice)

    }

    # data frame with results
    return(res_mpd_regional2)


  }
