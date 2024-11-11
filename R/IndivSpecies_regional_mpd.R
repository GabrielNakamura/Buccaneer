#' Compute regional time series of Individual Species Mean Pairwise trait distance
#'
#' @param df.TS.TE a data frame object containing at least four columns. Species names,
#'     origination time, extinction time and a trait value for each species.
#' @param time.slice scalar indicating the time interval between consecutive time slices
#' @param trait Numeric indicating the values of the traits for each species
#' @param round.digits scalar indicating the number of digits for time of origination and time for
#'     extinction
#' @param species character indicating the name of the column of the data frame
#'     containing the species name information
#' @param TS character indicating the name of the columns of the data frame
#'     containing the information on origination time
#' @param TE character indicating the name of the column of the data frame
#'     containing the information on extinction time
#'
#' @return A data frame containing the name of species and the mean value of mpd of each species
#'     in each time slice
#'
#' @export
#'
#' @examples
IndivSpecies_regional_mpd <-
  function(df.TS.TE,
           time.slice,
           trait,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE"){
    # subsetting columns with TS, TE, trait and species
    df.TS.TE <- df.TS.TE[, c(species, trait, TS, TE)]
    colnames(df.TS.TE) <- c("species", "trait", "TS", "TE")

    # creating time slices
    seq_interval <- seq(from = max(df.TS.TE[, "TS"]), to = min(df.TS.TE[, "TE"]), by = -time.slice)
    seq_interval <- round(seq_interval, digits = round.digits)

    # co-occurrence matrix
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = 1,
                               species = "species",
                               TS = "TS",
                               TE = "TE")


    # trait distance
    matrix_dist_trait <-
      as.matrix(dist(x = df.TS.TE[, "trait"], method = "euclidean", upper = T, diag = T))

    rownames(matrix_dist_trait) <- df.TS.TE$species
    colnames(matrix_dist_trait) <- df.TS.TE$species

    # calculating mpd for all species
    spp_coex_names <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    list_subset_coex <-
      lapply(seq_along(spp_coex_names), function(x){
        matrix_dist_trait[spp_coex_names[[x]], spp_coex_names[[x]]]
      })

    # mean mpd per species
    mean_per_spp <-
      lapply(list_subset_coex, function(x){
        if(!is.matrix(x) == TRUE){
          NA
        } else{
          rowSums(x)/(ncol(x) - 1)
        }
      })


    # organizing data frame with results

    df_IndivSpecies_mpd <-
      lapply(seq_along(mean_per_spp), function(x){
        data.frame(mpd = mean_per_spp[[x]],
                   slice = rep(seq_interval[[x]], length(mean_per_spp[[x]])))
      })

    df_IndivSpecies_mpd2 <-
      lapply(df_IndivSpecies_mpd, function(x){
        data.frame(species = rownames(x), mpd = x$mpd, time.slice = x$slice)
      })

    df_IndivSpecies_mpd3 <- do.call(rbind, df_IndivSpecies_mpd2)

    # mean coexistence per species
    df_mean_individual_coex <-
      df_IndivSpecies_mpd3 |>
      dplyr::group_by(species) |>
      dplyr::summarise(mean_mpd = mean(mpd))

    # list with results
    list_res <- vector(mode = "list")
    list_res$df_IndivSpp_mpd <- df_IndivSpecies_mpd3
    list_res$mean_species_mpd <- df_mean_individual_coex

    # data frame with results
    return(list_res)


  }
