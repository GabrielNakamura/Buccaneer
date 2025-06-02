#' Calculate individual species coexistence in regional scale
#'
#' @param df.TS.TE a data frame object containing at least three columns. Species names,
#'     origination time and extinction time for each species.
#' @param time.slice scalar indicating the time interval between consecutive time slices
#' @param round.digits scalar indicating the number of digits for time of origination and time for
#'     extinction
#' @param species character indicating the name of the column of the data frame
#'     containing the species name information
#' @param TS character indicating the name of the columns of the data frame
#'     containing the information on origination time
#' @param TE character indicating the name of the column of the data frame
#'     containing the information on extinction time
#'
#' @return A list with two elements:
#'     \item{df_IndivSpp_coexist}{A data frame with species name, the number of coexistence for each sepecies in each time slice}
#'     \item{mean_species_coexistence}{A data frame with three column. The species name, the mean number of coexistence for each species
#'         and the variance considering all time slices}
#' @export
#'
#' @examples
IndivSpecies_regional_coex <-
  function(df.TS.TE,
           time.slice,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE"){
    # reorganizing data frame TS and TE
    df.TS.TE <- df.TS.TE[, c(species, TS, TE)]

    # Generating time intervals used to compute temporal coexistence
    seq_interval <- seq(from = max(df.TS.TE[, "TS"]), to = min(df.TS.TE[, "TE"]), by = -time.slice)
    seq_interval <- round(seq_interval, digits = round.digits)

    # coexistence per time slice
    list_matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = 1,
                               species = "species",
                               TS = "TS",
                               TE = "TE")


    # calculating individual species metric
    individual_coex <-
      lapply(list_matrix_coex, function(x){
        richness <- rowSums(x)
        indiv_coex <- length(richness[which(richness >= 1)]) - 1 # removing the self coexistence
        spp_names <- names(which(rowSums(x) >= 1))
        df_indiv_coex <- data.frame(species = spp_names,
                                    n.coexistence = rep(indiv_coex, times = length(spp_names)))
        return(df_indiv_coex)
      })

    # naming data frames
    names(individual_coex) <- seq_interval
    individual_coex2 <-
      lapply(names(individual_coex), function(name) {
        if(nrow(individual_coex[[name]]) == 0){
          df <- NA
        } else{
          df <- individual_coex[[name]]
          df$time.slice <- name  # Add the name of the element as the new column
        }
        return(df)
      })

    df_long_indiv <- do.call(rbind, individual_coex2)

    # removing NAs
    df_long_indiv2 <-
      df_long_indiv |>
      tidyr::drop_na() |>
      dplyr::mutate(time.slice = as.numeric(time.slice))

    # mean coexistence per species
    df_mean_individual_coex <-
      df_long_indiv2 |>
      dplyr::group_by(species) |>
      dplyr::summarise(mean_coex = mean(n.coexistence))

    list_res <- vector(mode = "list")
    list_res$df_IndivSpp_coexist <- df_long_indiv2
    list_res$mean_species_coexistence <- df_mean_individual_coex
    return(list_res)
  }
