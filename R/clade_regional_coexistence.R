#' Calculates regional clade coexistence for lineages
#'
#' @param df.TS.TE Data frame object containing at least three columns. Species names,
#'     origination time, extinction time.
#' @param time.slice Scalar indicating the time interval between consecutive
#'     time slices in which coexistence will be calculated.
#' @param round.digits Scalar indicating the precision of time slices. This corresponds
#'     to the number of digits to be considered in the calculation of time
#'     coexistence in each time slice.
#' @param species Character indicating the name of the column of the data frame
#'     containing the species name information. The default is "species".
#' @param TS Character indicating the name of the columns of the data frame
#'     containing the information on origination time. The default is "TS".
#' @param TE Character indicating the name of the column of the data frame
#'     containing the information on extinction time. The default is "TE".
#'
#' @return A two column data frame object containing the number of coexisting
#'     clades (coexistence) at each time slice (time.slice).
#' @export
#'
#' @examples
#' data(longevities_canidae) # list containing longevities for canidae species
#' longs_1 <- data.frame(longs[[1]], species = rownames(longs[[1]])) # extracting only one replicate
#' res_regional_coexistence <- clade_regional_coexistence(df.TS.TE = longs_1,
#'                                                        time.slice,
#'                                                        round.digits = 1,
#'                                                        species = "species",
#'                                                        TS = "TS",
#'                                                        TE = "TE")
#'
#'
clade_regional_coexistence <-
  function(df.TS.TE,
           time.slice,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE"){


    df.TS.TE <- df.TS.TE[, c(species, TS, TE)]
    colnames(df.TS.TE) <- c("species", "TS", "TE")

    # creating time slices
    seq_interval <- seq(from = ceiling(max(df.TS.TE[, "TS"])),
                        to = ceiling(min(df.TS.TE[, "TE"])),
                        by = -time.slice)

    # co-occurrence matrix
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = round.digits,
                               species = "species",
                               TS = "TS",
                               TE = "TE")

    # community composition matrix

    spp_slice <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    nspp_slice <- lapply(spp_slice, function(x) length(x))

    # dataframe with richness
    df_richness_slice <- data.frame(coexistence = unlist(nspp_slice), time.slice = seq_interval)

    return(df_richness_slice)

  }
