#' Calculates regional clade richness for lineages
#'
#' @param df.TS.TE Data frame object containing at least four columns. Species names,
#'     origination time, extinction time and a trait value for each species.
#' @param time.slice Scalar indicating the time interval between consecutive time slices.
#' @param round.digits Scalar indicating the precision of time slices.
#' @param species Character indicating the name of the column of the data frame
#'     containing the species name information.
#' @param TS Character indicating the name of the columns of the data frame
#'     containing the information on origination time.
#' @param TE Character indicating the name of the column of the data frame
#'     containing the information on extinction time.
#'
#' @return A data frame with richness information in each timeslice
#' @export
#'
#' @examples
clade_regional_richness <-
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
    df_richness_slice <- data.frame(richness = unlist(nspp_slice), time.slice = seq_interval)

    return(df_richness_slice)

  }
