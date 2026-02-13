#' Calculate Regional Individual Species Coexistence Across Time Slices
#'
#' This function computes the number of coexisting species at each time slice
#' in a regional scale analysis. It divides the temporal range into discrete
#' intervals and counts how many species overlap in each interval based on
#' their origination (TS) and extinction (TE) times.
#'
#' @param df.TS.TE A data frame containing species temporal data with at least
#'     three columns: species names, origination times, and extinction times.
#' @param time.slice Numeric. The time interval (in the same units as TS and TE)
#'     between consecutive time slices.
#' @param round.digits Integer. The number of decimal places to round time slice
#'     values. Default is 1. This affects the precision of temporal binning.
#' @param species Character. The name of the column in \code{df.TS.TE} containing
#'     species identifiers. Default is "species".
#' @param TS Character. The name of the column in \code{df.TS.TE} containing
#'     origination (start) times for each species. Default is "TS".
#' @param TE Character. The name of the column in \code{df.TS.TE} containing
#'     extinction (end) times for each species. Default is "TE".
#'
#' @return A data frame with three columns:
#'   \item{time.slice}{Numeric. The time point representing each slice.}
#'   \item{species}{Character. The name of each species.}
#'   \item{n.coexistence}{Numeric. The number of species coexisting in each time
#'     slice for all species.}
#'
#'
#' @export
#'
#' @examples
#' # Example with fossil data
#' df <- data.frame(
#'   species = c("sp1", "sp2", "sp3"),
#'   TS = c(100, 90, 80),
#'   TE = c(50, 40, 30)
#' )
#' IndivSpec_regional_coex(df, time.slice = 10)
IndivSpec_regional_coex <-
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

    # number of coexistences per species

    list_n_coex_all <- lapply(matrix_coex, function(x) rowSums(x))
    list_n_coex_all2 <- lapply(list_n_coex_all, function(x) ifelse(x == 0, NA, x))
    list_n_coex_all3 <- lapply(list_n_coex_all2, function(x) x[-which(is.na(x) == TRUE)])
    list_n_coex_all4 <- lapply(list_n_coex_all3, function(x) x - 1) # removing self coexistence
    # naming list with individual species coexistence
    names(list_n_coex_all4) <- format(seq_interval, trim = TRUE, scientific = FALSE)
    list_n_coex_all5 <- list_n_coex_all4[ lengths(list_n_coex_all4) > 0 ]

    # dataframe with coexistence
    df_long_res <-
      do.call(rbind, lapply(names(list_n_coex_all5), function(nm) {
        data.frame(
          time.slice = as.numeric(nm),
          species = names(list_n_coex_all5[[nm]]),
          n.coexistence = unname(list_n_coex_all5[[nm]])
        )
      })
      )

    return(df_long_res)

  }
