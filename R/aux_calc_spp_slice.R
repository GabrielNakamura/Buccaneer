
#' Auxiliar function to compute the names of species in each timeslice
#'
#' @param df.TS.TE
#' @param time.slice
#' @param round.digits
#' @param species
#' @param TS
#' @param TE
#'
#' @returns a list with names of species in each timeslice
#' @export
#'
#' @examples
calc_spp_slice <-
  function(df.TS.TE,
           time.slice,
           seq.interval,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE"){
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

    names(spp_slice) <- seq.interval
    return(spp_slice)
  }
