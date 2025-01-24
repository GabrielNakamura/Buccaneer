#' Compute regional time series of Individual Species Mean Pairwise trait distance
#'
#' @param df.TS.TE a data frame object containing at least four columns. Species names,
#'     origination time, extinction time and a trait value for each species.
#' @param time.slice scalar indicating the time interval between consecutive time slices
#' @param trait Numeric indicating the values of the traits for each species
#' @param dist.trait A distance matrix containing pairwise distance between species
#' @param nearest.taxon scalar indicating the number of closest species that will be used to compute
#'     the mean pairwise distance. If "all", the default, all species coexisting in a time slice
#'     will be used to calculate mpd
#' @param round.digits scalar indicating the number of digits for time of origination and time for
#'     extinction
#' @param species character indicating the name of the column of the data frame
#'     containing the species name information
#' @param TS character indicating the name of the columns of the data frame
#'     containing the information on origination time
#' @param TE character indicating the name of the column of the data frame
#'     containing the information on extinction time
#'
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
           trait = NULL,
           dist.trait = NULL,
           nearest.taxon = "all",
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE"){

    # subsetting columns with TS, TE, trait and species
    df.TS.TE <- df.TS.TE[, c(species, trait, TS, TE)]
    if(is.null(trait) == TRUE){
      colnames(df.TS.TE) <- c("species", "TS", "TE")
    } else{
      colnames(df.TS.TE) <- c("species", "trait", "TS", "TE")
    }


    # creating time slices
    seq_interval <- c(seq(from = max(df.TS.TE[, "TS"]), to = min(df.TS.TE[, "TE"]), by = -time.slice), 0)
    seq_interval <- round(seq_interval, digits = round.digits)

    # co-occurrence matrix
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = 1,
                               species = "species",
                               TS = "TS",
                               TE = "TE")


    if(is.null(dist.trait) != TRUE){
      if(!inherits(dist.trait, "dist")){
        stop("dist.trait must be a distance object")
      }
      matrix_dist_trait <- as.matrix(dist.trait)
      rownames(matrix_dist_trait) <- df.TS.TE$species
      colnames(matrix_dist_trait) <- df.TS.TE$species

    } else{
      matrix_dist_trait <-
        as.matrix(dist(x = df.TS.TE[, "trait"], method = "euclidean", upper = T, diag = T))

      rownames(matrix_dist_trait) <- df.TS.TE$species
      colnames(matrix_dist_trait) <- df.TS.TE$species
    }

    # calculating mpd for all species
    spp_coex_names <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    list_subset_coex <-
      lapply(seq_along(spp_coex_names), function(x){
        matrix_dist_trait[spp_coex_names[[x]], spp_coex_names[[x]]]
      })
    # sorting distance matrix
    list_subset_coex <-
      lapply(list_subset_coex,
             function(x){
               if(!is.matrix(x) == TRUE){
                 NA
               } else{
                 apply(x, 1,
                       function(y) sort(y, decreasing = FALSE)
                 )
               }
             }
      )

    # removing self-distances between species
    list_subset_coex <-
      lapply(list_subset_coex,
             function(x){
               if(!is.matrix(x) == TRUE){
                 NA
               } else{
                 x[-1, ]
               }
             }
      )

    # for nearest taxon
    if(is.numeric(nearest.taxon) == TRUE){
      # ordering the matrix based on trait distances
      list_subset_coex <-
        lapply(list_subset_coex, function(x){
          if(!is.matrix(x) == TRUE){
            NA
          } else{
            x[1:nearest.taxon, 1:ncol(x)]
          }
        })

    }

    # mean mpd per species
    mean_per_spp <-
      lapply(list_subset_coex, function(x){
        if(!is.matrix(x) == TRUE & !is.vector(x) == TRUE){
          NA
        } else{
          if(is.vector(x) == TRUE){
            x
          } else{
            colSums(x)/(nrow(x))
          }
        }
      })


    # organizing data frame with results

    # naming data frames
    names(mean_per_spp) <- seq_interval
    mean_per_spp2 <-
      lapply(names(mean_per_spp), function(name) {
        if(length(mean_per_spp[[name]]) == 0){
          df <- NA
        } else{
          df <- mean_per_spp[[name]]
          df$time.slice <- name  # Add the name of the element as the new column
        }
        return(df)
      })

    df_IndivSpecies_mpd <-
      lapply(seq_along(mean_per_spp), function(x){
        data.frame(mpd = mean_per_spp[[x]],
                   slice = rep(seq_interval[[x]], length(mean_per_spp[[x]])))
      })

    df_IndivSpecies_mpd2 <-
      lapply(df_IndivSpecies_mpd, function(x){
        data.frame(spp = rownames(x), mpd = x$mpd, time.slice = x$slice)
      })

    df_IndivSpecies_mpd3 <- do.call(rbind, df_IndivSpecies_mpd2)

    df_IndivSpecies_mpd4 <-
      df_IndivSpecies_mpd3 |>
      filter(spp != 1)

    # data frame with results
    return(df_IndivSpecies_mpd4)

  }
