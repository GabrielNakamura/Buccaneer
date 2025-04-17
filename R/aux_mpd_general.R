#' Auxiliar function to calculate trait distance considering different distance thresholds
#'
#' @param df.TS.TE
#' @param dist.trait
#' @param nearest.taxon
#'
#' @returns
#' @export
#'
#' @examples
aux_mpd_general <-
  function(df.TS.TE, dist.trait, nearest.taxon){
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice = 0.1, round.digits = 1,
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
    return(list_subset_coex) # distance matrix for each timeslice
  }
