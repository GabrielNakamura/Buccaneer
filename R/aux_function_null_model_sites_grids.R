
#' Compute null model with fixed richness for sites and assemblage
#'
#' @param list_occurrence a data frame containing species occurrence in sites
#'     or assemblages/grids
#' @param dist_matrix_trait a distance matrix computed
#' @param nperm the number of permutations for the null model
#'
#' @returns a data frame with distance metric for all time slices
#' @export
#'
#' @examples
calc_null_model <-
  function(list_occurrence, dist_matrix_trait, nearest.taxon, nperm = 1000){

    # data frame to receive the results
    df_null_res <- matrix(NA, nrow = nperm, ncol = length(list_occurrence),
                          dimnames = list(paste("perm", 1:nperm, sep = "_"),
                                          paste("slice", 1:length(list_occurrence), sep = "_")
                          )
    )

    for(i in 1:length(list_occurrence)){
      # i = 310
      # compute null matrix for each timeslice composition
      null_occ_site <-
        vegan::permatfull(m = list_occurrence[[i]][, -1],
                          fixedmar = "rows",
                          mtype = "prab",
                          times = nperm)

      # computing species co-occurrence matrix
      list_matrix_coocccur_null_site <-
        lapply(null_occ_site$perm, function(x){
          mat <- as.matrix(x)  # remove site column
          cooccur_matrix_sites <- t(mat) %*% mat  # species-by-species co-occurrence
          cooccur_matrix_sites
        })

      # matching species co-occurrence matrix and species distance matrix
      list_dist_matrix_null <-
        lapply(list_matrix_coocccur_null_site, function(x){
          species_row <- intersect(rownames(x), rownames(dist_matrix_trait))
          species_col <- intersect(colnames(x), colnames(dist_matrix_trait))
          dist_matrix <- dist_matrix_trait[species_row, species_col]
          dist_mat_filter <- ifelse(x == 1, dist_matrix, NA)
          dist_mat_filter
        })

      list_dist_matrix_null2 <-
        lapply(list_dist_matrix_null, function(x){
          diag(x) <- NA
          x
        })

      # in case of nearest taxon
      if(is.numeric(nearest.taxon)){
        # ordering the null distance matrix
        list_order_null <-
          lapply(list_dist_matrix_null2, function(x){
            apply(x, 1, simplify = TRUE, function(y) sort(y))
          })

        # this is for nearest taxon
        list_mean_nearest <-
          lapply(list_order_null, function(x){
            mean(unlist(lapply(x, function(x) x[1:nearest.taxon])), na.rm = T)
          })

        # unlisting and storing the null values
        df_null_res[, i] <- unlist(list_mean_nearest)
      } else{ # this is the case when the user set trait distance as all


        # removing lower triangle of the matrix to get rid of the duplicated values
        list_dist_lower_matrix_null <-
          lapply(list_dist_matrix_null, function(x){
            x[lower.tri(x)]
          })

        vec_null_mpd_slice <-
          unlist(lapply(list_dist_lower_matrix_null, function(x){
            mean(x, na.rm = TRUE)
          }))

        # mean_null_slice <- mean(vec_null_mpd_slice)
        df_null_res[, i] <- vec_null_mpd_slice
      }
    } # end loop for each timeslice
    return(df_null_res)
  }


