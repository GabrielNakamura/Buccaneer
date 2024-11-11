#' Auxiliar function to calculate reach criteria 
#'
#' @param df.coords 
#'
#' @return
#' @export
#'
#' @examples
calc_reach <- 
  function(df.coords){
    
    # Using this data frame to calculate reach for each time interval
    names_slice <- unique(df.coords$time.slice)
    
    # Removing non numeric values
    df_long_slice_coord_spp2 <- 
      df.coords %>%
      filter(if_any(c(lat, lon), ~ !is.na(as.numeric(.)))) |> 
      mutate(lat = as.numeric(lat), lon = as.numeric(lon))
    
    # subsetting coordinate data frame to accommodate each bin at each element in a list
    list_coords_slice <- 
      lapply(names_slice, function(x){
        pos <- which(x == df_long_slice_coord_spp2$time.slice)
        df_long_slice_coord_spp2[pos, ]
      })
    names(list_coords_slice) <- names_slice
    
    # getting the name of unique species occurring in each bin 
    list_spp_slice <- lapply(list_coords_slice, function(x) unique(x$species))
    list_spp_slice_all <- lapply(list_coords_slice, function(x) x$species)
    
    
    # transforming coordinates in sf objects to calculate spherical distances 
    list_sf_coords <- 
      lapply(list_coords_slice, function(x){
        x %>% 
          sf::st_as_sf(coords = c("lon", "lat"))
      })
    
    # calculating distance matrix 
    matrix_all_dist <- 
      lapply(list_sf_coords, function(x){
        sf::st_distance(x)/1000
      })
    
    # naming columns and rows of all distance matrix
    for(i in 1:length(matrix_all_dist)){
      colnames(matrix_all_dist[[i]]) <- list_spp_slice_all[[i]]
      rownames(matrix_all_dist[[i]]) <- list_spp_slice_all[[i]]
    }
    
    # getting combination of species
    combination <- 
      lapply(list_sf_coords, function(x){
        if(length(unique(x$species)) == 1){
          df <- NA
        } else{
          df <- data.frame(combn(unique(x$species), m = 2, simplify = TRUE))
        }
        return(df)
      })
    
    # removing NAs - no co-occurrence, only one species
    combination2 <- combination[-which(is.na(combination) == TRUE)]
    matrix_all_dist2 <- matrix_all_dist[-which(is.na(combination) == TRUE)]
    
    # filtering the combinations and filling the co-occurrence matrix with zeroes and 1 accordingly to the reach criteria 
    list_res <- vector(mode = "list", length = length(matrix_all_dist2))
    for(i in 1:length(matrix_all_dist2)){
      # i = 5
      res <- matrix(0, nrow = length(unique(rownames(matrix_all_dist2[[i]]))), ncol = length(unique(colnames(matrix_all_dist2[[i]]))), 
                    dimnames = list(unique(rownames(matrix_all_dist2[[i]])), 
                                    unique(colnames(matrix_all_dist2[[i]])))
      )
      for(j in 1:ncol(combination2[[i]])){
        # j = 5
        
        min_between <- matrix_all_dist2[[i]][which(combination2[[i]][1, j] == rownames(matrix_all_dist2[[i]])), 
                                            which(combination2[[i]][2, j] == colnames(matrix_all_dist2[[i]]))]
        max_1 <- matrix_all_dist2[[i]][which(combination2[[i]][1, j] == rownames(matrix_all_dist2[[i]])), 
                                      which(combination2[[i]][1, j] == colnames(matrix_all_dist2[[i]]))]
        max_2 <- matrix_all_dist2[[i]][which(combination2[[i]][2, j] == rownames(matrix_all_dist2[[i]])), 
                                      which(combination2[[i]][2, j] == colnames(matrix_all_dist2[[i]]))]
        res[combination2[[i]][1, j], combination2[[i]][2, j]] <- ifelse(min(min_between) < sum(max(max_1), max(max_2)), 1, 0)
        
      }
      list_res[[i]] <- res
    }
    
    # naming coexistence matrices
    names_slice2 <- names_slice[-which(is.na(combination) == TRUE)]
    names(list_res) <- names_slice2
    
    return(list_res)
    
  }



