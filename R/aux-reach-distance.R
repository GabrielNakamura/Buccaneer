library(dplyr)
df.coords = df_can
df.TS.TE = longs[[1]]
df.coords <- 
  df.coords %>% 
  select(sp.name, lng, lat, midpoint) %>% 
  rename(species = sp.name, lon = lng, lat = lat, age = midpoint)

df.TS.TE <- 
  df.TS.TE %>% 
  mutate(species = rownames(df.TS.TE), TS = TS, TE = TE)

intervals = dog_NALMAs_age

source(here::here("R", "gap-filling.R"))
source(here::here("R", "temporal-coexist.R"))

#' Auxiliar function to compute reach distance coexistence matrix
#'
#' @param df.coords A four column data frame with coordinates of fossil occurrences (latitude and longitude) for each species and their respective age 
#' @param df.TS.TE A three column data frame with species' longevities. One column contains the name of species, one TS and one TE
#' @param intervals A numeric vector containing time intervals 
#' @param type.fill Character indicating which type of method must be used to calculate geographical gaps in fossil records between adjacent intervals
#'
#' @return A list with pairwise matrix containing the coexistence of each species in each time interval (each element of the matrix)
#' @export
#'
#' @examples
calc_reach <- 
  function(df.coords, df.TS.TE, intervals, type.fill = "centroid"){
    
    # fill all the coordinate gaps in time intervals
    df_gap_fill <- make_gap_filling(df.coords = df.coords,
                                    df.TS.TE = df.TS.TE,
                                    intervals = intervals,
                                    type.fill = "centroid")
    
    # Using this data frame to calculate reach for each time interval
    names_bin <- gtools::mixedsort(names(table(df_gap_fill$bin)), decreasing = TRUE) 
    
    # subsetting coordinate data frame to accommodate each bin at each element in a list
    list_coords_bin <- 
      lapply(names_bin, function(x){
      pos <- which(x == df_gap_fill$bin)
      df_gap_fill[pos, ]
    })
    names(list_coords_bin) <- names_bin
    
    # getting the name of unique species occurring in each bin 
    list_spp_bin <- lapply(list_coords_bin, function(x) unique(x$species))
    list_spp_bin_all <- lapply(list_coords_bin, function(x) x$species)
    
    
    # transforming coordinates in sf objects to calculate spherical distances 
    list_sf_coords <- 
      lapply(list_coords_bin, function(x){
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
      colnames(matrix_all_dist[[i]]) <- list_spp_bin_all[[i]]
      rownames(matrix_all_dist[[i]]) <- list_spp_bin_all[[i]]
    }
    
    # getting combination of species
    combination <- 
      lapply(list_sf_coords, function(x){
        data.frame(combn(unique(x$species), m = 2, simplify = TRUE))
      })
    
    # filtering the combinations and filling the co-occurrence matrix with zeroes and 1 accordingly to the reach criteria 
    list_res <- vector(mode = "list", length = length(matrix_all_dist))
    for(i in 1:length(matrix_all_dist)){
      # i = 5
      res <- matrix(0, nrow = length(unique(rownames(matrix_all_dist[[i]]))), ncol = length(unique(colnames(matrix_all_dist[[i]]))), 
                    dimnames = list(unique(rownames(matrix_all_dist[[i]])), 
                                    unique(colnames(matrix_all_dist[[i]])))
      )
      for(j in 1:ncol(combination[[i]])){
        # j = 5
        
        min_between <- matrix_all_dist[[i]][which(combination[[i]][1, j] == rownames(matrix_all_dist[[i]])), 
                                            which(combination[[i]][2, j] == colnames(matrix_all_dist[[i]]))]
        max_1 <- matrix_all_dist[[i]][which(combination[[i]][1, j] == rownames(matrix_all_dist[[i]])), 
                                      which(combination[[i]][1, j] == colnames(matrix_all_dist[[i]]))]
        max_2 <- matrix_all_dist[[i]][which(combination[[i]][2, j] == rownames(matrix_all_dist[[i]])), 
                                      which(combination[[i]][2, j] == colnames(matrix_all_dist[[i]]))]
        res[combination[[i]][1, j], combination[[i]][2, j]] <- ifelse(min(min_between) < sum(max(max_1), max(max_2)), 1, 0)
        
      }
      list_res[[i]] <- res
    }
    names(list_res) <- names_bin
    
    return(list_res)
    
  }



