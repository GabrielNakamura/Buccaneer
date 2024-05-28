#' Filling geographical gaps in fossil occurrence through time
#'
#' @param df.coords A data frame containing the species occurrence records. This must have a column 
#'     with the name of the species (species), the longitude (lon), latitude (lat) and the time in which the 
#'     register was sampled (age).
#' @param df.TS.TE A data frame containing the origination (TS), extinction (TE) times and the name of the species
#' @param intervals A numeric vector containing the time intervals in which the fossil occurrence register will be analyzed
#' @param type.fill Character indicating the method used to fill the geographical gap. This can be one of three options: 
#'     \item{centroid}: a new coordinate is created based on the centroid of all fossil occurrences available for the species
#'     \item{random}: a new coordinate is created based on a random sampling of a occurrence record
#'     \item{adjacent}: a new coordinate is created based on the closest occurrence record
#' @param range Logical. If TRUE the method will compute the occurrences in a bin based on the age range of fossil
#'     record
#' @param species Character indicating the name of df.coords column containing the species names 
#' @param lon Character indicating the name of the df.coords column containing the longitude
#' @param lat Character indicating the name of the df.coords column containing the latitude 
#' @param age Character indicating the name of the column of df.coords containing the age when the fossil
#'     record was sampled
#' @param min.age Character indicating the name of the column of df.coords containing the minimum range 
#'     estimate for the age of fossil occurrence record
#' @param max.age Character indicating the name of the column of df.coords containing the minimum range 
#'     estimate for the age of fossil occurrence record
#'
#' @return A data frame with six columns containing:
#'     \item{species}: column containing the name of the species corresponding to the fossil record
#'     \item{lon}: column containing the longitude of the fossil record
#'     \item{lat}: column containing the latitude of the fossil record
#'     \item{age}: column containing the age of the fossil record
#'     \item{bin}: column containing the time bin in which the fossil record was sampled
#'     \item{status}: column containing the information about the inputation of the record. occ correspond 
#'         to fossils that present true occurrence information. inputed correspond to occurrences that 
#'         were inputed based one of the methods
#' @export
#'
#' @examples
make_gap_filling <- 
  function(df.coords, 
           df.TS.TE,
           intervals,
           type.fill = c("centroid", "random", "adjacent"),
           range = FALSE, 
           species = "species",
           lon = "lon",
           lat = "lat",
           age = "age",
           min.age = "min.age",
           max.age = "max.age"){
    
    
    
    # fossil occ records in each bin
    list_spp <- vector(mode = "list", length = (length(intervals) - 1))
    for(i in 1:(length(intervals) - 1)){
      # i = 10
      list_spp[[i]] <- df.coords[(df.coords$age <= intervals[i]) & (df.coords$age >= intervals[i + 1]), ]
      list_spp[[i]]$bin <- paste(as.character(intervals[i]), as.character(intervals[i + 1]), sep = "-")
    }
    
    df.coords2 <- do.call(rbind, list_spp) 
    
    # species in each bin
    list_spp_fossil <- lapply(list_spp, function(x) unique(x$species))
    
    # just giving names to the bins
    names_interval <- vector(length = (length(intervals) - 1))
    for(i in 1:(length(intervals) - 1)){
      # i = 10
      names_interval[i] <- paste(as.character(intervals[i]), as.character(intervals[i + 1]), sep = "-")
    }
    names(list_spp_fossil) <- names_interval
    
    # calculating temporal coexistence to get species in each bin from TS and TE
    temp_coex <- calc_temp_only_coex(x = df.TS.TE, intervals = dog_NALMAs_age)
    
    # species in each bin based on TS and TE
    list_spp_long <-
      lapply(unique(temp_coex$long.time.occurence$time.interval), function(x){
        as.vector(temp_coex$long.time.occurence[which(x == temp_coex$long.time.occurence$time.interval), "species"]$species)
      })
    names(list_spp_long) <- names_interval
    
    # list of species with no fossil occ records
    list_no_fossil <- 
      lapply(1:length(list_spp_long), function(x){
        list_spp_long[[x]][which(is.na(match(list_spp_long[[x]], list_spp_fossil[[x]])) == TRUE)]
      })
    names(list_no_fossil) <- names_interval
    
    # in case of centroid method, calculate the centroid for the species with no fossil in a bin
    df_centroid_sf <- 
      lapply(unique(unlist(list_no_fossil)), function(x){
        df_coords <- df.coords[which(x == df.coords$species), ]
        df_sf <- sf::st_as_sf(df_coords, coords = c(2, 3))
        df_sf2 <- sf::st_set_crs(df_sf, 3857)
        df_sf3 <- sf::st_combine(df_sf2)
        df_sf_centroid <- sf::st_centroid(df_sf3)
        return(df_sf_centroid)
      })
    names(df_centroid_sf) <- unique(unlist(list_no_fossil))
    
    # list with data frames containing inputed coordinates for each species in each bin
    list_inputation_coords <- vector(mode = "list", length = length(names(df_centroid_sf)))
    df_list_no_fossil <- data.frame(species = unlist(list_no_fossil), bin = rep(names(list_no_fossil), sapply(list_no_fossil, length)))
    for(i in 1:length(names(df_centroid_sf))){
      r1 <- which(names(df_centroid_sf)[[i]] == df_list_no_fossil$species)
      which(names(df_centroid_sf)[[i]] == list_no_fossil)
      df_inputation_coords <- 
        data.frame(species = names(df_centroid_sf)[[i]],
                   lon = unlist(df_centroid_sf[[i]])[1],
                   lat = unlist(df_centroid_sf[[i]])[2], 
                   bin = df_list_no_fossil$bin[r1])
      list_inputation_coords[[i]] <- df_inputation_coords
    }
    
    df_inputation_coords2 <- do.call(rbind, list_inputation_coords)
    df_inputation_coords2$status <- "geo_input"
    
    # joining data frame with inputation with occ data frame
    df_all_coords <- dplyr::bind_rows(df.coords2, df_inputation_coords2)
    
    # identifying in final data frame the species that have been inputed
    df_coords_input <- 
      df_all_coords %>% 
      mutate(status = ifelse(is.na(status) == TRUE, "occ", status))
    
    # final table organized by bins
    names_bin <- unique(df_coords_input$bin)
    df_coords_input2 <- 
      df_coords_input %>% 
      group_by(bin) %>% 
      arrange(factor(bin, levels = gtools::mixedsort(names_bin, decreasing = T)))
    
    return(df_coords_input2)
    
  }