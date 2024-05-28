
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

calc_site <- 
  function(df.coords, df.TS.TE, intervals, site.name = FALSE){
    
   
    
    if(site.name == TRUE){
      
    }
    
    # first identifying sites based on their lon and lat
    if(site.name == FALSE){
      
      
      # using gap filling just to get the bins
      df_bin <- 
        make_gap_filling(df.coords = df.coords,
                         df.TS.TE = df.TS.TE,
                         intervals = intervals,
                         type.fill = "centroid")
      
      df_bin2 <- 
        df_bin %>% 
        filter(status == "occ") %>% 
        select(species, lon, lat, age, bin)
      
      # Using this data frame to calculate reach for each time interval
      names_bin <- gtools::mixedsort(names(table(df_bin2$bin)), decreasing = TRUE) 
      
      # subsetting coordinate data frame to accommodate each bin at each element in a list
      list_coords_bin <- 
        lapply(names_bin, function(x){
          pos <- which(x == df_bin2$bin)
          df_bin2[pos, ]
        })
      names(list_coords_bin) <- names_bin
      
      # getting the name of unique species occurring in each bin 
      list_spp_bin <- lapply(list_coords_bin, function(x) unique(x$species))
      list_spp_bin_all <- lapply(list_coords_bin, function(x) x$species)
      
      for(i in 1:length(list_spp_bin)){
        # i = 3
        for(j in 1:nrow(list_coords_bin[[i]])){
          # j = 1
          df1 <- 
            list_coords_bin[[i]] %>% 
            filter(species == "Archaeocyon_leptodus")
          df2 <- 
            list_coords_bin[[i]] %>% 
            filter(species != "Archaeocyon_leptodus")
          generics::intersect(df1, df2)
          which(df1$lon[1] == df2$lon & df1$lat[1] == df2$lat)
          df2[30, ]  
            
          all_equal(df1[, c("lon", "lat")], list_coords_bin[[i]][, c("lon", "lat")])
          which(list_coords_bin[[i]][j, c("lon", "lat")] == list_coords_bin[[i]][, c("lon", "lat")])
        }
      }
      
      
    }
    
    
    
  }