#' Mean richness by communities in a time series
#'
#' @param df.TS.TE 
#' @param df.occ 
#' @param time.slice 
#' @param grid.size 
#' @param round.digits 
#' @param species 
#' @param TS 
#' @param TE 
#' @param Max.age 
#' @param Min.age 
#' @param lat 
#' @param lon 
#' @param crs 
#'
#' @return
#' @export
#'
#' @examples
Assemblage_regional_richness <- 
  function(df.TS.TE,
           df.occ,
           time.slice, 
           grid.size,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE",
           Max.age = "Max.age", 
           Min.age = "Min.age",
           lat = "lat",
           lon = "lng",
           crs = 4326
  ){
    
    df.TS.TE <- 
      df.TS.TE[, c(species, TS, TE)] 
    vars <- list(species, TS, TE)
    name_vars <- c("species", "TS", "TE")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df.TS.TE) <- column.names
    
    # filtering data from occurrence data frame
    df_occ <- 
      df.occ[, c(species, Max.age, Min.age, lat, lon)]
    vars <- list(species, Max.age, Min.age, lat, lon)
    name_vars <- c("species", "Max.age", "Min.age", "lat", "lon")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names
    
    # converting lat lon to numeric
    
    df_occ2 <-
      df_occ |> 
      mutate(lat = as.numeric(lat), lon = as.numeric(lon)) |> 
      filter(!is.na(lat) & !is.na(lon))
    
    # transforming lat long to sf object
    
    sf_occ <- sf::st_as_sf(df_occ2, coords = c("lon", "lat"), crs = crs) # WGS84 CRS
    
    # Create the grid over the extent of the points
    grid <- sf::st_make_grid(sf_occ, cellsize = c(grid.size, grid.size), what = "polygons")
    sf_grid <- sf::st_sf(grid_id = 1:length(grid), geometry = grid)  # Assign IDs to grid cells
    
    # finding points in each grid
    
    df_grid_occ <- sf::st_join(sf_occ, sf_grid, join = sf::st_within)
    
    # summary information on grid
    
    df_mean_var_age_grid <- 
      df_grid_occ |> 
      dplyr::mutate(midpoint = (abs(Max.age + Min.age)/2)) |> 
      group_by(grid_id) |> 
      mutate(mean.age.grid = mean(midpoint), 
             var.age.grid = var(midpoint))
    
    
    # joining grid information with occurrence information
    
    sf_grid_occ <- 
      sf_grid |>
      sf::st_join(df_mean_var_age_grid)
    
    # summary calculation with grid
    
    # Generating time intervals used to compute temporal coexistence
    seq_interval <- seq(from = max(df.TS.TE[, "TS"]), to = min(df.TS.TE[, "TE"]), by = -time.slice)
    seq_interval <- c(round(seq_interval, digits = round.digits), 0)
    
    
    # defining species per bin and subsetting longevities data frame
    df_sub_slice <- vector(mode = "list", length = length(seq_interval))
    for (i in 1:length(seq_interval)){
      # i = 1
      df_sub_slice[[i]] <- df.TS.TE[which(df.TS.TE$TS >= seq_interval[i] & df.TS.TE$TE <= seq_interval[i]), ]
    }
    
    names(df_sub_slice) <- seq_interval
    df_sub_slice2 <- 
      lapply(names(df_sub_slice), function(name) {
        if(nrow(df_sub_slice[[name]]) == 0){
          df <- NA
        } else{
          df <- df_sub_slice[[name]]
          df$time.slice <- name  # Add the name of the element as the new column
        }
        return(df)
      })
    
    
    names(df_sub_slice2) <- seq_interval
    list_interval <- 
      lapply(as.character(seq_interval), function(x) {
        if(!is.data.frame(df_sub_slice2[[x]]) == TRUE){
          NA
        } else{
          sf_grid_occ[match(df_sub_slice2[[x]]$species, sf_grid_occ$species), ] 
        }
        
      })
    
    names(list_interval) <- seq_interval
    
    list_res_rich <- 
      lapply(1:length(list_interval), function(x){
        if(!is.data.frame(list_interval[[x]]) == TRUE){
          df_res <- NA
        } else{
          df_res <- 
            list_interval[[x]] |> 
            filter(Max.age >= seq_interval[x] & Min.age <= seq_interval[x]) |> 
            mutate(time.slice = seq_interval[x]) |> 
            group_by(grid_id.x) |> 
            add_count(grid_id.x, name = "rich.by.grid") |> 
            mutate(mean.rich.by.grid.slice = mean(rich.by.grid)) |> 
            as.data.frame()
        }
        return(df_res)
      })
    
    names(list_res_rich) <- seq_interval
    
    # binding list 
    
    
    list_df_res_rich2 <- list_df_res_rich[!is.na(list_res_rich)]
    
    df_long_rich_assemblage <- do.call(rbind, list_df_res_rich2)
    return(df_long_rich_assemblage)
    
    
    
  }