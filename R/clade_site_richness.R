#' Compute time series based on mean site richness 
#'
#' @param df.TS.TE 
#' @param df.occ 
#' @param time.slice 
#' @param round.digits 
#' @param species 
#' @param TS 
#' @param TE 
#'
#' @return
#' @export
#'
#' @examples
clade_site_richness <- 
  function(df.TS.TE,
           df.occ,
           time.slice, 
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE"){
    df.TS.TE <- 
      df.TS.TE[, c(species, TS, TE)] 
    vars <- list(species, TS, TE)
    name_vars <- c("species", "TS", "TE")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df.TS.TE) <- column.names
    
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
    remove_no_occ <- which(is.na(df_sub_slice2) == TRUE)
    df_sub_slice3 <- df_sub_slice2[-remove_no_occ]
    species_slices <- lapply(df_sub_slice3, function(x) x$species)
    
    # calculating the number of species per site in each time slice
    
    list_site_interval <- 
      lapply(1:length(species_slices), function(x){
      names_slices <- names(species_slices)[[x]]
      names_species <- species_slices[[x]]
      filtered_df <- 
        df.occ %>% filter(species %in% names_species)
      
      filtered_df_site <-  
        filtered_df |> 
        filter(Max.age >= names_slices & Min.age <= names_slices)
      
      spp_coex_site <- 
        filtered_df_site |> 
        distinct(species, .keep_all = TRUE) |> 
        select(species, TS, TE)
      
      filtered_df_site2 <- 
        filtered_df_site |> 
        group_by(site) |> 
        distinct(species, .keep_all = T) |> 
        add_count(site, name = "species.per.site") |> 
        select(species, site, TS, TE, species.per.site)
      
      filtered_df_site3 <- 
        filtered_df_site2 |> 
        ungroup(site) |> 
        distinct(site, .keep_all = T) |> 
        select(site, species.per.site)
      
      filtered_df_site4 <- 
        filtered_df_site3 |> 
        mutate(time.slice = as.character(names_slices))
      return(filtered_df_site4)
    })
    
    
    df_long_slice_site <- do.call(rbind, list_site_interval)
    
    df_mean_rich_site <- 
      df_long_slice_site |> 
      group_by(time.slice) |> 
      mutate(mean.richness.site = mean(species.per.site), 
             variance.richness.site = var(species.per.site)) |> 
      mutate(n.site.slice = n()) |> 
      distinct(time.slice, .keep_all = TRUE) |> 
      select(time.slice, mean.richness.site, variance.richness.site, n.site.slice)
    
    return(df_mean_rich_site)
    
  }