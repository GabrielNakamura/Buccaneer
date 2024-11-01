#' Compute mean pairwise distance for each species co-occurring in sites and time slices
#'
#' @param df.TS.TE a data frame object containing at least three columns. Species names, 
#'     origination time and extinction time for each species. 
#' @param df.occ  a data frame object containing the occurrence records for each species. 
#'     This must have at least a column indicating the name of species, its minimum and maximum age estimate,
#'     and its site location ID.
#' @param time.slice scalar indicating the time interval between consecutive time slices.
#' @param trait Numeric indicating the values of the traits for each species
#' @param round.digits scalar indicating the number of digits for time of origination and time for 
#'     extinction.
#' @param species character indicating the name of the column of the data frame 
#'     containing the species name information.
#' @param TS character indicating the name of the columns of the data frame
#'     containing the information on origination time.
#' @param TE character indicating the name of the column of the data frame 
#'     containing the information on extinction time.
#' @param Max.age character indicating the name of the column in df.occ containing the upper age limit for occurrence record.
#' @param Min.age character indicating the name of the column in df.occ containing the lower age limit for occurrence record.
#' @param site character indicating the name of the column in df.occ containing the information on site location.
#'
#' @return A data frame containing the name of species, its mean mpd value calculated by each time slice
#'     considering all sites in which the species occur in that time slice
#' 
#' @export
#'
#' @examples
IndivSpec_site_mpd <- 
  function(df.TS.TE,
           df.occ,
           time.slice, 
           trait,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE",
           Max.age = "Max.age", 
           Min.age = "Min.age", 
           site = "site"){
    
    # filtering data from longevities data frame
    df.TS.TE <- 
      df.TS.TE[, c(species, TS, TE, trait)] 
    vars <- list(species, TS, TE, trait)
    name_vars <- c("species", "TS", "TE", "trait")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df.TS.TE) <- column.names
    
    # filtering data from occurrence data frame
    df_occ <- 
      df.occ[, c(species, Max.age, Min.age, site)]
    vars <- list(species, Max.age, Min.age, site)
    name_vars <- c("species", "Max.age", "Min.age", "site")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names
    
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
    
    list_site_spp_interval <- 
      lapply(1:length(species_slices), function(x){
        names_slices <- as.numeric(names(species_slices)[[x]])
        names_species <- species_slices[[x]]
        filtered_df <- 
          df_occ %>% filter(species %in% names_species)
        
        filtered_df_site <-  
          filtered_df |> 
          filter(Max.age >= names_slices & Min.age <= names_slices)
        
        filtered_df_site2 <- 
          filtered_df_site |> 
          group_by(site) |> 
          distinct(species, .keep_all = T) |> 
          add_count(site, name = "species.per.site") |> 
          mutate(time.slice = as.character(names_slices)) |> 
          select(species, site, species.per.site, time.slice)
        return(filtered_df_site2)
      })
    
    
    df_long_slice_site_spp <- do.call(rbind, list_site_spp_interval)
    
    # calculating distance matrix for all species
    all_dist <- dist(x = df.TS.TE$trait, method = "euclidean", diag = T, upper = T) 
    matrix_all_dist <- as.matrix(all_dist)
    colnames(matrix_all_dist) <- df.TS.TE$species
    rownames(matrix_all_dist) <- df.TS.TE$species
    
    # computing a list with mean mpd per species per site per time slices
    
    slices_names <- unique(df_long_slice_site_spp$time.slice)
    site_names <- unique(df_long_slice_site_spp$site)
    list_slice_mean_mpd_site <- vector(mode = "list", length = length(slices_names))
    for(i in 1:length(slices_names)){
      # i = 33
      df_slice <- df_long_slice_site_spp[which(df_long_slice_site_spp$time.slice == slices_names[i]), ]
      sites_vec <- unique(df_slice$site)
      
      for(j in 1:length(sites_vec)){
        # j = 1
        species_slice_site <- df_slice[which(df_slice$site == sites_vec[j]), "species"]
        matrix_dist <- matrix_all_dist[rownames(matrix_all_dist) %in% species_slice_site$species, 
                                       colnames(matrix_all_dist) %in% species_slice_site$species]
        if(length(dim(matrix_dist)) < 2){
          mpd_res <- NA
        } else{
          mpd_res <- colMeans(matrix_dist)
        }
        list_sites[[j]] <- mpd_res
        
      }
      mpd_vector_slice_site <- unlist(list_sites)
      if(all(is.na(mpd_vector_slice_site)) == TRUE){
        list_slice_mean_mpd_site[[i]] <- NA
      } else{
        df_mpd_slice_site <- data.frame(species = names(mpd_vector_slice_site)[!is.na(mpd_vector_slice_site)], 
                                        mpd = as.numeric(mpd_vector_slice_site)[!is.na(mpd_vector_slice_site)])
        df_mpd_slice_site2 <- 
          df_mpd_slice_site |>
          group_by(species) |> 
          mutate(mean.mpd = mean(mpd), var.mpd = var(mpd)) |> 
          distinct(species, .keep_all = T) |> 
          select(species, mean.mpd, var.mpd)
        list_slice_mean_mpd_site[[i]] <- df_mpd_slice_site2
      }
      
    }
    
    names(list_slice_mean_mpd_site) <- slices_names
    
    # creating a column to identify the time slice for each species mean mpd
    
    list_df_mpd_site_slice <- 
      lapply(names(list_slice_mean_mpd_site), function(name) {
        if(all(is.na(list_slice_mean_mpd_site[[name]])) == TRUE){
          df <- NA
        } else{
          df <- list_slice_mean_mpd_site[[name]]
          df$time.slice <- name  # Add the name of the element as the new column
        }
        return(df)
      })
    
    # removing slices with NAs
    remove_no_occ <- which(is.na(list_df_mpd_site_slice) == TRUE)
    list_df_mpd_site_slice2 <- list_df_mpd_site_slice[-remove_no_occ]
    
    df_mpd_site_slice <- do.call(rbind, list_df_mpd_site_slice2)
    
    return(df_mpd_site_slice)
    
  }