#' Compute individual species coexistence based on reach criteria
#'
#' @param df.TS.TE a data frame object containing at least three columns. Species names,
#'     origination time and extinction time for each species.
#' @param df.occ a data frame object containing the occurrence records for each species.
#'     This must have at least a column indicating the name of species, its minimum and maximum age estimate,
#'     and its site location ID.
#' @param time.slice scalar indicating the time interval between consecutive time slices.
#' @param round.digits scalar indicating the number of digits for time of origination and time for
#'     extinction.
#' @param species character indicating the name of the column of the data frame
#'     containing the species name information.
#' @param TS character indicating the name of the columns of the data frame
#'     containing the information on origination time.
#' @param TE character indicating the name of the column of the data frame
#'     containing the information on extinction time.
#' @param Max.age character indicating the name of the column containing the upper age limit for occurrence record.
#' @param Min.age character indicating the name of the column containing the lower age limit for occurrence record.
#' @param lat character indicating the name of the column containing the latitude of occurrence record in df.occ.
#' @param lon character indicating the name of the column containing the longitude of occurrence record in df.occ.
#'
#' @return data frame containing the name of species, its mean coexistence value calculated by each time slice
#'     considering reach criteria of co-occurrence
#'
#' @export
#'
#' @examples
IndivSpec_reach_richness <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE",
           Max.age = "Max.age",
           Min.age = "Min.age",
           lat = "lat",
           lon = "lng"){

    # filtering data from longevities data frame
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
          mutate(time.slice = as.character(names_slices)) |>
          select(species, time.slice, lat, lon)
        return(filtered_df_site2)
      })

    df_long_slice_coord_spp <- do.call(rbind, list_site_spp_interval)

    # reach matrices - at this point use auxiliar function for reach calculation

    list_reach <- calc_reach(df.coords = df_long_slice_coord_spp)

    # computing coexistences

    n_coexistence_species <-
      lapply(list_reach, function(x){
        colSums(x)
      })

    df_coexistence_reach <-
      lapply(n_coexistence_species, function(x){
        data.frame(species = names(x), n.coexistence = as.numeric(x))
      })

    # transforming for long format and adding to the df the time slice column
    df_coexistence_reach2 <-
      lapply(names(df_coexistence_reach), function(name) {
        if(nrow(df_coexistence_reach[[name]]) == 0){
          df <- NA
        } else{
          df <- df_coexistence_reach[[name]]
          df$time.slice <- name  # Add the name of the element as the new column
        }
        return(df)
      })


    # results with individual species coexistence
    df_coexistence_reach3 <- do.call(rbind, df_coexistence_reach2)
    return(df_coexistence_reach3)

  }
