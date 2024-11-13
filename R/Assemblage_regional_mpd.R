#' Compute mean pairwise distance for each site in time slices
#'
#' @param df.TS.TE
#' @param df.occ
#' @param time.slice
#' @param grid.size
#' @param trait
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
Assemblage_regional_mpd <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           grid.size,
           trait,
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
      df.TS.TE[, c(species, TS, TE, trait)]
    vars <- list(species, TS, TE, trait)
    name_vars <- c("species", "TS", "TE", "trait")
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

    # trait data
    matrix_dist_trait <-
      as.matrix(dist(x = df.TS.TE[, "trait"], method = "euclidean", upper = T, diag = T))

    rownames(matrix_dist_trait) <- df.TS.TE$species
    colnames(matrix_dist_trait) <- df.TS.TE$species

    # computing ses mpd

    res_mpd <-
      lapply(1:length(list_interval), function(x){
        if(!is.data.frame(list_interval[[x]]) == TRUE){
          NA
        } else{
          # site names
          sites_test <- names(table(list_interval[[x]]$grid_id.x))
          geometry_site <- list_interval[[x]][match(sites_test, list_interval[[x]]$grid_id.x), "geometry"]

          # building community matrix
          comm_mat_test <-
            lapply(sites_test, function(y){
              spp <- as.data.frame(list_interval[[x]][which(y == list_interval[[x]]$grid_id.x), "species"])$species
              comm_mat <- matrix(1, nrow = 1, ncol = length(spp), dimnames = list("comm", spp))
              comm_mat
            })

          # calculating mpd
          list_sesmpd_test <-
            lapply(comm_mat_test, function(x){
              ses.mpd.modif(samp = x, dis = matrix_dist_trait, runs = 1, null.model = "taxa.labels")
            })

          # naming list
          names(list_sesmpd_test) <- sites_test

          # response data frame
          data.frame(do.call(rbind, list_sesmpd_test), site = sites_test, time.slice = names(list_interval)[x], geometry = geometry_site)

        }
      })

    names(res_mpd) <- seq_interval

    res <- do.call(rbind, res_mpd)

    res_mean_mpd_timeslice <-
      res |>
      group_by(time.slice) |>
      mutate(mean.mpd = mean(mpd.obs, na.rm = TRUE)) |>
      mutate(time.slice = as.numeric(time.slice))

    # calculating the mean mpd and variance for mpd at grid level


    df_grid_mean <-
      res_mean_mpd_timeslice |>
      ungroup() |>
      group_by(site) |>
      mutate(grid.mean = mean(mpd.obs, na.rm = TRUE)) |>
      rename(grid_id = site) |>
      select(grid_id, grid.mean, geometry) |>
      distinct(grid_id, .keep_all = TRUE)

    # joining with grid information

    #sf_grid_mean_mpd <-
    #  sf_grid2 |>
    #  left_join(df_grid_mean, by = c(grid_id = "grid_id")) |>
    #  distinct(grid_id, .keep_all = TRUE)


    # returning both mpd by time slice and mean value per grid

    list_res <- vector(mode = "list")

    list_res$mean_mpd_grid <- df_grid_mean
    list_res$mean_mpd_timeslice <- res_mean_mpd_timeslice

    return(list_res)


  }
