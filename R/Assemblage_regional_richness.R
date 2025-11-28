#' Calculate Species Richness Across Time Slices and Grid Cells
#'
#' This function computes species richness for fossil assemblages across spatial
#' grid cells and temporal bins. It divides geographic space into a regular grid
#' and temporal ranges into discrete time slices, then calculates the number of
#' species present in each grid cell during each time slice based on occurrence
#' records and species longevities.
#'
#' @param df.TS.TE A data frame containing species temporal data with at least
#'     three columns: species names, origination times (TS), and extinction
#'     times (TE).
#' @param df.occ A data frame containing fossil occurrence records with at least
#'     five columns: species names, minimum age, maximum age, latitude, and
#'     longitude. Each row represents a single occurrence record at a specific
#'     location.
#' @param time.slice Numeric. The time interval (in the same units as TS and TE)
#'     between consecutive time slices for temporal binning.
#' @param grid.size Numeric. The size of grid cells in degrees for spatial binning.
#'     For example, \code{grid.size = 5} creates 5° × 5° grid cells.
#' @param round.digits Integer. The number of decimal places to round time slice
#'     values. Default is 1. This affects temporal binning precision.
#' @param species Character. The name of the column in \code{df.TS.TE} and
#'     \code{df.occ} containing species identifiers. Default is "species".
#' @param TS Character. The name of the column in \code{df.TS.TE} containing
#'     origination (start) times for each species. Default is "TS".
#' @param TE Character. The name of the column in \code{df.TS.TE} containing
#'     extinction (end) times for each species. Default is "TE".
#' @param Max.age Character. The name of the column in \code{df.occ} containing
#'     the maximum (oldest) age estimate for each occurrence record. Default is "Max.age".
#' @param Min.age Character. The name of the column in \code{df.occ} containing
#'     the minimum (youngest) age estimate for each occurrence record. Default is "Min.age".
#' @param lat Character. The name of the column in \code{df.occ} containing
#'     latitude coordinates for occurrence records. Default is "lat".
#' @param lon Character. The name of the column in \code{df.occ} containing
#'     longitude coordinates for occurrence records. Default is "lng".
#' @param crs Numeric or character. The coordinate reference system (CRS) code
#'     for spatial data. Default is 4326 (WGS84 geographic coordinates).
#'
#' @return A list containing two elements:
#'   \item{grid_mean_age}{An sf object with grid cell geometries and temporal
#'       summaries. Contains columns:
#'       \itemize{
#'         \item \code{grid_id}: Unique identifier for each grid cell
#'         \item \code{mean.age.grid}: Mean age of all occurrences in the grid cell
#'         \item \code{var.age.grid}: Variance of ages in the grid cell
#'         \item \code{geometry}: Spatial geometry of the grid cell
#'       }}
#'   \item{time_series_rich}{A data frame with richness metric for each grid
#'       cell and time slice. Contains columns:
#'       \itemize{
#'         \item \code{grid_id}: Grid cell identifier
#'         \item \code{mean.age.grid}: Mean age of occurrences in the grid cell
#'         \item \code{var.age.grid}: Variance of ages in the grid cell
#'         \item \code{time.slice}: Time slice identifier
#'         \item \code{rich.by.grid}: Species richness in the grid cell during the time slice
#'         \item \code{mean.rich.by.grid.slice}: Mean richness across all grid cells in the time slice
#'         \item \code{geometry}: Spatial geometry of the grid cell
#'       }}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates a spatial grid covering the extent of occurrence records
#'   \item Assigns each occurrence to a grid cell based on coordinates
#'   \item Calculates temporal statistics (mean and variance of ages) for each grid cell
#'   \item Divides the temporal range into discrete time slices
#'   \item For each time slice, identifies species present based on longevities
#'   \item Counts species richness in each grid cell for each time slice
#'   \item Computes mean richness across grid cells within each time slice
#' }
#'
#' Missing or invalid coordinates (NA values) are automatically removed before
#' spatial analysis. Grid cells without occurrences in a time slice will not
#' appear in the results for that time slice.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example data
#' df_temporal <- data.frame(
#'   species = c("sp1", "sp2", "sp3", "sp4"),
#'   TS = c(100, 95, 90, 85),
#'   TE = c(50, 45, 40, 35)
#' )
#'
#' df_occurrences <- data.frame(
#'   species = c("sp1", "sp1", "sp2", "sp3", "sp4"),
#'   Max.age = c(100, 95, 95, 90, 85),
#'   Min.age = c(90, 85, 85, 80, 75),
#'   lat = c(10.5, 15.2, 20.1, 25.3, 30.0),
#'   lng = c(-50.0, -55.5, -60.2, -65.0, -70.5)
#' )
#'
#' # Calculate regional richness with 10 Ma time slices and 5° grid cells
#' results <- Assemblage_regional_richness(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10,
#'   grid.size = 5
#' )
#'
#' # Access grid-level summaries
#' grid_summary <- results$grid_mean_age
#'
#' # Access time series of richness
#' richness_timeseries <- results$time_series_rich
#'
#' # Plot richness for a specific time slice
#' library(ggplot2)
#' slice_data <- subset(results$time_series_rich, time.slice == 90)
#' ggplot(slice_data) +
#'   geom_sf(aes(fill = rich.by.grid)) +
#'   scale_fill_viridis_c() +
#'   theme_minimal()
#' }
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

    # grid information with age
    sf_grid_mean_age <-
      sf::st_join(sf_grid, df_mean_var_age_grid)
    sf_grid_mean_age2 <-
      sf_grid_mean_age |>
      select(grid_id.x, mean.age.grid, var.age.grid, geometry) |>
      rename(grid_id = "grid_id.x")


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


    list_df_res_rich2 <- list_res_rich[!is.na(list_res_rich)]

    df_long_rich_assemblage <- do.call(rbind, list_df_res_rich2)
    df_long_rich_assemblage2 <-
      df_long_rich_assemblage |>
      select(grid_id.x, mean.age.grid, var.age.grid, time.slice, rich.by.grid, mean.rich.by.grid.slice, geometry) |>
      rename(grid_id = grid_id.x)


    # list with results

    list_res <- vector(mode = "list")

    # joining grid with mean age

    list_res$grid_mean_age <-
      sf_grid_mean_age2 |>
      mutate(grid_id = as.character(grid_id))# metrics at grid level
    list_res$time_series_rich <- df_long_rich_assemblage2 |>
      mutate(grid_id = as.character(grid_id))
    # metrics at time slice level
    return(list_res)



  }



