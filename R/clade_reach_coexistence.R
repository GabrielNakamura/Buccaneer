#' Compute species coexistence for each time slice based on reach metric
#'
#' @param df.TS.TE A data frame containing species temporal data with at least
#'     three columns: species names, origination times (TS), and extinction
#'     times (TE). Additional columns may include group assignments.
#' @param df.occ A data frame containing fossil occurrence records with at least
#'     four columns: species names, minimum age, maximum age, and site location ID.
#'     Each row represents a single occurrence record at a specific site.
#' @param time.slice Numeric. The time interval (in the same units as TS and TE)
#'     between consecutive time slices for temporal binning.
#' @param round.digits Integer. The number of decimal places to round time slice
#'     values. Default is 1. This affects temporal binning precision.
#' @param species Character. The name of the column in \code{df.TS.TE} and
#'     \code{df.occ} containing species identifiers. Default is "species".
#' @param TS Character. The name of the column in \code{df.TS.TE} containing
#'     origination (first appearance) times for each species. Default is "TS".
#' @param TE Character. The name of the column in \code{df.TS.TE} containing
#'     extinction (last appearance) times for each species. Default is "TE".
#' @param lat Numeric. The latitude coordinate of the occurrence record in
#'     \code{df.occ}.
#' @param lon Numeric. The longitude coordinate of the occurrence records in
#'     \code{df.occ}.
#' @param Max.age Character. The name of the column in \code{df.occ} containing
#'     the maximum (oldest) age estimate for each occurrence record. Default is "Max.age".
#' @param Min.age Character. The name of the column in \code{df.occ} containing
#'     the minimum (youngest) age estimate for each occurrence record. Default is "Min.age".
#' @param crs Numeric. The code indicating the coordinate reference system
#'     to be used for latitude and longitude of occurrence records in \code{df.occ}
#' @param remove.singletons Logical. Should singleton species (species occurring
#'     alone at a site with no co-occurring species) be excluded from mean and
#'     variance calculations? Default is TRUE. When TRUE, singletons are treated
#'     as NA; when FALSE, they contribute 0 to the mean.
#' @param group Character. The name of the column in \code{df.TS.TE} containing
#'     group assignments for species (e.g., clade, family). Required if using
#'     \code{group.focal.compare}. Default is NULL.
#' @param group.focal.compare Character vector of length 2. The first element
#'     specifies the focal group and the second specifies the comparison group.
#'     If NULL (default), coexistence is calculated across all species regardless
#'     of group membership.
#' @param type.comparison Character. Specifies the type of coexistence comparison:
#'     \itemize{
#'       \item \code{"between"}: Count only co-occurrences between species from
#'             the focal and comparison groups.
#'       \item \code{"within"}: Count only co-occurrences among species within
#'             the focal group.
#'       \item NULL (default): Count all co-occurrences regardless of group.
#'     }
#'
#' @returns A data frame with three columns:
#'   \item{mean.coexistence}{Numeric. The mean number of co-occurring species in
#'       each time slice. This represents average local coexistence
#'       (excluding or including singletons based on \code{remove.singletons}).}
#'   \item{var.coexistence}{Numeric. The variance in the number of co-occurring
#'       species in each time slice.}
#'   \item{time.slice}{Numeric. The time point representing each slice, typically
#'       the upper (older) boundary of the time bin.}
#' @examples
#' \dontrun{
#' # Create example fossil data
#' df_temporal <- data.frame(
#'   species = c("sp1", "sp2", "sp3", "sp4"),
#'   TS = c(100, 95, 90, 85),
#'   TE = c(50, 45, 40, 35),
#'   group = c("A", "A", "B", "B")
#' )
#'
#' df_occurrences <- data.frame(
#'   species = c("sp1", "sp1", "sp2", "sp3", "sp4", "sp4"),
#'   Max.age = c(100, 95, 95, 90, 85, 85),
#'   Min.age = c(90, 85, 85, 80, 75, 75),
#'   site = c("site1", "site2", "site1", "site1", "site2", "site3")
#' )
#'
#' # Calculate mean site reach through time
#' result <- clade_reach_coexistence(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10
#' )
#'
#' # View results
#' head(result)
#'
#' # Calculate reach between groups
#' result_between <- clade_reach_coexistence(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10,
#'   group = "group",
#'   group.focal.compare = c("A", "B"),
#'   type.comparison = "between"
#' )
#'
#' # Calculate reach within a single group
#' result_within <- clade_reach_coexistence(
#'   df.TS.TE = df_temporal,
#'   df.occ = df_occurrences,
#'   time.slice = 10,
#'   group = "group",
#'   group.focal.compare = c("A", "B"),
#'   type.comparison = "within"
#' )
#' }
#' @export

clade_reach_coexistence <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           round.digits,
           species = "species",
           TS = "TS",
           TE = "TE",
           lat = "lat",
           lon = "lon",
           Max.age = "Max.age",
           Min.age = "Min.age",
           crs = 4326,
           remove.singletons = TRUE,
           group = NULL,
           group.focal.compare = NULL,
           type.comparison = NULL){

    if(!is.null(group)){
      df.TS.TE <- df.TS.TE[, c(species, TS, TE, group)]
      colnames(df.TS.TE) <- c("species", "TS", "TE", "group")
    } else{
      df.TS.TE <- df.TS.TE[, c(species, TS, TE)]
      colnames(df.TS.TE) <- c("species", "TS", "TE")
    }

    df_occ <- df.occ[, c(species, lat, lon, Max.age, Min.age)]
    vars <- list(species, lat, lon, Max.age, Min.age)
    name_vars <- c("species", "lat", "lon", "Max.age", "Min.age")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names

    names_slice <- seq(from = ceiling(max(df.TS.TE[, "TS"])),
                       to   = ceiling(min(df.TS.TE[, "TE"])),
                       by   = -time.slice)

    df_long_slice_coord_spp2 <-
      df_occ |>
      dplyr::filter(dplyr::if_any(c(lat, lon), ~ !is.na(as.numeric(.)))) |>
      dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon))

    df_long_slice_coord_spp3 <-
      sf::st_as_sf(x = df_long_slice_coord_spp2,
                   coords = c("lon", "lat"),
                   crs = crs)

    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = round.digits,
                               species = "species", TS = "TS", TE = "TE")

    if(!is.null(group.focal.compare)){
      focal    <- group.focal.compare[1]
      compare  <- group.focal.compare[2]
      spp_focal    <- df.TS.TE[which(df.TS.TE$group == focal), "species"]
      spp_compare  <- df.TS.TE[which(df.TS.TE$group == compare), "species"]

      if(type.comparison == "between"){
        matrix_coex <- lapply(matrix_coex, function(x) x[spp_focal, spp_compare])
      }
      if(type.comparison == "within"){
        matrix_coex <- lapply(matrix_coex, function(x) x[spp_focal, spp_focal])
      }
    }

    spp_slice <- lapply(matrix_coex, function(x) names(which(rowSums(x) >= 1)))
    names(spp_slice) <- format(names_slice, trim = TRUE, scientific = FALSE)

    df_occ_occurr <-
      comp_slice_occ_cooccurr(spp_slice = spp_slice,
                              df.occ = df_long_slice_coord_spp3,
                              species = "species",
                              Max.age = "Max.age",
                              Min.age = "Min.age")

    list_matrix_dist_occ <-
      lapply(df_occ_occurr, function(x){
        dist_occ <- sf::st_distance(x = x)
        colnames(dist_occ) <- x$species
        rownames(dist_occ) <- x$species
        dist_occ
      })

    combination <-
      lapply(df_occ_occurr, function(x){
        spp_names <- unique(x$species)
        if(length(spp_names) <= 1){
          NA
        } else{
          data.frame(utils::combn(spp_names, m = 2, simplify = TRUE))
        }
      })

    combination2     <- combination[!is.na(combination)]
    matrix_all_dist2 <- list_matrix_dist_occ[!is.na(combination)]
    matrix_coex2     <- matrix_coex[!is.na(combination)]

    list_res <- vector(mode = "list", length = length(matrix_coex2))
    for(i in seq_along(matrix_all_dist2)){
      res <- matrix_coex2[[i]]
      for(j in 1:ncol(combination2[[i]])){
        min_between <- matrix_all_dist2[[i]][combination2[[i]][1, j],
                                             combination2[[i]][2, j], drop = FALSE]
        max_1 <- matrix_all_dist2[[i]][combination2[[i]][1, j],
                                       combination2[[i]][1, j], drop = FALSE]
        max_2 <- matrix_all_dist2[[i]][combination2[[i]][2, j],
                                       combination2[[i]][2, j], drop = FALSE]
        res[combination2[[i]][1, j], combination2[[i]][2, j]] <-
          ifelse(min(min_between) <= sum(max(max_1), max(max_2)), 1, 0)
      }
      list_res[[i]] <- res
    }

    names_slice2 <- names_slice[!is.na(combination)]
    names(list_res) <- format(names_slice2, trim = TRUE, scientific = FALSE)

    res_coex_slice <-
      lapply(list_res, function(x){
        coex_slice <- rowSums(x) - 1
        coex_slice[coex_slice > 0]
      })

    res_mean_coex_slice <- lapply(res_coex_slice, mean)
    res_var_coex_slice  <- lapply(res_coex_slice, var)

    df_res <-
      data.frame(mean.coexistence = unlist(res_mean_coex_slice),
                 var.distance     = unlist(res_var_coex_slice),
                 time.slice       = names_slice2)

    return(df_res)
  }
