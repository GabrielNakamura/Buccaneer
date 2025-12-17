#' Compute Individual Species Coexistence based on reach criteria
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
#'
#' @returns A data frame with three columns:
#'   \item{time.slice}{Numeric. The time point representing each slice, typically
#'       the upper (older) boundary of the time bin.}
#'   \item{species}{Character. The name of all species with at least one
#'       cooccurrence with another species.}
#'   \item{n.coexistence}{Numeric. The number of co-occurring species in each
#'       time slice per species.}
#
#' @export
#'
#' @examples
IndivSpecies_reach_coexistence <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           round.digits,
           species,
           TS,
           TE,
           lat,
           lon,
           Max.age = "Max.age",
           Min.age = "Min.age",
           crs = 4326,
           remove.singletons = TRUE,
           group = NULL,
           group.focal.compare = NULL,
           type.comparison = NULL){

    # renaming longevity data frame
    if(!is.null(group) == TRUE){
      df.TS.TE <- df.TS.TE[, c(species, TS, TE, group)]
      colnames(df.TS.TE) <- c("species", "TS", "TE", "group")
    } else{
      df.TS.TE <- df.TS.TE[, c(species, TS, TE)]
      colnames(df.TS.TE) <- c("species", "TS", "TE")
    }

    # renaming occurrence data frame
    df_occ <-
      df.occ[, c(species, lat, lon, Max.age, Min.age)]
    vars <- list(species, lat, lon, Max.age, Min.age)
    name_vars <- c("species", "lat", "lon", "Max.age", "Min.age")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names

    # time slices
    names_slice <- seq(from = ceiling(max(df.TS.TE[, "TS"])),
                       to = ceiling(min(df.TS.TE[, "TE"])),
                       by = -time.slice)

    # Removing non numeric values in lat long in occurrence data frame
    df_long_slice_coord_spp2 <-
      df_occ |>
      filter(if_any(c(lat, lon), ~ !is.na(as.numeric(.)))) |>
      mutate(lat = as.numeric(lat), lon = as.numeric(lon))

    # transforming to a coordinate system  in occurrence data frame
    df_long_slice_coord_spp3 <-
      sf::st_as_sf(x = df_long_slice_coord_spp2,
                   coords = c("lon", "lat"),
                   crs = crs)


    # Time coexistence matrix for all species
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = round.digits,
                               species = "species",
                               TS = "TS",
                               TE = "TE")

    # modified coexistence matrix containing group comparison, in rows are focal species, in colums comparison
    # the "within" argument place focal species in both columns and rows
    if(!is.null(group.focal.compare) == TRUE){
      focal <- group.focal.compare[1]
      compare <- group.focal.compare[2]
      spp_focal <- df.TS.TE[which(df.TS.TE$group == focal), "species"]$species
      spp_compare <- df.TS.TE[which(df.TS.TE$group == compare), "species"]$species

      if(type.comparison == "between"){# comparison between species of two groups
        matrix_coex <- lapply(matrix_coex, function(x) x[spp_focal, spp_compare]) # focal species in lines and comparison in columns
      }
      if(type.comparison == "within"){ # comparison only within the focal group
        matrix_coex <- lapply(matrix_coex, function(x) x[spp_focal, spp_focal])
      }
    } else{
      matrix_coex <- matrix_coex
    }

    # species composition at each time slice
    spp_slice <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    # naming list with time slices
    names(spp_slice) <- format(names_slice, trim = TRUE, scientific = FALSE)


    # obtaining occurrence records for each time slice based on species composition in time slices
    df_occ_occurr <-
      comp_slice_occ_cooccurr(spp_slice = spp_slice,
                              df.occ = df_long_slice_coord_spp3,
                              species = "species",
                              Max.age = "Max.age",
                              Min.age = "Min.age")

    # list with matrices containing geographical distances for all occurrence records in all time slices
    list_matrix_dist_occ <-
      lapply(df_occ_occurr, function(x){
        dist_occ <- sf::st_distance(x = x)
        colnames(dist_occ) <- x$species
        rownames(dist_occ) <- x$species
        return(dist_occ)
      })

    # getting pairwise combination between all species cooccurring in all time slices
    combination <-
      lapply(df_occ_occurr, function(x){
        spp_names <- unique(x$species)
        if(length(spp_names) == 1 | length(spp_names) == 0){
          df <- NA
        } else{
          df <- data.frame(combn(spp_names, m = 2, simplify = TRUE))
        }
        return(df)
      })

    # removing NAs - no co-occurrence, only one species occurrence in the timeslice
    combination2 <- combination[-which(is.na(combination) == TRUE)]
    matrix_all_dist2 <- list_matrix_dist_occ[-which(is.na(combination) == TRUE)]
    matrix_coex2 <- matrix_coex[-which(is.na(combination) == TRUE)]
    # filtering the combinations and filling the co-occurrence matrix with zeroes and 1 accordingly with reach criteria
    list_res <- vector(mode = "list", length = length(matrix_coex2))
    for(i in 1:length(matrix_all_dist2)){
      # i = 100
      #res <-
      #  matrix(0, nrow = length(unique(rownames(matrix_all_dist2[[i]]))),
      #         ncol = length(unique(colnames(matrix_all_dist2[[i]]))),
      #         dimnames = list(unique(rownames(matrix_all_dist2[[i]])),
      #                         unique(colnames(matrix_all_dist2[[i]]))))
      res <- matrix_coex2[[i]] # regional cooccurrence matrix
      for(j in 1:ncol(combination2[[i]])){
        # j = 5

        min_between <- matrix_all_dist2[[i]][which(combination2[[i]][1, j] == rownames(matrix_all_dist2[[i]])),
                                             which(combination2[[i]][2, j] == colnames(matrix_all_dist2[[i]])),
                                             drop = FALSE]
        max_1 <- matrix_all_dist2[[i]][which(combination2[[i]][1, j] == rownames(matrix_all_dist2[[i]])),
                                       which(combination2[[i]][1, j] == colnames(matrix_all_dist2[[i]])),
                                       drop = FALSE]
        max_2 <- matrix_all_dist2[[i]][which(combination2[[i]][2, j] == rownames(matrix_all_dist2[[i]])),
                                       which(combination2[[i]][2, j] == colnames(matrix_all_dist2[[i]])),
                                       drop = FALSE]
        res[combination2[[i]][1, j], combination2[[i]][2, j]] <- ifelse(min(min_between) <= sum(max(max_1), max(max_2)), 1, 0)

      }
      list_res[[i]] <- res
    }

    # naming coexistence matrices
    names_slice2 <- names_slice[-which(is.na(combination) == TRUE)]
    names(list_res) <- format(names_slice2, trim = TRUE, scientific = FALSE)

    # calculating coexistence
    res_coex_slice <-
      lapply(list_res, function(x){
        coex_slice <- rowSums(x) - 1
        coex_slice2 <- coex_slice[which(coex_slice > 0)]
        return(coex_slice2)
      })

    # organizing in a data frame
    df_long_res <-
      do.call(rbind, lapply(names(res_coex_slice), function(nm) {
        data.frame(
          time.slice = nm,
          species = names(res_coex_slice[[nm]]),
          n.coexistence = unname(res_coex_slice[[nm]])
        )
      })
      )

    # removing or not singleton species
    if(remove.singletons != TRUE){
      df_long_res2 <- df_long_res
    } else{
      df_long_res2  <-
        df_long_res |>
        filter(n.coexistence != 0)
    }

    # output
    return(df_long_res2)

  }
