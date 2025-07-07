#' Compute mean distance for individual species cooccurring in sites
#'
#' @param df.TS.TE Data frame object containing at least four columns. Species names,
#'     origination time, extinction time and a trait value for each species.
#' @param df.occ a data frame object containing the occurrence records for each species.
#'     This must have at least a column indicating the name of species, its minimum and maximum age estimate,
#'     and its site location ID.
#' @param time.slice Scalar indicating the time interval between consecutive time slices.
#' @param dist.trait A dist object containing the clade pairwise distance matrices.
#'     The name of the clades in this object must be equal to the name of the
#'     clades in df.TS.TE data frame.
#' @param nearest.taxon A scalar indicating the number of nearest species/genus that will be used.
#'     1 computes mnnd metric and the option "all" computes mpd.
#' @param group Character indicating the name of the column that contain the groups that will be used
#'     in comparison.
#' @param group.focal.compare Character vector indicating the focal (first element) and comparison (second element)
#'     groups used in the calculation. If NULL, the default, the metrics  will be calculated
#'     using all  clades.
#' @param type.comparison Character. It can be "between" to compute distances only between species/genus of two groups
#'     or "within" to calculate distance only inside the focal group. If null the distance is computed
#'     considering all clades together
#' @param trait Character indicating the name of the column containing values of the traits for each species. If NULL,
#'     the default, the user must provide a distance matrix.
#' @param round.digits Scalar indicating the precision of time slices.
#' @param species Character indicating the name of the column of the data frame
#'     containing the species name information.
#' @param TS haracter indicating the name of the columns of the data frame
#'     containing the information on origination time.
#' @param TE Character indicating the name of the column of the data frame
#'     containing the information on extinction time.
#' @param Max.age Character indicating the name of the column containing the upper age limit for occurrence record.
#' @param Min.age Character indicating the name of the column containing the lower age limit for occurrence record.
#' @param site Character indicating the name of the column containing the information on site location.
#'
#' @returns
#' @export
#'
#' @examples
IndivSpec_site_distance <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           dist.trait,
           nearest.taxon,
           group = NULL,
           group.focal.compare = NULL,
           type.comparison = NULL,
           trait = NULL,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE",
           Max.age = "Max.age",
           Min.age = "Min.age",
           site = "site"){
    # subseting columns
    if(!is.null(group) == TRUE){
      df.TS.TE <- df.TS.TE[, c(species, trait, TS, TE, group)]
      if(is.null(trait) == TRUE){
        colnames(df.TS.TE) <- c("species", "TS", "TE", "group")
      } else{
        colnames(df.TS.TE) <- c("species", "trait", "TS", "TE", "group")
      }
    } else{
      df.TS.TE <- df.TS.TE[, c(species, trait, TS, TE)]
      if(is.null(trait) == TRUE){
        colnames(df.TS.TE) <- c("species", "TS", "TE")
      } else{
        colnames(df.TS.TE) <- c("species", "trait", "TS", "TE")
      }
    }

    df_occ <-
      df.occ[, c(species, Max.age, Min.age, site)]
    vars <- list(species, Max.age, Min.age, site)
    name_vars <- c("species", "Max.age", "Min.age", "site")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names

    # Generating time intervals used to compute temporal coexistence
    seq_interval <- seq(from = max(df.TS.TE[, "TS"]), to = min(df.TS.TE[, "TE"]), by = -time.slice)
    seq_interval <- c(round(seq_interval, digits = round.digits))

    # coexistence matrix
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = 1,
                               species = "species",
                               TS = "TS",
                               TE = "TE")

    # species composition at each timeslice
    spp_slice <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    names(spp_slice) <- seq_interval

    # calculating trait distances for all clades

    # trait distance
    if(!is.null(dist.trait) == TRUE){
      matrix_dist_trait <- as.matrix(dist.trait)
    } else{
      matrix_dist_trait <-
        as.matrix(dist(x = df.TS.TE[, "trait"], method = "euclidean", upper = T, diag = T))
    }
    rownames(matrix_dist_trait) <- df.TS.TE$species
    colnames(matrix_dist_trait) <- df.TS.TE$species


    # modified trait matrix containing group comparison
    if(!is.null(group.focal.compare) == TRUE){
      focal <- group.focal.compare[1]
      compare <- group.focal.compare[2]
      spp_focal <- df.TS.TE[which(df.TS.TE$group == focal), "species"]$species
      spp_compare <- df.TS.TE[which(df.TS.TE$group == compare), "species"]$species

      if(type.comparison == "between"){# comparison between species of two groups
        matrix_dist_trait_comp <- matrix_dist_trait[spp_focal, spp_compare] # focal speices in lines and comparison in columns
      }
      if(type.comparison == "within"){ # comparison only within the focal group
        matrix_dist_trait_comp <- matrix_dist_trait[spp_focal, spp_focal]
      }
    }


    # calculating matrix of species cooccurrence for site

    list_matrix_cooccur_site <-
      comp_site_cooccurr(spp_slice = spp_slice, df.occ = df_occ)


    # Ensure same species and same order in both cooccurrence matrix and distance matrix
    list_dist_spp <-
      lapply(list_matrix_cooccur_site, function(x){
        species_row <- intersect(rownames(x), rownames(matrix_dist_trait_comp))
        species_col <- intersect(colnames(x), colnames(matrix_dist_trait_comp))
        cooccur_matrix <- x[species_row, species_col]
        dist_matrix <- matrix_dist_trait_comp[species_row, species_col]

        # Create matrix to store the results of mean pairwise distances
        mean_distances <- matrix(NA, nrow = length(species_row), ncol = 1,
                                 dimnames = list(species_row, paste("mean.dist.to.cooccur", nearest.taxon, sep = ".")))

        # calculating distances for all species
        for (sp in species_row) {

          #checking if there are no species in the slice or
          if(length(species_row) == 0 | length(species_col) == 0){
            mean_distances[sp, 1] <- NA
            if(length(species_row) != 0){
              mean_distances[sp, 1] <- "NA_singleton"
            }
          } else{
            cooccur_species <- names(which(cooccur_matrix[sp, ] > 0 & names(cooccur_matrix[sp, ]) != sp))

            # If there are co-occurring species, compute mean distance
            if (length(cooccur_species) > 0) {
              distance_sorted <- sort(dist_matrix[sp, cooccur_species], decreasing = FALSE)
              if(nearest.taxon == "all"){ # calculating for all taxon
                mean_distances[sp, 1] <- mean(as.numeric(dist_matrix[sp, cooccur_species]), na.rm = TRUE)
              } else{ # using the threshold distance set by the user
                mean_distances[sp, 1] <- mean(distance_sorted[1:nearest.taxon], na.rm = TRUE)
              }
            }
          }
        }
        return(mean_distances)
      })

    mean_dist_timeslice <- lapply(list_dist_spp, function(x) mean(x, na.rm = TRUE))
    var_dist_timeslice <- lapply(list_dist_spp, function(x) var(x, na.rm = TRUE))

    lapply(list_dist_spp, function(x){
      if(!is.matrix(x) == TRUE){
        NA
      } else{
        data.frame(species = rownames(x), mean.dist = colnames(x))
      }
    })
    names(list_dist_spp) <- seq_interval


    # Suppose your list is called `mean_dist_list`
    # Each element is a matrix or data frame

    df_dist_spp <-
      do.call(rbind, lapply(names(list_dist_spp), function(age) {
      element <- list_dist_spp[[age]]

      if ((is.matrix(element) || is.data.frame(element)) && nrow(element) > 0) {
        data.frame(
          species = rownames(element),
          age = age,
          mean_dist_to_cooccur = element[, 1],
          row.names = NULL
        )
      } else {
        # Return NA row if element is not a matrix/data.frame or has 0 rows
        data.frame(
          species = NA,
          age = age,
          mean_dist_to_cooccur = NA
        )
      }
    }))

    # Ensure the column is character first
    df_dist_spp$mean_dist_to_cooccur <- as.character(df_dist_spp$mean_dist_to_cooccur)

    # Create flag column for singleton entries
    df_dist_spp$is_singleton <- df_dist_spp$mean_dist_to_cooccur == "NA_singleton"

    # Replace "NA_singleton" with "0"
    df_dist_spp$mean_dist_to_cooccur[df_dist_spp$mean_dist_to_cooccur == "NA_singleton"] <- "0"

    # Replace "<NA>" strings with real NA
    df_dist_spp$mean_dist_to_cooccur[df_dist_spp$mean_dist_to_cooccur == "<NA>"] <- NA

    # Convert mean_dist_to_cooccur to numeric
    df_dist_spp$mean_dist_to_cooccur <- as.numeric(df_dist_spp$mean_dist_to_cooccur)


    return(df_dist_spp)

  }
