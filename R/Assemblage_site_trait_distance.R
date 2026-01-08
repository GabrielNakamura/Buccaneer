#' Compute mean distance for species co occurring in sites/grids
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
#' @param TS Character indicating the name of the columns of the data frame
#'     containing the information on origination time.
#' @param TE Character indicating the name of the column of the data frame
#'     containing the information on extinction time.
#' @param Max.age Character indicating the name of the column containing the upper age limit for occurrence record.
#' @param Min.age Character indicating the name of the column containing the lower age limit for occurrence record.
#' @param site Character indicating the name of the column containing the information on site location.
#'
#' @return A data frame with one row per site (assemblage) per time slice, containing:
#' \itemize{
#'   \item \code{site}: Site/assemblage identifier.
#'   \item \code{time.slice}: Time-slice label
#'   \item \code{mean.distance}: Mean trait distance among co-occurring species
#'       in that site and time slice (MPD if \code{nearest.taxon = "all"},
#'       MNND if \code{nearest.taxon = 1}, or based on the specified neighbor count).
#'   \item \code{var.distance}: Variance of those distances (NA if fewer than two species).
#'   \item \code{n.species} (if returned): Number of species in the site for that time slice.
#'   \item \code{n.pairs} (if returned): Number of pairwise comparisons used in
#'       the distance calculation.
#' }
#' @export
#'
#' @examples
assemblage_site_trait_distance <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           dist.trait,
           nearest.taxon = TRUE,
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
    df_occ$site <- as.factor(df_occ$site)

    # Generating time intervals used to compute temporal coexistence
    seq_interval <- seq(from = ceiling(max(df.TS.TE[, "TS"])),
                        to = ceiling(min(df.TS.TE[, "TE"])),
                        by = -time.slice)

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

    # calculating occurrence data frame
    list_matrix_occur_site <-
      comp_site_occurrence(spp_slice = spp_slice, df.occ = df_occ)

    # transforming in a matrix and removing the site column
    list_mat_comp_site <-
      lapply(list_matrix_occur_site, function(x){
        mat_comp <- x[, -1]
        mat_comp_matrix <- as.matrix(mat_comp)
        rownames(mat_comp_matrix) <- x$site
        return(mat_comp_matrix)
      })

    # calculating trait distance metrics
    if(nearest.taxon == TRUE){
      list_mpd_site <-
        lapply(list_mat_comp_site, function(x){
          if(length(x) == 0){
            NA
          } else{
            res_mpd_vector <-
              picante::ses.mntd(samp = x,
                                dis = matrix_dist_trait,
                                abundance.weighted = F)
            names(res_mpd_vector) <- rownames(x)
            return(res_mpd_vector)
          }
        })
      names(list_mpd_site) <- format(seq_interval, trim = TRUE, scientific = FALSE)
    } else{
      list_mpd_site <-
        lapply(list_mat_comp_site, function(x){
          if(length(x) == 0){
            NA
          } else{
            res_mpd_vector <-
              picante::ses.mpd(samp = x,
                               dis = matrix_dist_trait,
                               abundance.weighted = F)
            names(res_mpd_vector) <- rownames(x)
            return(res_mpd_vector)
          }
        })
      names(list_mpd_site) <- format(seq_interval, trim = TRUE, scientific = FALSE)
    }

    # organizing the long df with mpd results
    df_dist_site <-
      do.call(rbind, lapply(names(list_mpd_site), function(age) {
        element <- list_mpd_site[[age]]

        if (is.vector(element) && length(element) > 1) {
          data.frame(
            sites = names(element),
            time.slice = as.numeric(age),
            mean_dist_to_cooccur = element,
            row.names = NULL
          )
        } else {
          # Return NA row if element is not a vector or has length 0
          data.frame(
            sites = NA,
            time.slice = as.numeric(age),
            mean_dist_to_cooccur = NA
          )
        }
      }))

    return(df_dist_site)

  }


