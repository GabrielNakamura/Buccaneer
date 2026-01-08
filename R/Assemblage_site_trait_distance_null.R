#' Compute Mean Pairwise Distances Between Species Cooccurring in Local Assemblages
#'
#' @param df.TS.TE Data frame object containing at least four columns. Species names,
#'     origination time, extinction time and a trait value for each species.
#' @param time.slice Numeric. The time interval (in the same units as TS and TE)
#'     between consecutive time slices for temporal binning.
#' @param dist_matrix_trait A distance matrix object (class \code{dist} or \code{matrix})
#'     containing pairwise trait distances between species. Row and column names
#'     must match species names in \code{df.TS.TE}. If NULL, distances will be
#'     computed from the trait column using Euclidean distance.
#' @param round.digits Integer. The number of decimal places to round time slice
#'     values. Default is 1. This affects temporal binning precision.
#' @param species Character. The name of the column in \code{df.TS.TE} and
#'     \code{df.occ} containing species identifiers. Default is "species".
#' @param TS Character. The name of the column in \code{df.TS.TE} containing
#'     origination (first appearance) times for each species. Default is "TS".
#' @param TE Character. The name of the column in \code{df.TS.TE} containing
#'     extinction (last appearance) times for each species. Default is "TE".
#' @param nearest.taxon Numeric or character. The number of nearest neighbors to
#'     consider when calculating mean distances. Use \code{1} for mean nearest
#'     neighbor distance (MNND), or \code{"all"} for mean pairwise distance (MPD).
#' @param nperm A scalar, indicating the number of permutations to be used in the
#'     null model
#' @param fixedmar Character, stating which row/column sums should be preserved.
#'     Options are "none", "rows", "columns" or "both"
#' @param mtype Character, indicanting the type of null model to be performed.
#'     \code{"prab"} is the option for presence-absence (detection/non-detection)
#'     matrix and \code{"count"} for matrix with count data.
#' @param df.occ A data frame containing fossil occurrence records with at least
#'     four columns: species names, minimum age, maximum age, and site location ID.
#'     Each row represents a single occurrence record at a specific site.
#' @param Max.age Character. The name of the column in \code{df.occ} containing
#'     the maximum (oldest) age estimate for each occurrence record. Default is "Max.age".
#' @param Min.age Character. The name of the column in \code{df.occ} containing
#'     the minimum (youngest) age estimate for each occurrence record. Default is "Min.age".
#' @param site Character. The name of the column in \code{df.occ} containing
#'     site location identifiers. Default is "site".
#'
#' @returns
#' @export
#'
#' @examples
assemblage_site_trait_distance_null <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           dist_matrix_trait,
           round.digits = 10,
           species = "species",
           TS = "TS",
           TE = "TE",
           Max.age = "Max.age",
           Min.age = "Min.age",
           site = "site",
           nearest.taxon = c("mpd", "mnnd"),
           nperm = 1000,
           fixedmar = "rows",
           mtype = "prab"){

    # subsetting TS TE matrix and checking procedures

    df.TS.TE <- df.TS.TE[, c(species, TS, TE)]


    colnames(df.TS.TE) <- c("species", "TS", "TE")

    # subsetting occurrence matrix
    df_occ <-
      df.occ[, c(species, Max.age, Min.age, site)]
    vars <- list(species, Max.age, Min.age, site)
    name_vars <- c("species", "Max.age", "Min.age", "site")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names
    df_occ$site <- as.factor(df_occ$site)


    # matrix coexistence at each timeslice
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE = df_longevities_canidae,
                               time.slice = 0.1,
                               round.digits = 1,
                               species = "species",
                               TS = "TS",
                               TE = "TE")

    # species composition at each timeslice
    spp_slice <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    seq_interval <- seq(from = ceiling(max(df_longevities_canidae[, "TS"])),
                        to = ceiling(min(df_longevities_canidae[, "TE"])),
                        by = -0.1)

    # naming time slices
    names(spp_slice) <- format(seq_interval, trim = TRUE, scientific = FALSE)

    # calculating occurrence in local assemblages
    list_occurrence <-
      comp_site_occurrence(spp_slice = spp_slice,
                           df.occ = df.occ,
                           species = "species",
                           Max.age = "max_T",
                           Min.age = "min_T",
                           site = "site.char")

    dist_matrix_trait <- as.matrix(dist_matrix_trait)
    rownames(dist_matrix_trait) <- df.TS.TE$species
    colnames(dist_matrix_trait) <- df.TS.TE$species


    # data frame to receive the results
    list_res_timeslice <-
      vector(mode = "list", length = length(list_occurrence))

    # beginning the computation of null model and ses metrics
    for(i in 1:length(list_occurrence)){
      # i = 12
      # compute null matrix for each timeslice composition
      if(dim(list_occurrence[[i]])[1] <= 1){
        list_res_timeslice[[i]] <- NA
      } else{
        null_occ_site <-
          vegan::permatfull(m = list_occurrence[[i]][, -1],
                            fixedmar = fixedmar,
                            mtype = mtype,
                            times = nperm)

        # computing distances
        if(nearest.taxon == "mpd"){
          list_matrix_coocccur_null_site <-
            lapply(null_occ_site$perm, function(x){
              picante::mpd(samp = x,
                           dis = dist_matrix_trait,
                           abundance.weighted = FALSE)
            })
        }
        if(nearest.taxon == "mnnd"){
          list_matrix_coocccur_null_site <-
            lapply(null_occ_site$perm, function(x){
              picante::mntd(samp = x,
                           dis = dist_body_mass,
                           abundance.weighted = FALSE)
            })
        }


        # joining all values of null mpd, columns are assemblages, rows are reps, only for one timeslice
        matrix_null_mpd <- do.call(rbind, list_matrix_coocccur_null_site)

        # getting only the observed value for that slice
        dist_obs_slice <-
          res_obs_assemblage_trait_site |>
          filter(time.slice == format(seq_interval[i],
                                      trim = TRUE,
                                      scientific = FALSE)) |>
          pull(mean_dist_to_cooccur) # the same timeslice used to calculate null matrix

        # mean value for all communities in the slice
        mean_null_slice <- apply(matrix_null_mpd, 2, mean)

        # variance for all communities in slice
        sd_null_slice <- apply(matrix_null_mpd, 2, sd)

        # standardized effect size for all slices - this is ses.mpd
        z_score_slice <- ((dist_obs_slice - mean_null_slice) / sd_null_slice)

        # absolute deviations of null
        dev_null <- abs(sweep(matrix_null_mpd, 2, mean_null_slice, "-"))

        # absolute deviations of observed
        dev_obs <- abs(dist_obs_slice - mean_null_slice)

        # count number of null >= observed for each column
        counts <- colSums(dev_null >= rep(dev_obs, each = nrow(matrix_null_mpd)))

        # two-tailed p-value
        p_two <- (counts + 1) / (nrow(matrix_null_mpd) + 1)

        # joining results for one timeslice in a dataframe
        df_res <-
          data.frame(dist.obs = dist_obs_slice,
                     dist.obs.z = z_score_slice,
                     p.value = p_two,
                     time.slice = names(list_occurrence[i]))

        # list to receive all results

        list_res_timeslice[[i]] <- df_res

      } # end of conditional
      print(i)
    }# end loop for each timeslice

    return(list_res_timeslice)
  }


