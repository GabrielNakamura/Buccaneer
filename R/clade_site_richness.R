#' Compute time series based on mean site richness
#'
#' @param df.TS.TE Data frame object containing at least four columns. Species names,
#'     origination time, extinction time and a trait value for each species.
#' @param df.occ a data frame object containing the occurrence records for each species.
#'     This must have at least a column indicating the name of species, its minimum and maximum age estimate,
#'     and its site location ID.
#' @param time.slice Scalar indicating the time interval between consecutive time slices.
#' @param round.digits Scalar indicating the precision of time slices.
#' @param species Character indicating the name of the column of the data frame
#'     containing the species name information.
#' @param TS Character indicating the name of the columns of the data frame
#'     containing the information on origination time.
#' @param TE Character indicating the name of the column of the data frame
#'     containing the information on extinction time.
#' @param group Character indicating the name of the column that contain the groups that will be used
#'     in comparison.
#' @param group.focal.compare Character vector indicating the focal (first element) and comparison (second element)
#'     groups used in the calculation. If NULL, the default, the metrics  will be calculated
#'     using all  clades.
#' @param type.comparison Character. It can be "between" to compute distances only between species/genus of two groups
#'     or "within" to calculate distance only inside the focal group. If null the distance is computed
#'     considering all clades together
#' @param Max.age Character indicating the name of the column containing the upper age limit for occurrence record.
#' @param Min.age Character indicating the name of the column containing the lower age limit for occurrence record.
#' @param site Character indicating the name of the column containing the information on site location.
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
           TE = "TE",
           Max.age = "Max.age",
           Min.age = "Min.age",
           site = "site",
           group = NULL,
           group.focal.compare = NULL,
           type.comparison = NULL){

    # sub-setting the dataframe to keep or remove the group information

    if(!is.null(group) == TRUE){
      df.TS.TE <- df.TS.TE[, c(species, TS, TE, group)]
      colnames(df.TS.TE) <- c("species", "TS", "TE", "group")
    }
    else{
      df.TS.TE <- df.TS.TE[, c(species, TS, TE)]
      colnames(df.TS.TE) <- c("species", "TS", "TE")
    }

    # sub-setting and correcting the names of occurrence dataframe
    df_occ <-
      df.occ[, c(species, Max.age, Min.age, site)]
    vars <- list(species, Max.age, Min.age, site)
    name_vars <- c("species", "Max.age", "Min.age", "site")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names
    df_occ$site <- as.factor(df_occ$site)

    # Generating time intervals used to compute temporal coexistence
    seq_interval <- seq(from = max(df.TS.TE[, "TS"]), to = min(df.TS.TE[, "TE"]), by = -time.slice)
    seq_interval <- c(round(seq_interval, digits = round.digits))

    # Time coexistence matrix for all species
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = 1,
                               species = "species",
                               TS = "TS",
                               TE = "TE")

    # modified coexistence matrix containing group comparison
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

    # species composition at each timeslice
    spp_slice <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    names(spp_slice) <- seq_interval

    # calculating matrix of species cooccurrence for site

    list_matrix_cooccur_site <-
      comp_site_cooccurr(spp_slice = spp_slice, df.occ = df_occ)

    # calculating mean coexistence and variance in coexistence for each species in each timeslice

    mean_species_site_coex <-
      lapply(list_matrix_cooccur_site,
             function(x){
               mean(which((rowSums(x) -1) != 0))
             })

    var_species_site_coex <-
      lapply(list_matrix_cooccur_site,
             function(x){
               var(which((rowSums(x) -1) != 0))
             })

    # function output

    df_res <-
      data.frame(mean.coexistence = unlist(mean_species_site_coex),
                 var.distance = unlist(var_species_site_coex),
                 time.slice = seq_interval)

    return(df_res)


  }
