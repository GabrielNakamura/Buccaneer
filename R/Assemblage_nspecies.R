#' Calculate the number of species in each grid/assemblage/site
#'
#' @param df.TS.TE Data frame object containing at least four columns. Species names,
#'     origination time, extinction time and a trait value for each species.
#' @param df.occ Data frame object containing the occurrence records for each species.
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
#' @param Max.age Character indicating the name of the column containing the upper
#'     age limit for occurrence record.
#' @param Min.age Character indicating the name of the column containing the
#'     lower age limit for occurrence record.
#' @param site Character indicating the name of the column containing the
#'     information on assemblage location.
#'
#' @returns
#' @export
#'
#' @examples
assemblage_nspecies <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE",
           Max.age = "Max.age",
           Min.age = "Min.age",
           site = "site"){
    # subseting columns
    df.TS.TE <- df.TS.TE[, c(species, TS, TE)]
    colnames(df.TS.TE) <- c("species", "TS", "TE")


    df_occ <-
      df.occ[, c(species, Max.age, Min.age, site)]
    vars <- list(species, Max.age, Min.age, site)
    name_vars <- c("species", "Max.age", "Min.age", "site")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names
    df_occ$site <- as.factor(df_occ$site)

    # Generating time intervals used to compute temporal coexistence
    seq_interval <- c(seq(from = ceiling(max(df.TS.TE[, "TS"])),
                          to = ceiling(min(df.TS.TE[, "TE"])),
                          by = -time.slice))

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

    names(spp_slice) <- format(seq_interval, trim = TRUE, scientific = FALSE)


    # calculating data frame of species cooccurrence for site

    list_occur_site <-
      comp_site_occurrence(spp_slice = spp_slice, df.occ = df_occ)

    # calculating the number of species by site in each timeslice

    list_n_species_site <-
      lapply(list_occur_site, function(x){
        nspecies <- rowSums(x[, -1])
        names(nspecies) <- x$site
        nspecies
      })

    names(list_n_species_site) <- format(seq_interval, trim = TRUE, scientific = FALSE)

    df_long_res <-
      do.call(rbind, lapply(names(list_n_species_site), function(nm) {
        if(length(list_n_species_site[[nm]]) == 0){
          NA
        } else{
          data.frame(
            time.slice = nm,
            sites = names(list_n_species_site[[nm]]),
            n.species = unname(list_n_species_site[[nm]])
          )
        }
      })
      )

    return(df_long_res)
  }
