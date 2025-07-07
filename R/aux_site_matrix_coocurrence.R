#' Auxiliary function to compute site cooccurrence matrix per timeslice
#'
#' @param spp_slice A list with species name in each timeslice
#' @param df.occ a occurrence data frame with the name of species and the sites in which
#'     they occur
#'
#' @returns
#' @export
#'
#' @examples
comp_site_cooccurr <-
  function(spp_slice, df.occ){

    # data frame with species per slice
    list_site_interval <-
      lapply(1:length(spp_slice), function(x){
        names_slices <- names(spp_slice)[[x]]
        names_species <- spp_slice[[x]]
        filtered_df <-
          df.occ %>% filter(species %in% names_species)

        filtered_df_site <-
          filtered_df |>
          filter(Max.age >= names_slices & Min.age <= names_slices)

        spp_coex_site <-
          filtered_df_site |>
          distinct(species, .keep_all = TRUE) |>
          select(species, TS, TE)

        filtered_df_site2 <-
          filtered_df_site |>
          group_by(site) |>
          distinct(species, .keep_all = T) |>
          add_count(site, name = "species.per.site") |>
          select(species, site, TS, TE, species.per.site)

        return(filtered_df_site2)
      })

    # calculate species cooccurrence matrix for each timeslice besed on site cooccurrence
    list_matrix_cooccur_site <-
      lapply(list_site_interval, function(x){

        # Step 1: Remove duplicates to avoid counting same species twice in a site
        df_unique <-
          x |>
          distinct(site, species)

        # Step 2: Create a presence/absence matrix (sites x species)
        site_species_matrix <-
          df_unique  |>
          mutate(present = 1) |>
          pivot_wider(names_from = species, values_from = present, values_fill = 0)

        # Step 3: Convert to a matrix and compute cross-product
        mat <- as.matrix(site_species_matrix[,-1])  # remove site column
        cooccur_matrix_sites <- t(mat) %*% mat  # species-by-species co-occurrence
        cooccur_matrix_sites
      })

    return(list_matrix_cooccur_site)
  }
