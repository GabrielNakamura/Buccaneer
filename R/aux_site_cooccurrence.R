#' Calculate site cooccurrence
#'
#' @param list_site_interval
#'
#' @returns
#' @export
#'
#' @examples
calc_cooccur_site_matrix <-
  function(spp_slice, df.occ){

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

    species_matrix <-
      list_site_interval |>
      distinct(species, site) |>
      mutate(present = 1) |>
      tidyr::pivot_wider(names_from = site, values_from = present, values_fill = 0)

    # Remove species column to get a numeric matrix
    pa_matrix <- as.matrix(species_matrix[,-1])
    rownames(pa_matrix) <- species_matrix$species

    # Co-occurrence matrix: number of shared sites between species
    cooccurrence_matrix_sites <- pa_matrix %*% t(pa_matrix)
    cooccurrence_matrix_sites_pa <- ifelse(cooccurrence_matrix_sites >= 1, 1, 0)
    return(cooccurrence_matrix_sites_pa)
  }
