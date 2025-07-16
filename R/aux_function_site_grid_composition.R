#' Auxiliar function to compute site and grid composition
#'
#' @param spp_slice
#' @param df.occ
#' @param species
#' @param Max.age
#' @param Min.age
#' @param site
#'
#' @returns
#' @export
#'
#' @examples
comp_site_occurrence <-
  function(spp_slice,
           df.occ,
           species = "species",
           Max.age = "Max.age",
           Min.age = "Min.age",
           site = "site"){

    df_occ <-
      df.occ[, c(species, Max.age, Min.age, site)]
    vars <- list(species, Max.age, Min.age, site)
    name_vars <- c("species", "Max.age", "Min.age", "site")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names
    df_occ$site <- as.factor(df_occ$site)

    # data frame with species per slice
    list_site_interval <-
      # x = 1
      lapply(1:length(spp_slice), function(x){
        names_slices <- names(spp_slice)[[x]]
        names_species <- spp_slice[[x]]
        filtered_df <-
          df_occ |> filter(species %in% names_species)

        filtered_df_site <-
          filtered_df |>
          filter(Max.age >= as.numeric(names_slices) & Min.age <= as.numeric(names_slices))

        spp_coex_site <-
          filtered_df_site |>
          distinct(species, .keep_all = TRUE) |>
          select(species, Max.age, Min.age)

        filtered_df_site2 <-
          filtered_df_site |>
          group_by(site) |>
          distinct(species, .keep_all = T) |>
          add_count(site, name = "species.per.site") |>
          select(species, site, Max.age, Min.age, species.per.site)

        return(filtered_df_site2)
      })

    # calculate species cooccurrence matrix for each timeslice besed on site cooccurrence
    list_matrix_occurrence_site <-
      lapply(list_site_interval, function(x){

        # Step 1: Remove duplicates to avoid counting same species twice in a site
        df_unique <-
          x |>
          distinct(site, species)

        # Step 2: Create a presence/absence matrix (sites x species)
        site_species_matrix <-
          df_unique  |>
          mutate(present = 1) |>
          tidyr::pivot_wider(names_from = species, values_from = present, values_fill = 0)

        # Step 3: Convert to a matrix and compute cross-product
        # mat <- as.matrix(site_species_matrix[,-1])  # remove site column
        # cooccur_matrix_sites <- t(mat) %*% mat  # species-by-species co-occurrence
        # cooccur_matrix_sites
        site_species_matrix
      })

    return(list_matrix_occurrence_site)
  }
