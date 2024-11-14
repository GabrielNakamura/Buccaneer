#' Cleaning and flagging fossil occurrence records according different criteria
#'
#' @param df.occ.fossil
#' @param method.ages
#' @param thresh.age.range
#' @param species
#' @param Max.age
#' @param Min.age
#' @param TS
#' @param TE
#' @param lat
#' @param lng
#' @param site
#' @param group
#' @param trait
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
clean_occ_fossil <-
  function(df.occ.fossil,
           method.ages = c("midpoint", "upper", "lower"),
           thresh.age.range = 10,
           species = "species",
           Max.age = "Maximum_Age",
           Min.age = "Minimum_Age",
           remove.sub.species = TRUE,
           comp.TS.TE = TRUE,
           lat = NULL,
           lng = NULL,
           site = NULL,
           group = NULL,
           trait = NULL){
    # subsetting the data frame with the variables provided
    df.occ.fossil <-
      df.occ.fossil[, c(species, Max.age, Min.age, lat, lng, site, group, trait)]
    vars <- list(species, Max.age, Min.age, lat, lng, site, group, trait)
    name_vars <- c("species", "Max.age", "Min.age", "lat", "lng", "site", "group", "trait")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    df.occ.fossil2 <- df.occ.fossil
    colnames(df.occ.fossil2) <- column.names

    # checking for numeric values in min and max age
    df.occ.fossil2 <-
      df.occ.fossil2 |>
      dplyr::mutate(Min.age = as.numeric(Min.age), Max.age = as.numeric(Max.age))

    # adding midpoint to data and flagging occurrence ranges
    df.occ.fossil3 <-
      df.occ.fossil2 |>
      dplyr::mutate(midpoint = (abs(Max.age + Min.age)/2)) |>
      dplyr::mutate(age.range = abs(Max.age - Min.age)) |>
      dplyr::mutate(flag.age.range = ifelse(age.range >= thresh.age.range, "TRUE", "FALSE"))

    # adding TS and TE based on midpoint for each species
    if(comp.TS.TE == TRUE){
      df.occ.fossil4 <-
        df.occ.fossil3 |>
        dplyr::group_by(species) |>
        dplyr::mutate(TS = max(Max.age), TE = min(Min.age))

      # looking for the presence of NA
      if(any(is.na(df.occ.fossil4$TS)) == TRUE | any(is.na(df.occ.fossil4$TE)) == TRUE){
        warning("There are NAs in TS and/or TE columns. These values should be
                  numeric and will be automatically removed")
        na_ts_te <- which(is.na(df.occ.fossil4$TS) == TRUE | is.na(df.occ.fossil4$TE) == TRUE)
        df.occ.fossil4 <- df.occ.fossil4[-na_ts_te, ]
      }
    } else{
      df_occ.fossil4 <- df.occ.fossil3
    }

    # detecting subspecies
    if(remove.sub.species == TRUE){
      df.occ.fossil5 <-
        df.occ.fossil4 |>
        dplyr::mutate(n.words = stringr::str_count(species, "\\S+")) |>
        dplyr::mutate(subspecies = ifelse(n.words >= 3, "subspecies", "species")) |>
        dplyr::select(-n.words)
    }

    return(df.occ.fossil5)

  }
