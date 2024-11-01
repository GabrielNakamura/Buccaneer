
#' Getting summary information for time bins
#'
#' @param df.TS.TE
#' @param timeframe
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
get_summary_slices <-
  function(df.TS.TE,
           timeframe,
           method.ages = c("midpoint", "upper", "lower"),
           thresh.age.range = 10,
           species = "species",
           Max.age = "Maximum_Age",
           Min.age = "Minimum_Age",
           TS = NULL,
           TE = NULL,
           lat = NULL,
           lng = NULL,
           site = NULL,
           group = NULL,
           trait = NULL, ...){

    # Some checking functions
    if(!is.list(df.TS.TE) | !is.data.frame(df.TS.TE) == TRUE){
      stop("Data must be a list with a set of data frames or a single data frame")
    }

    # filtering columns
    df.TS.TE <-
      df.TS.TE[, c(species, Max.age, Min.age, lat, lng, site, group, trait)]
    df.TS.TE2 <- df.TS.TE
    colnames(df.TS.TE2) <- c("species", "Max.age", "Min.age")

    # adding midpoint to data and flagging occurrence ranges
    df.TS.TE3 <-
      df.TS.TE2 |>
      dplyr::mutate(midpoint = (abs(Max.age + Min.age)/2)) |>
      dplyr::mutate(age.range = abs(Max.age - Min.age)) |>
      dplyr::mutate(flag.age.range = ifelse(age.range >= thresh.age.range, "TRUE", "FALSE"))

    # adding TS and TE based on midpoint for each species
    df.TS.TE4 <-
      df.TS.TE3 |>
      dplyr::group_by(species) |>
      dplyr::mutate(TS = max(midpoint), TE = min(midpoint))

    # looking for the presence of NA
    if(any(is.na(df.TS.TE4$TS)) == TRUE | any(is.na(df.TS.TE4$TE)) == TRUE){
      warning("There are NAs in TS and/or TE columns. These values should be
                  numeric and will be automatically removed")
      na_ts_te <- which(is.na(df.TS.TE4$TS) == TRUE | is.na(df.TS.TE4$TE) == TRUE)
      df.TS.TE4 <- df.TS.TE4[-na_ts_te, ]
    }

    # filtering for only numeric
    if(!is.numeric(df.TS.TE4$TS) == TRUE  | !is.numeric(df.TS.TE4$TE) == TRUE){
      stop("TS and/or TE are not numeric values")
    }


    # Generating time intervals used to compute temporal coexistence
    seq_slice <- seq(from = max(df.TS.TE4[, "TS"]), to = min(df.TS.TE4[, "TE"]), by = -timeframe)
    seq_slice <- round(seq_interval, digits = round.digits)

    # defining species per slice and subsetting longevities data frame
    df_sub_slice <- vector(mode = "list", length = length(seq_slice))
    for (i in 1:length(seq_slice)){
      # i = 450
      df_sub_slice[[i]] <- df.TS.TE4[which(df.TS.TE4$TS >= seq_slice[i] & df.TS.TE4$TE <= seq_slice[i]), ]
    }

    # naming slices
    names(df_sub_slice) <- paste("slice", seq_slice, sep = "_")

    # adding occurrence counts for each bin
    df_sub_slice2 <-
      lapply(df_sub_slice, function(x){
        x |>
          dplyr::group_by(species) |>
          dplyr::add_count(species, name = "count_records")
      })

    # creating a long data frame with occurrence counts per bin
    df_all_slice <- purrr::imap(df_sub_slice2, ~ dplyr::mutate(.x, slice = .y))
    df_long_slice <- do.call(rbind, df_all_slice)
    df_count_record_bin <-
      df_long_slice |>
      dplyr::ungroup() |>
      dplyr::distinct(slice, count_records)

    # creating a long data frame with richness per bin
    df_rich_slice <-
      df_long_slice |>
      dplyr::group_by(slice) |>
      dplyr::distinct(species) |>
      dplyr::add_count(slice, name = "n_species_slice") |>
      dplyr::distinct(slice, n_species_slice)

    # pairwise matrix with co-occurrences
    list_co_occur <- vector(mode = "list", length = length(df_sub_slice2))
    for(i in 1:length(list_co_occur)){
      # i = 10
      matrix_tmp <- matrix(0, nrow = length(unique(df.TS.TE4$species)), ncol = length(unique(df.TS.TE4$species)),
                           dimnames = list(unique(df.TS.TE4$species), unique(df.TS.TE4$species))
      )
      matrix_tmp[unique(match(df_sub_slice2[[i]]$species, rownames(matrix_tmp))),
                 unique(match(df_sub_slice2[[i]]$species, colnames(matrix_tmp)))] <- 1
      list_co_occur[[i]] <- matrix_tmp
    }

    names(list_co_occur) <- paste("slice", seq_interval, sep = "_")

    # decomposing the output in two. One with a data frame and another with a coexistence matrix
    list_res <- vector(mode = "list", length = 2)
    list_res$df_sub_bin <- df_sub_slice
    list_res$coex_matrix <- list_co_occur

    return(list_res)
  }
