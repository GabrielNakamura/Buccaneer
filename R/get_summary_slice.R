#' Get Summary Information for Species Across Time Slices
#'
#' This function processes fossil occurrence data to generate summaries across
#' discrete time slices. It calculates species temporal ranges (origination and
#' extinction times) from occurrence records, identifies which species are present
#' in each time slice, counts occurrence records per species, and creates
#' species-by-species co-occurrence matrices. The function also flags occurrences
#' with large age uncertainties based on a user-defined threshold.
#'
#' @param df.TS.TE A data frame containing fossil occurrence records with at least
#'     three columns: species names, maximum (oldest) age estimates, and minimum
#'     (youngest) age estimates. Additional columns may include spatial coordinates,
#'     site information, group assignments, and trait values.
#' @param timeframe Numeric. The time interval (in millions of years or appropriate
#'     time units) between consecutive time slices. Negative values will create
#'     backwards intervals from the oldest to youngest occurrences.
#' @param method.ages Character. The method used to estimate species ages from
#'     occurrence records. Options include:
#'     \itemize{
#'       \item \code{"midpoint"} (default): Use the midpoint between max and min ages
#'       \item \code{"upper"}: Use the maximum (oldest) age
#'       \item \code{"lower"}: Use the minimum (youngest) age
#'     }
#' @param thresh.age.range Numeric. The threshold for flagging occurrence records
#'     with large age uncertainties. Records with age ranges (Max.age - Min.age)
#'     greater than or equal to this value are flagged. Default is 10.
#' @param species Character. The name of the column in \code{df.TS.TE} containing
#'     species identifiers. Default is "species".
#' @param Max.age Character. The name of the column in \code{df.TS.TE} containing
#'     the maximum (oldest) age estimate for each occurrence record. Default is "Maximum_Age".
#' @param Min.age Character. The name of the column in \code{df.TS.TE} containing
#'     the minimum (youngest) age estimate for each occurrence record. Default is "Minimum_Age".
#' @param TS Character. The name of the column containing origination (first
#'     appearance) times. Default is NULL. If NULL, TS is calculated as the
#'     maximum midpoint age for each species.
#' @param TE Character. The name of the column containing extinction (last
#'     appearance) times. Default is NULL. If NULL, TE is calculated as the
#'     minimum midpoint age for each species.
#' @param lat Character. The name of the column containing latitude coordinates.
#'     Default is NULL. If provided, latitude information is retained in output.
#' @param lng Character. The name of the column containing longitude coordinates.
#'     Default is NULL. If provided, longitude information is retained in output.
#' @param site Character. The name of the column containing site location identifiers.
#'     Default is NULL. If provided, site information is retained in output.
#' @param group Character. The name of the column containing group assignments
#'     for species (e.g., clade, family). Default is NULL. If provided, group
#'     information is retained in output.
#' @param trait Character. The name of the column containing trait values for
#'     species. Default is NULL. If provided, trait information is retained in output.
#' @param ... Additional arguments (currently not used but reserved for future extensions).
#'
#' @return A list containing two elements:
#'   \item{df_sub_bin}{A named list of data frames, one for each time slice.
#'       Each data frame contains all occurrence records present in that slice,
#'       with added columns:
#'       \itemize{
#'         \item \code{midpoint}: Numeric. The midpoint age of each occurrence
#'         \item \code{age.range}: Numeric. The age uncertainty (Max.age - Min.age)
#'         \item \code{flag.age.range}: Character. "TRUE" if age.range >= thresh.age.range,
#'               "FALSE" otherwise
#'         \item \code{TS}: Numeric. Species origination time (maximum midpoint)
#'         \item \code{TE}: Numeric. Species extinction time (minimum midpoint)
#'         \item \code{count_records}: Integer. Number of occurrence records for
#'               each species in that time slice
#'       }
#'       List names follow the format "slice_X" where X is the time slice value.}
#'   \item{coex_matrix}{A named list of binary co-occurrence matrices, one for
#'       each time slice. Each matrix:
#'       \itemize{
#'         \item Has dimensions: n × n, where n = total number of unique species
#'               in the dataset
#'         \item Contains 1 if a species is present in the time slice, 0 if absent
#'         \item Has row and column names as species identifiers
#'       }
#'       List names follow the format "slice_X".}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input data structure (must be a data frame)
#'   \item Calculates midpoint ages for each occurrence: (Max.age + Min.age) / 2
#'   \item Calculates age uncertainties: Max.age - Min.age
#'   \item Flags occurrences with age uncertainties >= thresh.age.range
#'   \item Computes species temporal ranges (TS and TE) from occurrence midpoints
#'   \item Creates time slices from oldest to youngest occurrence
#'   \item For each time slice, identifies species present (TS >= slice AND TE <= slice)
#'   \item Counts occurrence records per species per time slice
#'   \item Generates binary co-occurrence matrices for each time slice
#' }
#'
#' Age uncertainty handling:
#' \itemize{
#'   \item \strong{flag.age.range}: Identifies potentially problematic occurrences
#'         with large temporal uncertainties
#'   \item \strong{Threshold}: Adjustable via \code{thresh.age.range} parameter
#'   \item \strong{Use case}: Helps assess data quality and decide whether to
#'         filter or weight occurrences differently
#' }
#'
#' Missing value handling:
#' \itemize{
#'   \item NA values in TS or TE columns trigger a warning and are removed
#'   \item Non-numeric TS or TE values trigger an error
#' }
#'
#' This function is useful for preparing fossil occurrence data for downstream
#' analyses such as diversity dynamics, trait evolution, or biogeographic patterns
#' through time.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example fossil occurrence data
#' df_fossils <- data.frame(
#'   species = c("sp1", "sp2", "sp3", "sp1", "sp2", "sp4"),
#'   Maximum_Age = c(100, 95, 90, 88, 85, 80),
#'   Minimum_Age = c(95, 90, 85, 83, 80, 75),
#'   lat = c(10, 15, 20, 12, 18, 25),
#'   lng = c(-50, -55, -60, -52, -57, -65),
#'   site = c("A", "B", "C", "A", "B", "D"),
#'   group = c("G1", "G1", "G2", "G1", "G1", "G2"),
#'   trait = c(1.2, 2.5, 3.1, 1.2, 2.5, 4.0)
#' )
#'
#' # Get summary for 5 Ma time slices
#' results <- get_summary_slices(
#'   df.TS.TE = df_fossils,
#'   timeframe = 5,
#'   thresh.age.range = 10
#' )
#'
#' # View occurrences in the first time slice
#' head(results$df_sub_bin[[1]])
#'
#' # View co-occurrence matrix for the first time slice
#' results$coex_matrix[[1]]
#'
#' # Count species richness per time slice
#' richness_per_slice <- sapply(results$df_sub_bin, function(df) {
#'   length(unique(df$species))
#' })
#' richness_per_slice
#'
#' # Count total occurrences per time slice
#' records_per_slice <- sapply(results$df_sub_bin, nrow)
#' records_per_slice
#'
#' # Identify occurrences with large age uncertainties
#' flagged_occurrences <- do.call(rbind, results$df_sub_bin) %>%
#'   filter(flag.age.range == "TRUE")
#' head(flagged_occurrences)
#'
#' # Calculate sampling intensity through time
#' sampling_intensity <- data.frame(
#'   slice = names(results$df_sub_bin),
#'   n_records = records_per_slice,
#'   n_species = richness_per_slice,
#'   records_per_species = records_per_slice / richness_per_slice
#' )
#'
#' # Plot richness through time
#' library(ggplot2)
#' ggplot(sampling_intensity, aes(x = slice, y = n_species)) +
#'   geom_col() +
#'   labs(x = "Time Slice", y = "Species Richness") +
#'   theme_minimal() +
#'   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#'
#' # Use custom age method and threshold
#' results_upper <- get_summary_slices(
#'   df.TS.TE = df_fossils,
#'   timeframe = 5,
#'   method.ages = "upper",
#'   thresh.age.range = 5
#' )
#' }
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
    seq_interval <- seq(from = ceiling(max(df.TS.TE[, "TS"])),
                        to = ceiling(min(df.TS.TE[, "TE"])),
                        by = -time.slice)

    # defining species per slice and subsetting longevities data frame
    df_sub_slice <- vector(mode = "list", length = length(seq_interval))
    for (i in 1:length(seq_interval)){
      # i = 450
      df_sub_slice[[i]] <- df.TS.TE4[which(df.TS.TE4$TS >= seq_interval[i] & df.TS.TE4$TE <= seq_interval[i]), ]
    }

    # naming slices
    names(df_sub_slice) <- paste("slice", seq_interval, sep = "_")

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
