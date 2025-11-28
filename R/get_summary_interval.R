#' Compute Summary Statistics from Fossil Records Across Time Intervals
#'
#' This function calculates summary statistics for fossil occurrence records
#' grouped into discrete time intervals (bins). For each time interval, it
#' computes the number of occurrence records and the number of unique species
#' (taxonomic richness). The function bins occurrences based on their age
#' estimates and provides both record-level and species-level summaries.
#'
#' @param df.occ.fossil A data frame containing fossil occurrence records with
#'     columns for species names, age estimates (maximum and minimum ages, or
#'     midpoint), and optionally spatial coordinates, site information, group
#'     assignments, and trait values.
#' @param interval Numeric vector. A sequence of time bin boundaries in descending
#'     order (from oldest to youngest). For example, \code{c(100, 90, 80, 70)}
#'     creates three time bins: 100-90, 90-80, and 80-70. Must have at least
#'     two elements.
#' @param age.occ Character. The name of the column in \code{df.occ.fossil}
#'     containing the age estimate used for binning occurrences. Typically
#'     "midpoint" for the midpoint between maximum and minimum ages. Default is "midpoint".
#' @param species Character. The name of the column in \code{df.occ.fossil}
#'     containing species identifiers. Default is "species".
#' @param Max.age Character. The name of the column in \code{df.occ.fossil}
#'     containing the maximum (oldest) age estimate for each occurrence record.
#'     Default is "Maximum_Age".
#' @param Min.age Character. The name of the column in \code{df.occ.fossil}
#'     containing the minimum (youngest) age estimate for each occurrence record.
#'     Default is "Minimum_Age".
#' @param TS Character. The name of the column containing origination (first
#'     appearance) times. Default is NULL. This parameter is retained for future
#'     functionality but not currently used in calculations.
#' @param TE Character. The name of the column containing extinction (last
#'     appearance) times. Default is NULL. This parameter is retained for future
#'     functionality but not currently used in calculations.
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
#' @return A list containing two data frames:
#'   \item{df_summary_occurrences}{A data frame with all occurrence records
#'       augmented with two columns:
#'       \itemize{
#'         \item \code{time.interval}: Character. The time bin identifier (e.g., "100-90")
#'         \item \code{n.records.bin}: Integer. The total number of occurrence
#'               records in that time interval
#'       }
#'       This data frame has one row per occurrence record.}
#'   \item{df_summary_richness}{A data frame with unique species per time interval
#'       augmented with two columns:
#'       \itemize{
#'         \item \code{time.interval}: Character. The time bin identifier
#'         \item \code{n.species.bin}: Integer. The number of unique species
#'               (taxonomic richness) in that time interval
#'       }
#'       This data frame has one row per unique species per time interval.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates that input is a data frame
#'   \item Extracts specified columns from the occurrence data
#'   \item For each time interval (consecutive pairs in \code{interval}):
#'         \itemize{
#'           \item Filters occurrences where \code{age.occ} falls within the interval
#'           \item Assigns a time interval label (e.g., "100-90")
#'         }
#'   \item Calculates summary statistics:
#'         \itemize{
#'           \item \strong{n.records.bin}: Total number of occurrences per interval
#'           \item \strong{n.species.bin}: Number of unique species per interval
#'         }
#' }
#'
#' Binning logic:
#' \itemize{
#'   \item An occurrence is assigned to an interval if: \code{interval[i] >= age.occ >= interval[i+1]}
#'   \item Intervals are labeled using the format "upper-lower" (e.g., "100-90")
#'   \item Occurrences falling outside all intervals are excluded
#' }
#'
#' The two output data frames serve different purposes:
#' \itemize{
#'   \item \strong{df_summary_occurrences}: Retains all individual occurrence
#'         records for detailed analyses (e.g., spatial patterns, within-interval variation)
#'   \item \strong{df_summary_richness}: Provides species-level summaries for
#'         diversity analyses (e.g., richness through time, turnover rates)
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example fossil occurrence data
#' df_fossils <- data.frame(
#'   species = c("sp1", "sp2", "sp3", "sp1", "sp2", "sp4"),
#'   midpoint = c(95, 92, 88, 85, 82, 78),
#'   Maximum_Age = c(100, 95, 90, 90, 85, 80),
#'   Minimum_Age = c(90, 89, 86, 80, 79, 76),
#'   lat = c(10, 15, 20, 25, 30, 35),
#'   lng = c(-50, -55, -60, -65, -70, -75),
#'   site = c("A", "B", "C", "A", "B", "C"),
#'   group = c("G1", "G1", "G2", "G1", "G1", "G2"),
#'   trait = c(1.2, 2.5, 3.1, 1.2, 2.5, 4.0)
#' )
#'
#' # Define time intervals
#' time_bins <- c(100, 90, 80, 70)
#'
#' # Get summary statistics
#' results <- get_summary_interval(
#'   df.occ.fossil = df_fossils,
#'   interval = time_bins,
#'   age.occ = "midpoint"
#' )
#'
#' # View occurrence-level summary
#' head(results$df_summary_occurrences)
#'
#' # View species-level summary (richness)
#' head(results$df_summary_richness)
#'
#' # Plot richness through time
#' library(ggplot2)
#' richness_summary <- results$df_summary_richness %>%
#'   group_by(time.interval) %>%
#'   summarise(richness = first(n.species.bin))
#'
#' ggplot(richness_summary, aes(x = time.interval, y = richness)) +
#'   geom_col() +
#'   labs(x = "Time Interval", y = "Species Richness") +
#'   theme_minimal()
#'
#' # Calculate sampling intensity (records per species)
#' sampling_summary <- results$df_summary_occurrences %>%
#'   group_by(time.interval) %>%
#'   summarise(
#'     total_records = first(n.records.bin),
#'     n_species = n_distinct(species),
#'     records_per_species = total_records / n_species
#'   )
#'
#' # Use custom column names
#' df_custom <- data.frame(
#'   taxon = c("sp1", "sp2", "sp3"),
#'   age = c(95, 85, 75),
#'   oldest = c(100, 90, 80),
#'   youngest = c(90, 80, 70)
#' )
#'
#' results_custom <- get_summary_interval(
#'   df.occ.fossil = df_custom,
#'   interval = c(100, 90, 80, 70),
#'   age.occ = "age",
#'   species = "taxon",
#'   Max.age = "oldest",
#'   Min.age = "youngest"
#' )
#' }
get_summary_interval <-
  function(df.occ.fossil,
           interval,
           age.occ = "midpoint",
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
    if(!is.list(df.occ.fossil) | !is.data.frame(df.occ.fossil) == TRUE){
      stop("Data must be a list with a set of data frames or a single data frame")
    }

    # filtering columns
    df.occ.fossil <-
      df.occ.fossil[, c(age.occ, species, Max.age, Min.age, lat, lng, site, group, trait)]
    df.occ.fossil2 <- df.occ.fossil
    colnames(df.occ.fossil2) <- c("age.occ", "species", "Max.age", "Min.age")


    # Computing species coexisting at each time interval

    list_record_intervals <- vector(mode = "list", length = (length(interval) - 1))
    for(i in 1:(length(interval) - 1)){
      # i = 1
      df_1 <- df.occ.fossil2[(df.occ.fossil2$age.occ <= interval[i]) & (df.occ.fossil2$age.occ >= interval[i + 1]), ]
      df_res <-
        df_1 |>
        dplyr::mutate(time.interval = paste(interval[i], interval[i + 1], sep = "-"))
      list_record_intervals[[i]] <- df_res
    }

    # binding lists by row
    df_record_intervals <- do.call(rbind, list_record_intervals)

    # summary stats - number of records and number of species per interval
    df_record_intervals2 <-
      df_record_intervals |>
      dplyr::group_by(time.interval) |>
      dplyr::add_count(time.interval, name = "n.records.bin")

    df_record_intervals_richness <-
      df_record_intervals2 |>
      dplyr::ungroup() |>
      dplyr::group_by(time.interval) |>
      dplyr::distinct(species) |>
      dplyr::add_count(time.interval, name = "n.species.bin") |>
      dplyr::distinct(species, .keep_all = T)

    # Organizing results
    list_res <- vector(mode = "list", length = 2)
    names(list_res) <- c("df_summary_occurrences", "df_summary_richness")
    list_res$df_summary_occurrences <- df_record_intervals2
    list_res$df_summary_richness <- df_record_intervals_richness

    return(list_res)
  }

