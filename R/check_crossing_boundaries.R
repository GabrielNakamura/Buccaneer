#' Check Occurrence Records Crossing Time Bin Boundaries
#'
#' This function analyzes how fossil occurrence records align with predefined
#' time bins by classifying each record based on whether its age uncertainty
#' (minimum and maximum ages) falls completely within a bin, crosses one boundary,
#' or spans multiple bins. This is useful for assessing data quality and understanding
#' the temporal resolution of fossil occurrences relative to time bins used in analyses.
#'
#' @param df.occ.fossil A data frame containing fossil occurrence records with at
#'     least three columns: species names, maximum (oldest) age estimates, and
#'     minimum (youngest) age estimates. Each row represents a single occurrence
#'     record with age uncertainty.
#' @param interval Numeric vector. A sequence of time bin boundaries in descending
#'     order (from oldest to youngest). For example, \code{c(100, 90, 80, 70)}
#'     creates three time bins: 100-90, 90-80, and 80-70. Must have at least
#'     two elements.
#' @param species Character. The name of the column in \code{df.occ.fossil}
#'     containing species identifiers. Default is "species".
#' @param Max.age Character. The name of the column in \code{df.occ.fossil}
#'     containing the maximum (oldest) age estimate for each occurrence record.
#'     Default is "Maximum_Age".
#' @param Min.age Character. The name of the column in \code{df.occ.fossil}
#'     containing the minimum (youngest) age estimate for each occurrence record.
#'     Default is "Minimum_Age".
#' @param remove.sub.species Logical. Should subspecies-level identifications be
#'     removed, keeping only species-level records? Default is TRUE. This parameter
#'     is passed to \code{clean_occ_fossil()}.
#'
#' @return A list containing two elements:
#'   \item{summary_table}{A named numeric vector showing the percentage of
#'       occurrence records in each category:
#'       \itemize{
#'         \item \code{within}: Records completely contained within a time bin
#'         \item \code{upper_only}: Records crossing only the upper (older) boundary
#'         \item \code{lower_only}: Records crossing only the lower (younger) boundary
#'         \item \code{both}: Records spanning an entire time bin (crossing both boundaries)
#'       }}
#'   \item{df_records_interval}{A data frame containing all occurrence records
#'       with added columns:
#'       \itemize{
#'         \item \code{boundaries}: Classification of boundary crossing (one of the
#'               four categories above)
#'         \item \code{time.interval}: The time bin identifier (e.g., "100-90")
#'         \item All original columns from the cleaned occurrence data
#'       }}
#'
#' @details
#' The function classifies each occurrence record into one of four categories:
#' \enumerate{
#'   \item \strong{within}: Both Max.age and Min.age fall within the time bin
#'         boundaries (Max.age <= upper boundary AND Min.age >= lower boundary)
#'   \item \strong{upper_only}: The record crosses the upper (older) boundary
#'         (Max.age > upper boundary BUT Min.age is within the bin)
#'   \item \strong{lower_only}: The record crosses the lower (younger) boundary
#'         (Min.age < lower boundary BUT Max.age is within the bin)
#'   \item \strong{both}: The record spans the entire time bin, crossing both
#'         boundaries (Max.age >= upper boundary AND Min.age <= lower boundary)
#' }
#'
#' Before classification, the function calls \code{clean_occ_fossil()} to:
#' \itemize{
#'   \item Calculate midpoint ages
#'   \item Optionally remove subspecies-level identifications
#'   \item Standardize column names
#' }
#'
#' Records crossing boundaries indicate age uncertainty that exceeds the
#' temporal resolution of the time bins. High percentages of crossing records
#' may suggest that finer time bins or different binning strategies should be
#' considered for the analysis.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example fossil occurrence data with age uncertainties
#' df_fossils <- data.frame(
#'   species = c("sp1", "sp2", "sp3", "sp4", "sp5"),
#'   Maximum_Age = c(95, 105, 88, 92, 110),
#'   Minimum_Age = c(85, 95, 82, 78, 70)
#' )
#'
#' # Define time bins (100-90, 90-80, 80-70)
#' time_bins <- c(100, 90, 80, 70)
#'
#' # Check boundary crossing patterns
#' results <- check_crossing_boundaries(
#'   df.occ.fossil = df_fossils,
#'   interval = time_bins
#' )
#'
#' # View summary of boundary crossing
#' results$summary_table
#'
#' # View detailed classification for all records
#' head(results$df_records_interval)
#'
#' # Filter records that span entire bins
#' spanning_records <- subset(
#'   results$df_records_interval,
#'   boundaries == "both"
#' )
#'
#' # Visualize boundary crossing patterns
#' library(ggplot2)
#' ggplot(results$df_records_interval,
#'        aes(x = time.interval, fill = boundaries)) +
#'   geom_bar(position = "fill") +
#'   labs(y = "Proportion", x = "Time Interval") +
#'   theme_minimal()
#'
#' # Calculate statistics by time interval
#' library(dplyr)
#' interval_summary <- results$df_records_interval %>%
#'   group_by(time.interval, boundaries) %>%
#'   summarise(count = n(), .groups = "drop")
#' }
check_crossing_boundaries <-
  function(df.occ.fossil,
           interval,
           species = "species",
           Max.age = "Maximum_Age",
           Min.age = "Minimum_Age",
           remove.sub.species = TRUE){

    # cleaning and processing data
    df.occ.fossil2 <-
      clean_occ_fossil(df.occ.fossil = df.occ.fossil,
                       method.ages = "midpoint",
                       species = species,
                       Max.age = Max.age,
                       Min.age = Min.age,
                       remove.sub.species = remove.sub.species)

    # identifying the temporal position of occurrence records
    list_record_intervals <- vector(mode = "list", length = length(interval)-1)
    for(i in 1:(length(interval) - 1)){
      # i = 1
      df_1 <- df.occ.fossil2[(df.occ.fossil2$Max.age <= interval[i]) & (df.occ.fossil2$Max.age >= interval[i + 1])
                             & (df.occ.fossil2$Min.age <= interval[i])
                             & (df.occ.fossil2$Min.age >= interval[i + 1]), ] # within the interval
      df_2 <- df.occ.fossil2[(df.occ.fossil2$Max.age > interval[i])
                             & (df.occ.fossil2$Min.age < interval[i])
                             & (df.occ.fossil2$Min.age > interval[i + 1]), ] # cross upper boundary
      df_3 <- df.occ.fossil2[df.occ.fossil2$Max.age < interval[i]
                             & df.occ.fossil2$Max.age > interval[i+1]
                             & df.occ.fossil2$Max.age < interval[i+1], ] # cross lower boundary
      df_4 <- df.occ.fossil2[df.occ.fossil2$Max.age >= interval[i]
                             & df.occ.fossil2$Min.age <= interval[i+1], ] # cross upper and lower boundary
      boundaries <- rep(x = c("within", "upper_only", "lower_only", "both"),
                        times = c(dim(df_1)[1], dim(df_2)[1], dim(df_3)[1], dim(df_4)[1]))
      df_all <- rbind(df_1, df_2, df_3, df_4)
      df_all2 <- data.frame(df_all, boundaries)
      df_all_res <-
        df_all2 |>
        dplyr::mutate(time.interval = paste(interval[i], interval[i + 1], sep = "-"))
      list_record_intervals[[i]] <- df_all_res
    }

    # binding lists by row
    df_record_intervals <- do.call(rbind, list_record_intervals)

    # Percentage of crossing occurrences
    crossing_boundaries <- table(df_record_intervals$boundaries)/sum(table(df_record_intervals$boundaries)) *100

    # list of results
    list_res <- vector(mode = "list")

    list_res$summary_table <- crossing_boundaries # table with summary data
    list_res$df_records_interval <- df_record_intervals # long data frame with all records and its classification

    return(list_res)

  }
