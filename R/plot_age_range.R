#' Plot Fossil Occurrence Records and Species Age Ranges
#'
#' This function creates a visualization of fossil occurrence records showing
#' the temporal distribution of species across geological time. Each species is
#' represented by horizontal line segments spanning from maximum (oldest) to
#' minimum (youngest) ages for each occurrence record. The plot can optionally
#' color-code species by number of occurrences and overlay geological time scale
#' boundaries.
#'
#' @param df.occ.fossil A data frame containing fossil occurrence records with at
#'     least three columns: species names, maximum (oldest) age estimates, and
#'     minimum (youngest) age estimates. Each row represents a single occurrence
#'     record with age uncertainty.
#' @param species Character. The name of the column in \code{df.occ.fossil}
#'     containing species identifiers. Default is "accepted_name".
#' @param max.age Character. The name of the column in \code{df.occ.fossil}
#'     containing the maximum (oldest) age estimate for each occurrence record.
#'     Default is "max_ma".
#' @param min.age Character. The name of the column in \code{df.occ.fossil}
#'     containing the minimum (youngest) age estimate for each occurrence record.
#'     Default is "min_ma".
#' @param occ.gradient Logical. Should line segments be colored by the number of
#'     occurrences for each species? If TRUE (default), species with more occurrence
#'     records will have different colors based on the chosen color scale.
#' @param log Logical. Should the number of occurrences be log-transformed for
#'     color scaling? Default is TRUE. This can help visualize species with highly
#'     variable occurrence counts.
#' @param name.scale.occ Character. The name of the viridis color scale to use
#'     for coloring by occurrence count. Options include "viridis" (default),
#'     "magma", "plasma", "inferno", "cividis", "mako", "rocket", and "turbo".
#' @param x.axis.name Character. Label for the x-axis. Default is "age range".
#' @param y.axis.name Character. Label for the y-axis. Default is "species.name".
#' @param show.species.name Logical. Should species names be displayed on the
#'     y-axis? Default is FALSE. When FALSE, y-axis labels are hidden to reduce
#'     clutter when plotting many species.
#' @param age.scheme Logical. Should geological time period boundaries be overlaid
#'     as vertical dashed lines? Default is TRUE. Currently uses North American
#'     Land Mammal Ages (NALMAs) boundaries.
#'
#' @return A ggplot2 object displaying the occurrence records. Species are arranged
#'     vertically (y-axis) sorted by their maximum age, with time (age) on the
#'     reversed x-axis (oldest to the left, youngest to the right). Each horizontal
#'     segment represents one occurrence record spanning its age uncertainty.
#'
#' @details
#' The function performs the following visualization steps:
#' \enumerate{
#'   \item Orders species by their oldest occurrence (maximum age)
#'   \item Counts the number of occurrence records per species
#'   \item Optionally log-transforms occurrence counts for better visualization
#'   \item Creates horizontal segments for each occurrence spanning min to max age
#'   \item Colors segments by occurrence count using viridis color scales
#'   \item Optionally overlays geological time boundaries
#' }
#'
#' Visualization features:
#' \itemize{
#'   \item \strong{X-axis}: Time in millions of years (Ma), reversed so older
#'         ages appear on the left
#'   \item \strong{Y-axis}: Species names (ordered by first appearance)
#'   \item \strong{Color gradient}: Number of occurrences per species (with or
#'         without log transformation)
#'   \item \strong{Vertical lines}: Geological time period boundaries (if enabled)
#' }
#'
#' The default geological time scheme uses North American Land Mammal Ages (NALMAs)
#' boundaries. To use custom time boundaries, set \code{age.scheme = FALSE} and
#' modify the returned ggplot object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create example fossil occurrence data
#' df_fossils <- data.frame(
#'   accepted_name = rep(c("sp1", "sp2", "sp3"), each = 3),
#'   max_ma = c(100, 95, 90, 85, 80, 75, 70, 65, 60),
#'   min_ma = c(95, 90, 85, 80, 75, 70, 65, 60, 55)
#' )
#'
#' # Basic plot with default settings
#' plot_occ_bins(df.occ.fossil = df_fossils)
#'
#' # Plot with species names visible
#' plot_occ_bins(
#'   df.occ.fossil = df_fossils,
#'   show.species.name = TRUE
#' )
#'
#' # Plot without log transformation of occurrence counts
#' plot_occ_bins(
#'   df.occ.fossil = df_fossils,
#'   log = FALSE
#' )
#'
#' # Plot without geological time boundaries
#' plot_occ_bins(
#'   df.occ.fossil = df_fossils,
#'   age.scheme = FALSE
#' )
#'
#' # Use different color scale
#' plot_occ_bins(
#'   df.occ.fossil = df_fossils,
#'   name.scale.occ = "plasma"
#' )
#'
#' # Customize the plot further
#' p <- plot_occ_bins(df.occ.fossil = df_fossils)
#' p +
#'   ggplot2::labs(title = "Fossil Occurrences Through Time") +
#'   ggplot2::theme_minimal()
#'
#' # Use custom column names
#' df_custom <- data.frame(
#'   taxon = rep(c("sp1", "sp2"), each = 2),
#'   oldest = c(100, 95, 90, 85),
#'   youngest = c(95, 90, 85, 80)
#' )
#'
#' plot_occ_bins(
#'   df.occ.fossil = df_custom,
#'   species = "taxon",
#'   max.age = "oldest",
#'   min.age = "youngest"
#' )
#' }
plot_occ_bins <-
  function(df.occ.fossil,
           species = "accepted_name",
           max.age = "max_ma",
           min.age = "min_ma",
           occ.gradient = TRUE,
           log = TRUE,
           name.scale.occ = "viridis",
           x.axis.name = "age range",
           y.axis.name = "species.name",
           show.species.name = FALSE,
           age.scheme = TRUE
  ){
    df_occ_fossil <- df.occ.fossil[, c(species, max.age, min.age)]
    names(df_occ_fossil) <- c("species", "max.age", "min.age")
    species <- df_occ_fossil$species


    if(occ.gradient == TRUE){
      df_occ_fossil <-
        df_occ_fossil |>
        dplyr::group_by(species) |>
        dplyr::arrange(desc(max.age)) |>                     # Step 1: Arrange within groups
        dplyr::mutate(max_value = max(max.age)) |>           # Step 2: Calculate max within groups
        dplyr::ungroup() |>
        dplyr::arrange(desc(max_value), species, desc(max.age)) |>
        dplyr::mutate(sequence = row_number()) |>
        dplyr::add_count(species, name = "n.occurrences")
    }

    if(log == TRUE){
      df_occ_fossil <-
        df_occ_fossil |>
        dplyr::group_by(species) |>
        dplyr::arrange(desc(max.age)) |>                     # Step 1: Arrange within groups
        dplyr::mutate(max_value = max(max.age)) |>           # Step 2: Calculate max within groups
        dplyr::ungroup()  |>
        dplyr::arrange(desc(max_value), species, desc(max.age)) |>
        dplyr::mutate(sequence = row_number()) |>
        dplyr::add_count(species, name = "n.occurrences") |>
        dplyr::mutate(n.occurrences = log(n.occurrences))
    }

    # ploting basic plot
    basic_plot <-
      ggplot2::ggplot(df_occ_fossil) +
      ggplot2::scale_x_reverse() +
      ggplot2::scale_y_discrete(limits = unique(df_occ_fossil$species)) +
      ggplot2::geom_segment(aes(x = max.age, xend = min.age, y = species, colour = n.occurrences), position = position_dodge(width = 0.05)) +
      ggplot2::scale_color_viridis_c(option = name.scale.occ, name = "n.occ")


    if(show.species.name == TRUE){
      basic_plot2 <-
        basic_plot +
        labs(title = " age ranges",
             x = "range",
             y = "species") +
        theme(axis.text.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
        )   # Adjust the size of y-axis title

    } else{
      basic_plot2 <-
        basic_plot +
        labs(title = "",
             x = "age",
             y = "occurrences") +
        theme(axis.text.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
        )
    }

    if(age.scheme == TRUE){
      dog_NALMAs_age <- c(37.2,33.9,33.3,30.8,20.43,15.97,13.6,10.3,4.9,1.8,0.3,0.0117,0)
      basic_plot2 +
        geom_vline(xintercept = dog_NALMAs_age, linetype = "dashed", color = "black", size = 0.6)
    } else{
      basic_plot2
    }
  }
