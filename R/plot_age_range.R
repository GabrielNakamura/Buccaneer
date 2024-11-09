#' Plot occurrence records and age ranges
#'
#' @param df.occ.fossil 
#' @param species 
#' @param max.age 
#' @param min.age 
#' @param occ.gradient 
#' @param log 
#' @param name.scale.occ 
#' @param x.axis.name 
#' @param y.axis.name 
#' @param show.species.name 
#'
#' @return
#' @export
#'
#' @examples
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