
#' Compute summary statistics from fossil record
#'
#' @param df.occ.fossil 
#' @param interval 
#' @param age.occ 
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

