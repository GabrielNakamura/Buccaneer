#' Check occurrence records with crossing boundaries in time bins
#'
#' @param df.occ.fossil 
#' @param interval 
#' @param species 
#' @param Max.age 
#' @param Min.age 
#' @param remove.sub.species 
#'
#' @return
#' @export
#'
#' @examples
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