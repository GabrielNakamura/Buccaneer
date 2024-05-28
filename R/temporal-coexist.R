#' Calculation of temporal coexistence based on time of speciation and time of extinction
#' 
#' @param x Data frame with three columns. One with the name of species (species). Speciation time (TS) and Extinction time (TE).
#'     If columns names are different from the ones presented here, please, inform the names in column.names argument.
#' @param timeframe Scalar, indicating the time interval in which temporal coexistence will be calculated
#' @param round.digits Scalar indicating the rounding digits of ages used in intervals
#' @param intervals Numeric vector indicating the time slices in which the temporal coexistence will be computed
#' @param species Character indicating the name of the column containing species names, default is "species"
#' @param TS Character indicating the name of the column containing the speciation times for species, default is "TS"
#' @param TE Character indicating the name of the column containing the extinction times for species, default is "TE"
#'
#' @return
#' @export
#'
#' @examples
#' df_test <- 
#' sim.df.pyrate(n0 = 1,
#'               lambda = 0.1,
#'               mu = 0.04,
#'               tMax = 50, 
#'               seed = 1,
#'               rho = 1,
#'               bins = seq(100, 0, -1),
#'               returnAll = TRUE
#' )
#' df_test2 <- 
#' df_test %>% 
#'   distinct(TS, TE, .keep_all = TRUE) %>% 
#'   select(Species, TS, TE)
#' 
#' calc_temp_only_coex(x, timeframe = 10, species = "Species") 
#' 
calc_temp_only_coex <- 
  function(x, intervals, timeframe = NULL, round.digits = 1, species = "species", TS = "TS", TE = "TE"){
    # Some checking functions
    if(!is.list(x) | !is.data.frame(x) == TRUE){
      stop("Data must be a list with a set of data frames or a single data frame")
    }
    
    x <- x[, c(species, TS, TE)]
    # Generating time intervals used to compute temporal coexistence
    if(!is.null(timeframe) == TRUE){
      intervals <- seq(from = max(x[, "TS"]), to = min(x[, "TE"]), by = -timeframe)
      intervals <- round(intervals, round.digits)
    }
    
    # Computing species coexisting at each time interval
    
    list_living_intervals <- vector(mode = "list", length = (length(intervals) - 1))
    for(i in 1:(length(intervals) - 1)){
      # i = 12
      df_1 <- x[(x$TS <= intervals[i]) & (x$TS >= intervals[i + 1]) & (x$TE <= intervals[i]) & (x$TE >= intervals[i + 1]), ] 
      df_2 <- x[(x$TS > intervals[i]) & (x$TE < intervals[i]) & (x$TE > intervals[i + 1]), ]
      df_3 <- x[x$TS < intervals[i] & x$TS > intervals[i+1] & x$TE < intervals[i+1], ]
      df_4 <- x[x$TS >= intervals[i] & x$TE <= intervals[i+1], ]
      df_all <- rbind(df_1, df_2, df_3, df_4)
      df_all_res <- 
        df_all %>% 
        distinct(species, .keep_all = T) %>% 
        mutate(time.interval = paste(intervals[i], intervals[i + 1], sep = "-"))
      list_living_intervals[[i]] <- df_all_res
    }
    
    df_long_intervals <- do.call(rbind, list_living_intervals)
    
    # long format for temporal data
    df_long_intervals2 <- 
      df_long_intervals %>% 
      group_by(time.interval) %>% 
      add_count(name = "n.species.interval") %>% 
      ungroup() %>% 
      group_by(species) %>% 
      add_count(name = "species.freq.all") %>% 
      ungroup() %>% 
      mutate(longevities = TS - TE)
    
    # dense matrix of species occurrence by each time interval 
    matrix_dense <- as.matrix(phyloregion::long2sparse(x = df_long_intervals2, grids = "time.interval", species = "species")) 
    
    # coexistence matrix for each time interval
    list_all_coex <- lapply(1:nrow(matrix_dense), function(f) as.matrix(matrix_dense[f, ]) %*% t(matrix_dense[f, ]))
    names(list_all_coex) <- sort(rownames(matrix_dense), decreasing = T)
    
    # summary table for time 
    df_summary_species <- 
      df_long_intervals2 %>% 
      distinct(species, .keep_all = T) %>% 
      select(species, TS, TE, n.species.interval, species.freq.all, longevities)
    
    # Organizing results
    list_res <- vector(mode = "list", length = 4)
    names(list_res) <- c("long.time.occurence", "summary.time.occurence", "dense.time.occurrence", "cooccurence.time.matrix")
    list_res$long.time.occurence <- df_long_intervals2
    list_res$summary.time.occurence <- df_summary_species
    list_res$dense.time.occurrence <- matrix_dense
    list_res$cooccurence.time.matrix <- list_all_coex
    
    return(list_res)
  }