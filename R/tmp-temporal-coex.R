df_longevities <- data.frame(species = rownames(longs[[1]]), TS = longs[[1]]$TS, TE = longs[[1]]$TE)
x = df_longevities
timeframe = NULL
intervals = dog_NALMAs_age
calc_temp_only_coex <- 
  function(x, intervals, timeframe = NULL, round.digits = 1){
    # Some checking functions
    if(!is.list(x) | !is.data.frame(x) == TRUE){
      stop("Data must be a list with a set of data frames or a single data frame")
    }
    
    
    # Generating time intervals used to compute temporal coexistence
    if(!is.null(timeframe) == TRUE){
      intervals <- seq(from = max(x[, "TS"]), to = min(x[, "TE"]), by = -timeframe)
      intervals <- round(intervals, round.digits)
    }
    
    # Computing species coexisting at each time interval
    
    list_living_intervals <- vector(mode = "list", length = (length(intervals) - 1))
    for(i in 1:(length(intervals) - 1)){
      # i = 1
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
      distinct(species, .keep_all = T) 
      
    # Organizing results
    list_res <- vector(mode = "list", length = 4)
    names(list_res) <- c("long.time.occurence", "summary.time.occurence", "dense.time.occurrence", "cooccurence.time.matrix")
    list_res$long.time.occurence <- df_long_intervals2
    list_res$summary.time.occurence <- df_summary_species
    list_res$dense.time.occurrence <- matrix_dense
    list_res$cooccurence.time.matrix <- list_all_coex
    
    return(list_res)
  }