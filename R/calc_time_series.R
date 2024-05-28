#' Get time series of species coexistence and species number in specific time frame 
#'
#' @param df.TS.TE 
#' @param intervals 
#' @param timeframe 
#' @param round.digits 
#' @param species 
#' @param TS 
#' @param TE 
#'
#' @return
#' @export
#'
#' @examples
calc_div_time_series <- 
  function(df.TS.TE, intervals, timeframe, round.digits = 1, species = "species", TS = "TS", TE = "TE"){
    
    # reorganizing data frame TS and TE
    df.TS.TE <- df.TS.TE[, c(species, TS, TE)]
    
    # Generating time intervals used to compute temporal coexistence
    seq_interval <- seq(from = max(df.TS.TE[, "TS"]), to = min(df.TS.TE[, "TE"]), by = -timeframe)
    seq_interval <- round(seq_interval, digits = round.digits)
    
    # defining species per bin and subsetting longevities data frame
    df_sub_bin <- vector(mode = "list", length = length(seq_interval))
    for (i in 1:length(seq_interval)){
      # i = 1
      df_sub_bin[[i]] <- df.TS.TE[which(df.TS.TE$TS >= seq_interval[i] & df.TS.TE$TE <= seq_interval[i]), ]
    }
    matrix_coex_bin <- vector(mode = "list", length = length(df_sub_bin))
    
    # list with species pairwise matrix
    list_matrix_coex <- 
      lapply(matrix_coex_bin, function(x){
      matrix(0, nrow = length(df.TS.TE$species), ncol = length(df.TS.TE$species), 
             dimnames = list(df.TS.TE$species, df.TS.TE$species)
             )
    })
    
    # pairwise matrix with co-occurrences
    for(i in 1:length(list_matrix_coex)){
      # i = 10
      list_matrix_coex[[i]][match(df_sub_bin[[i]]$species, rownames(list_matrix_coex[[i]])), 
                            match(df_sub_bin[[i]]$species, colnames(list_matrix_coex[[i]]))] <- 1
    }
    
    # decomposing the output in two. One with a data frame and another with a coexistence matrix
    list_res <- vector(mode = "list", length = 2)
    list_res$df_sub_bin <- df_sub_bin
    list_res$coex_matrix <- list_matrix_coex
    
    return(list_res)
  }