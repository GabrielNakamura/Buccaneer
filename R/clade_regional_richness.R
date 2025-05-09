#' Calculates regional clade richness for lineages
#'
#' @param df.TS.TE 
#' @param time.slice 
#' @param round.digits 
#' @param species 
#' @param TS 
#' @param TE 
#'
#' @return
#' @export
#'
#' @examples
regional_clade_richness <- 
  function(df.TS.TE,
           time.slice, 
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE"){
    df.TS.TE <- df.TS.TE[, c(species, TS, TE)]
    
    # Generating time intervals used to compute temporal coexistence
    seq_interval <- seq(from = max(df.TS.TE[, "TS"]), to = min(df.TS.TE[, "TE"]), by = -time.slice)
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
    
    # just naming columns with slice names
    names(df_sub_bin) <- seq_interval
    df_sub_bin2 <- 
      lapply(names(df_sub_bin), function(name) {
        df <- df_sub_bin[[name]]
        df$time.slice <- name  # Add the name of the element as the new column
        return(df)
      })
    
    names(df_sub_bin2) <- names(df_sub_bin)
    
    # calculating lineage metrics 
    
    list_richness <- 
      lapply(df_sub_bin2, function(x){
        length(unique(x$species))
      })
    
    # dataframe with richness
    df_richness_slice <- data.frame(richness = unlist(list_richness), time.slice = names(df_sub_bin))
    
    return(df_richness_slice)
    
  }