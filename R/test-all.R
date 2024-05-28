
source(here::here("R", "aux_sim_df_pyrate.R"))

tMax <- 100
df_test <- 
  sim.df.pyrate(n0 = 1,
           lambda = 0.1,
           mu = 0.05,
           tMax = tMax, 
           seed = 1,
           rho = 1,
           bins = seq(tMax, 0, -3),
           returnAll = F
  )

df_test2 <- 
  df_test %>% 
  distinct(TS, TE, .keep_all = TRUE) %>% 
  select(Species, TS, TE)

#x = df_test2
#timeframe = 0.1
#intervals = NULL

# calculating temporal coexistence only 
source(here::here("R", "temporal-coexist.R"))
test_res <- calc_temp_only_coex(x = df_test2, timeframe = 0.1)



# testing gaps and time intervals composition -----------------------------

time_intervals <- c(seq(105, 0, by = - 10), 0)

x = df_test2
timeframe
intervals = time_intervals
round.digits = 1

calc_temp_only_coex <- 
  function(x, timeframe, intervals = NULL, round.digits = 1){
    # Some checking functions
    if(!is.list(x) | !is.data.frame(x) == TRUE){
      stop("Data must be a list with a set of data frames or a single data frame")
    }
    
    # When the time intervals are provided by the user
    if(!is.null(intervals) == TRUE){
      intervals <- sort(intervals, decreasing = TRUE)
    }
    
    # Converting NAs to zeroes
    x2 <- 
      x %>% 
      mutate(TE = ifelse(is.na(TE) == TRUE, 0, TE))
    
    
    # Generating time intervals used to compute temporal coexistence
    if(is.null(intervals) == TRUE){
      intervals <- seq(from = max(x2[, "TS"]), to = min(x2[, "TE"]), by = -timeframe)
      intervals <- round(intervals, round.digits)
    }
    
    # Computing species coexisting at each time interval
    calc_temp_only_coex <- 
  function(x, timeframe, intervals = NULL, round.digits = 1){
    # Some checking functions
    if(!is.list(x) | !is.data.frame(x) == TRUE){
      stop("Data must be a list with a set of data frames or a single data frame")
    }
    
    # When the time intervals are provided by the user
    if(!is.null(intervals) == TRUE){
      intervals <- sort(intervals, decreasing = TRUE)
    }
    
    # Converting NAs to zeroes
    x2 <- 
      x %>% 
      mutate(TE = ifelse(is.na(TE) == TRUE, 0, TE))
    
    
    # Generating time intervals used to compute temporal coexistence
    if(is.null(intervals) == TRUE){
      intervals <- seq(from = max(x2[, "TS"]), to = min(x2[, "TE"]), by = -timeframe)
      intervals <- round(intervals, round.digits)
    }
    
    # Computing species coexisting at each time interval
    
    vec <- list(type = "list")
    for(i in 1:length(intervals)){
      # i = 8
      FL <- (x2$TS <= intervals[i]) & (x2$TS >= intervals[i+1]) & (x2$TE <= intervals[i]) & (x2$TE >= intervals[i+1])
      bL <- (df_test$TS > intervals[i]) & (df_test$TE < intervals[i]) & (df_test$TE > intervals[i+1])
      Ft <- (df_test$TS < intervals[i]) & (df_test$TS > intervals[i+1]) & (df_test$TE < intervals[i+1])
      bt <- (df_test$TS >= intervals[i]) & (df_test$TE <= intervals[i+1])
      vec[[i]] <- res1
    }
    
    
    list_spp_coex <- 
      lapply(intervals, function(f){
        sp_names_coex <- x2$Species[which(x2$TS >= f & x2$TE <= f)]
        return(sp_names_coex)
      })
    
    # naming time slices
    names(list_spp_coex) <- paste("t", as.character(round(intervals, digits = round.digits)), sep = "_")
    names_interval_long <- 
      as.character(rep(round(intervals, digits = round.digits),
                       unlist(lapply(list_spp_coex, function(f) length(f)))))
    # data frame in long format with species per time interval
    df_long_coex <- data.frame(species = unlist(list_spp_coex), time.interval = names_interval_long)
    
    # just counting the number of species per time slice
    df_summary_coex <- 
      df_long_coex %>% 
      group_by(time.interval) %>% 
      dplyr::mutate(n.spp = dplyr::n()) %>% 
      dplyr::distinct(n.spp)
    
    # dense matrix of species occurrence by each time interval 
    matrix_dense <- as.matrix(phyloregion::long2sparse(x = df_long_coex, grids = "time.interval", species = "species")) 
    
    # coexistence matrix for each time interval
    list_all_coex <- lapply(1:nrow(matrix_dense), function(f) as.matrix(matrix_dense[f, ]) %*% t(matrix_dense[f, ]))
    names(list_all_coex) <- paste("interval", rownames(matrix_dense), sep = "_")
    
    # Organizing results
    list_res <- vector(mode = "list", length = 3)
    names(list_res) <- c("summary.time.occurence", "dense.time.occurrence", "cooccurence.time.matrix")
    list_res$summary.time.occurence <- df_summary_coex
    list_res$dense.time.occurrence <- matrix_dense
    list_res$cooccurence.time.matrix <- list_all_coex
    
    return(list_res)
  }
    
    list_spp_coex <- 
      lapply(intervals, function(f){
        sp_names_coex <- x2$Species[which(x2$TS >= f & x2$TE <= f)]
        return(sp_names_coex)
      })
    
    # naming time slices
    names(list_spp_coex) <- paste("t", as.character(round(intervals, digits = round.digits)), sep = "_")
    names_interval_long <- 
      as.character(rep(round(intervals, digits = round.digits),
                       unlist(lapply(list_spp_coex, function(f) length(f)))))
    # data frame in long format with species per time interval
    df_long_coex <- data.frame(species = unlist(list_spp_coex), time.interval = names_interval_long)
    
    # just counting the number of species per time slice
    df_summary_coex <- 
      df_long_coex %>% 
      group_by(time.interval) %>% 
      dplyr::mutate(n.spp = dplyr::n()) %>% 
      dplyr::distinct(n.spp)
    
    # dense matrix of species occurrence by each time interval 
    matrix_dense <- as.matrix(phyloregion::long2sparse(x = df_long_coex, grids = "time.interval", species = "species")) 
    
    # coexistence matrix for each time interval
    list_all_coex <- lapply(1:nrow(matrix_dense), function(f) as.matrix(matrix_dense[f, ]) %*% t(matrix_dense[f, ]))
    names(list_all_coex) <- paste("interval", rownames(matrix_dense), sep = "_")
    
    # Organizing results
    list_res <- vector(mode = "list", length = 3)
    names(list_res) <- c("summary.time.occurence", "dense.time.occurrence", "cooccurence.time.matrix")
    list_res$summary.time.occurence <- df_summary_coex
    list_res$dense.time.occurrence <- matrix_dense
    list_res$cooccurence.time.matrix <- list_all_coex
    
    return(list_res)
  }

# testing inputation method for geographical fossil occurrences -----------------------------------------------

df.coords = df_can
df.TS.TE = longs[[1]]
df.coords <- 
  df.coords %>% 
  select(sp.name, lng, lat, midpoint) %>% 
  rename(species = sp.name, lon = lng, lat = lat, age = midpoint)

df.TS.TE <- 
  df.TS.TE %>% 
  mutate(species = rownames(df.TS.TE), TS = TS, TE = TE)

intervals = dog_NALMAs_age
type.fill = "centroid"
range = FALSE
species = "species"
lon = "lon"
lat = "lat"
age = "age"
min.age = "min.age"
max.age = "max.age"


df_coords_fill <- make_gap_filling(df.coords = df.coords, df.TS.TE = df.TS.TE, intervals = intervals, type.fill = "centroid")



# testing reach function --------------------------------------------------

library(dplyr)
df.coords = df_can
df.TS.TE = longs[[1]]
df.coords <- 
  df.coords %>% 
  select(sp.name, lng, lat, midpoint) %>% 
  rename(species = sp.name, lon = lng, lat = lat, age = midpoint)

df.TS.TE <- 
  df.TS.TE %>% 
  mutate(species = rownames(df.TS.TE), TS = TS, TE = TE)

intervals = dog_NALMAs_age

source(here::here("R", "gap-filling.R"))
source(here::here("R", "temporal-coexist.R"))
load(file = "Modified_occurrence_dataset_with_gaps_filled.Rdata")
load("df_can_tot.Rdata")

res_test_reach <- calc_reach(df.coords = df.coords, df.TS.TE = df.TS.TE, intervals = intervals)

lapply(res_test_reach, function(x) dim(x))
