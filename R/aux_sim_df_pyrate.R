#' Auxiliary function to simulate data in a pyrate format
#'
#' @param n0 
#' @param lambda 
#' @param mu 
#' @param tMax 
#' @param seed 
#' @param rho 
#' @param bins 
#' @param returnAll 
#'
#' @return
#' @export
#'
#' @examples

library(dplyr)
sim.df.pyrate <- 
  function(n0,
           lambda,
           mu,
           tMax, 
           seed,
           rho,
           bins,
           returnAll = TRUE
  ){
    set.seed(seed)
    # simulating TS and TE
    sim <- paleobuddy::bd.sim(n0, lambda, mu, tMax)
    
    # simulating fossil samples
    fossils <- paleobuddy::sample.clade(sim = sim, rho = rho, tMax = tMax, bins = bins, returnAll = T)
    
    # naming species and binding in a data frame
    names(sim$TS) <- paste("t", 1:length(sim$TS), sep = "")
    names(sim$TE) <- paste("t", 1:length(sim$TE), sep = "")
    df_pyrate <- data.frame(TS = sim$TS, TE = sim$TE, species = names(sim$TE))
    
    # joining fossil occurrence and TS and TE data
    df_pyrate2 <- 
      fossils %>% 
      dplyr::left_join(df_pyrate, by = c(Species = "species"))
    
    return(df_pyrate2) # df output
    
  }