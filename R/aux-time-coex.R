#' Auxiliar function to compute matrix of temporal coexistence
#'
#' @param data_alive A list with a set of data frames indicating the species present at a given time interval
#'
#' @return
#' @export
#'
#' @examples
aux.compute.time.coex <- 
  function(data_alive){
  for(k in seq_along(data_alive)){
    print(k)
    mat_time_all <- 
      lapply(data_alive[[k]], function(x){
      mat_time[,] <- 0
      mat_time[x, x] <- 1
      return(mat_time)
    }
    )
    mat.time.replicas[[k]] <- mat_time_all
  }
}
