# Compiling results from PyRateContinuous after steps 05 from "PyRate_Commands.txt" ####
rm(list = ls());gc()

library(gtools)
library(plyr)

# load all the data ####
replica <- c(1:50)

# Here you can customize a different subset of time series to read and combine in a single object. Here we report the ones used in the main text.
vars <- c("div", "reach_diversity", "site_diversity", "regional_mpd", "reach_mpd", "site_mpd", "regional_mnnd", "reach_mnnd", "site_mnnd", "temperature")

# Reads data from your PyRateContinuous results directory.
# When setting different time windows in a PyRateContinuous analysis, the output .log files are named accordingly. Here you can customize the file names according to the time windows you provided.
# Here we provide the scales used in the article, "epoch" time scale for competition time series and temperature, and also an analyses through the "cenozoic" as a whole for temperature.
# It is also advised to keep all log files resulting from the different analyses in the same folder.

# Time series and temperature for epoch scale
epoch_results <- vector(mode = 'list', length(vars))
for (k in 1:length(vars)){
  message(paste("Reading variable:",vars[[k]]))
  files <- c(paste0("./Continuous/logs/Canidae_N_america_",replica,"rep_",replica,"_Grj_se_est_",vars[k],"_0_s_34.0_23.0_5.0_2.0linSp_linEx_HP.log"))
    epoch_results[[k]] <- do.call(rbind.fill, lapply(files, read.table, header = T))
    epoch_results[[k]] <- epoch_results[[k]][-which(epoch_results[[k]]$it %in% 0:2999000),]
    epoch_results[[k]] <- epoch_results[[k]][c(1,grep(pattern = "Gl", colnames(epoch_results[[k]])), grep(pattern = "Gm", colnames(epoch_results[[k]])), grep(pattern = "l0", colnames(epoch_results[[k]])) # reads the Gl and Gm, correlation parameters
, grep(pattern = "m0", colnames(epoch_results[[k]])))] # read l0 and m0, baseline rates
    # rename the windows accordingly
    epoch_results[[k]]$Gl_t38.34 <- dplyr::coalesce(epoch_results[[k]]$Gl_t38.34, epoch_results[[k]]$Gl_t37.34)
    epoch_results[[k]]$Gm_t38.34 <- dplyr::coalesce(epoch_results[[k]]$Gm_t38.34, epoch_results[[k]]$Gm_t37.34)
    epoch_results[[k]]$l0_t38.34 <- dplyr::coalesce(epoch_results[[k]]$l0_t38.34, epoch_results[[k]]$l0_t37.34)
    epoch_results[[k]]$m0_t38.34 <- dplyr::coalesce(epoch_results[[k]]$m0_t38.34, epoch_results[[k]]$m0_t37.34)
    epoch_results[[k]] <- epoch_results[[k]][,-grep(pattern = "37.34", x = colnames(epoch_results[[k]]))]
    epoch_results[[k]] <- epoch_results[[k]][,c(1,6,2,3,4,5,11,7,8,9,10,16,12,13,14,15,21,17,18,19,20)]
    epoch_results[[k]]$replica <- rep(1:50, each = 2700)
}
names(epoch_results) <- vars

    
# Temperature for the cenozoic    
files <- c(paste0("./Continuous/logs/Canidae_N_america_",replica,"rep_",replica,"_Grj_se_est_temperature_0_linSp_linEx_HP.log"))
ceno_results <- vector(mode = 'list', 1)
  ceno_results <- do.call(rbind.fill, lapply(files, read.table, header = T))
  ceno_results<- ceno_results[-which(ceno_results$it %in% 0:2999000),]
  ceno_results <- ceno_results[c(1,grep(pattern = "Gl", colnames(ceno_results)), grep(pattern = "Gm", colnames(ceno_results)), grep(pattern = "l0", colnames(ceno_results)),grep(pattern = "m0", colnames(ceno_results)))]
  # single time window: 38 myr to 0
  ceno_results$Gl_t38.0 <- dplyr::coalesce(ceno_results$Gl_t38.0, ceno_results$Gl_t37.0)
  ceno_results$Gm_t38.0 <- dplyr::coalesce(ceno_results$Gm_t38.0, ceno_results$Gm_t37.0)
  ceno_results$l0_t38.0 <- dplyr::coalesce(ceno_results$l0_t38.0, ceno_results$l0_t37.0)
  ceno_results$m0_t38.0 <- dplyr::coalesce(ceno_results$m0_t38.0, ceno_results$m0_t37.0)
  ceno_results <- ceno_results[,-grep(pattern = "37.0", x = colnames(ceno_results))]
  ceno_results$replica <- rep(1:50, each = 2700)


# After reading the files, adds the function to coalesce to ggplot format
# epoch scale
  epoch_results_tables <- vector("list", length(epoch_results))
for (k in 1:length(vars)){
    epoch_results_tables[[k]][[1]] <- stack(epoch_results[[k]][,grep(pattern = "Gl", colnames(epoch_results[[k]]))])
    epoch_results_tables[[k]][[2]] <- stack(epoch_results[[k]][,grep(pattern = "Gm", colnames(epoch_results[[k]]))])
    epoch_results_tables[[k]][[3]] <- stack(epoch_results[[k]][,grep(pattern = "l0", colnames(epoch_results[[k]]))])
    epoch_results_tables[[k]][[4]] <- stack(epoch_results[[k]][,grep(pattern = "m0", colnames(epoch_results[[k]]))])
    names(epoch_results_tables[[k]]) <- c("Gl", "Gm", "l0", "m0")
}
    
# Cenozoic temperature
ceno_results_table <- stack(ceno_results)

# Save the result_table format in .RData, rename accordingly to your time windows
save(epoch_results_tables, file = "./Continuous/logs/epoch_results_tables.RData", compress = "xz")

save(ceno_results_table, file = "./Continuous/logs/ceno_results_table.RData", compress = "xz")