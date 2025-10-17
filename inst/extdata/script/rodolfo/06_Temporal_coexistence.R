# analyses of 50 replicas ####

rm(list = ls());gc()
library(gtools)

load("./Spatio_temporal_analyses/Time_coex/longevities.RData")

spp

# determine the duration of the bins and the scale
# in order to avoid overestimating species true geographical data, used an approximation of root age for each replica. the formula sets the root age as the next interval in 0.1 myr scale from the speciation time of the first species
root_ages <- sapply(longs, function (x) floor(max(x$TS) * 10) /10) #

for (k in 1:length(longs)){
  longs[[k]]$TS[which.max(longs[[k]]$TS)] <- root_ages[[k]]
}

# create empty lists to store the information of which species are present in each bin
# long.alive = repeat the process for each replica
# alive.in.bin = which species are present in each 0.1 scale
long.alive <- list(NA)
alive.in.bin <- list(NA)
intervals <- list(NA)

for (k in 1:length(longs)){
  intervals[[k]] <- round(seq(from = root_ages[[k]], to = 0, by = -0.1),1)
  alive.in.bin[[k]] <- vector(mode = "list", length = length(seq(from = root_ages[[k]], to = 0, by = -0.1)))
  names(alive.in.bin[[k]]) <- intervals[[k]]
}

for (k in 1:length(longs)){
  for (i in 1:length(intervals[[k]])){
    alive.in.bin[[k]][[i]] <- rownames(longs[[k]][which(longs[[k]]$TS >= intervals[[k]][i] & longs[[k]]$TE <= intervals[[k]][i]),])
    long.alive[[k]] <- alive.in.bin[[k]]
  }
}

# creates the matrix of coexistence
# the numbers of cols and rows manually depend of the numbers of species
# starting with the base matrix template
mat <- matrix(0, ncol = length(spp), nrow = length(spp))
# each row and column named after the species
rownames(mat) <- spp
colnames(mat) <- spp

# empty repeated matrices for the different replicas
mat.time <- mat
mat.time.all <- vector(mode = "list", length = length(longs))
mat.time.replicas <- vector(mode = "list", length = length(longs))
for (k in 1:length(longs)){
  mat.time.all[[k]] <- rep(list(mat.time), length(intervals[[k]]))
  names(mat.time.all[[k]]) <- intervals[[k]]
  mat.time.replicas[[k]] <- vector(mode = "list", length = length(intervals[[k]]))
}

# from the previous vectors, fills the empty matrices with coexistence data, for each replica
for(k in 1:length(longs)){
  # print(k)
  for(i in 1:length(mat.time.all[[k]])){
    mat.time.all[[k]][[i]][,] <- 0
    mat.time.all[[k]][[i]][long.alive[[k]][[i]], long.alive[[k]][[i]]] <- 1
  }
  mat.time.replicas[[k]] <- mat.time.all[[k]]
}

# you can save the .Rdata with the objects
# save(mat.time.replicas, file = "./Spatio_temporal_analyses/Time_coex/mat_time_replicas.Rdata", compress = "xz")
save(mat.time.replicas, file = here::here("inst", "extdata", "script", "rodolfo", "mat_time_replicas.Rdata"), compress = "xz")
# save(long.alive, file = "./Spatio_temporal_analyses/Time_coex/long_alive.Rdata", compress = "xz")
save(long.alive, file = here::here("inst", "extdata", "script", "rodolfo", "long_alive.Rdata"), compress = "xz")

# you can extract diversity curves from the vectors to use in different diversity analyses

div.curves <- list(NA)
for (k in 1:length(long.alive)){
  div.curves[[k]]  <- sapply(long.alive[[k]], length)
}

# write div.curves as data for continuous analyses

div.curves.frame <- lapply(div.curves, FUN = function (x) data.frame(time = as.numeric(names(x)), div = x, row.names = NULL))
div.curves.frame <- lapply(div.curves.frame, function (x) x[order(x$time),])

# Create diversity curves directory
dir.create("./Spatio_temporal_analyses/Diversity_curves/", recursive = T)

for (k in 1:length(div.curves)){
  write.table(x = div.curves.frame[[k]], file = paste0("./Spatio_temporal_analyses/Diversity_curves/diversity", k, ".txt"), row.names = F, quote = F)
}
# save(div.curves, file = "./Spatio_temporal_analyses/Diversity_curves/div_curves.Rdata", compress = "xz")
save(div.curves, file = here::here("inst", "extdata", "script", "rodolfo", "div_curves.Rdata"), compress = "xz")
