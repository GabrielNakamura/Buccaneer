# Morphospace competition metrics, time_space_multiplication and PyRateContinuous setup  ####

rm(list = ls());gc()

# loading ####
load("./Spatio_temporal_analyses/Time_coex/long_alive.Rdata")
load("./Spatio_temporal_analyses/Time_coex/longevities.RData")
root_ages <- sapply(longs, function (x) floor(max(x$TS) * 10) /10) 

ms_data <- read.table("./final_data.csv", header = T, sep = ",")

# factoring diet_categories
ms_data$diet_cat = factor(ms_data$diet, labels = c("hyper", "meso", "hypo"))

# Matrices prep ####
# Here we define the basis of the following analyses, creating the morphospace distance matrices
# Euclidean distance of both axes #

# mass and carnivory

mass_ld <- data.frame(ms_data$log_mass, ms_data$LD1)
rownames(mass_ld) <- ms_data$species

# distance matrix
mat_ms <- as.matrix(dist(mass_ld, method = "euclidean"))

diag(mat_ms) <- NA

# test if there are species exactly on the same spot of the morphospace:
summary(c(mat_ms[lower.tri(mat_ms)])) # no distance should be equal to 0
sum(duplicated(c(mat_ms[lower.tri(mat_ms)]))) # there should be no equal values

# diet categories matrix

diet <- ms_data$diet_cat
names(diet) <- ms_data$species
diet

# remember, we only want to limit competition to not occur between hyper and hypo 
# carnivorous, but meso can compete with both
mat_diet <- sapply(diet, function(x) sapply(diet, function(y) x == y)) # determining which species share the same diet category

# combing meso with the other plausibe competitors
mat_diet[diet == "hyper", diet == "meso"] <- 1
mat_diet[diet == "meso", diet == "hyper"] <- 1
mat_diet[diet == "hypo", diet == "meso"] <- 1
mat_diet[diet == "meso", diet == "hypo"] <- 1

# now we multiply the distances by the species allowed to compete under this scenario
mat_ms_diet <- mat_ms * mat_diet 

# time_space_multiplications ####
# once again, this section requires a computer with a good amount of RAM

# Reach ####

load("./Spatio_temporal_analyses/Time_space/reach/reach_final.Rdata")

# testing if the species names are consistent
identical(colnames(mat_ms), colnames(mat.reach.nalma[[1]][[1]]))


# Here you can navigate between the preferred method

# MNND ####

# this loop will automatically multiply the morphological datasets by the preferred geographical
# filter and then apply the competition metric and save the different 10 dataset files

mat_ms_ts <- vector("list", 50) # here the "50" is the number of replicas you have, "ts" for "time_space"

for (k in 1:50){ # here the "50" is the number of replicas you have
  mat_ms_ts[[k]] <- lapply(mat.reach.nalma[[k]], function (x) x * mat_ms)
}

# here you can also change the size depending on the points in time scale

nnd_ms_ts <- vector("list", 50)
mnnd_ms_ts <- vector("list", 50) 

for (k in 1:50){
  nnd_ms_ts[[k]] <- vector("list", length(long.alive[[k]]))
  mnnd_ms_ts[[k]] <- vector("list", length(long.alive[[k]]))
}

# MNND loop - some workaround was required when there are 0 or only 1 species in a bin, so as to not produce false 0 distances
for (k in 1:50){ # here the "50" is the number of replicas you have
  for (i in 1:length(long.alive[[k]])){ 
    for (j in 1:length(long.alive[[k]][[i]])) {
      if (rlang::is_empty(long.alive[[k]][[i]])) {
      nnd_ms_ts[[k]][[i]][j] <- NaN
      }
      else if (length(long.alive[[k]][[i]]) == 1) {
      nnd_ms_ts[[k]][[i]][j] <- NaN
      }
      else {
      nnd_ms_ts[[k]][[i]][j] <- min(mat_ms_ts[[k]][[i]][long.alive[[k]][[i]],long.alive[[k]][[i]]][j,][mat_ms_ts[[k]][[i]][long.alive[[k]][[i]],long.alive[[k]][[i]]][j,] > 0], na.rm = T) # here is where we identify for each species its nearest neighbor
      }
    }
  nnd_ms_ts[[k]][[i]][nnd_ms_ts[[k]][[i]] == Inf] <- NA
  mnnd_ms_ts[[k]][[i]] <- mean(unique(nnd_ms_ts[[k]][[i]]), na.rm = T) # and then average those distances in the pool
  mnnd_ms_ts[[k]] <- unlist(mnnd_ms_ts[[k]])
  }
}

# the warning resulting from the loop above does not influence the calculation. It happens because a given species in a given  time interval does not coexist in space with any other species. Hence it is not possible to calculate an morphological distance. in those cases an "inf" is first added which is later replaced by an "NA".

mnnd_ms_table <- vector("list", 50) # here the "50" is the number of replicas you have

for (k in 1:50){ # here the "50" is the number of replicas you have
  mnnd_ms_table[[k]] <- data.frame(time = round(seq(from = root_ages[[k]], to = 0, by = -0.1),1), val = mnnd_ms_ts[[k]])
}


# save the data in .Rdata
dir.create("./Ecomorphological_data/competition/reach/mnnd/", recursive = T)
save(mnnd_ms_table, file = "./Ecomorphological_data/competition/reach/mnnd/reach_mnnd.Rdata", compress = "xz")

# save the data in txt file
for (k in 1:50){ # here the "50" is the number of replicas you have
  aux <- mnnd_ms_table[[k]]
  names(aux) <- c("time", paste0("reach_mnnd"))
  aux <- aux[order(aux$time, decreasing = F),]
  aux[,2][which(is.nan(aux[,2]))] <- 0
  write.table(x = aux, file = paste0("./Ecomorphological_data/competition/reach/mnnd/reach_mnnd", "_", k, ".txt"), row.names = F, quote = F)
}

# MPD ####

# now, for the MPD (and also when we do for site coexistence) is the exact same rationale

mat_ms_diet_ts <- vector("list", 50) # here the "50" is the number of replicas you have

for (k in 1:50){ # here the "50" is the number of replicas you have
  mat_ms_diet_ts[[k]] <- lapply(mat.reach.nalma[[k]], function (x) x * mat_ms_diet)
}

# Applying the methods

mpd_ms_diet_ts <- vector("list", 50) # here the "50" is the number of replicas you have

# you can see the MPD algorithm is a lot simpler because it averages between all possible combinations of competing species, but it adds the diet filter
for (k in 1:50){ # here the "50" is the number of replicas you have
  for (i in 1:length(mat_ms_diet_ts[[k]])){ 
    mpd_ms_diet_ts[[k]][[i]] <- mean(mat_ms_diet_ts[[k]][[i]][lower.tri(mat_ms_diet_ts[[k]][[i]])][mat_ms_diet_ts[[k]][[i]][lower.tri(mat_ms_diet_ts[[k]][[i]])] > 0], na.rm = T)
  }
}

mpd_ms_diet_table <- vector("list", 50) # here the "50" is the number of replicas you have


for (k in 1:50){ # here the "50" is the number of replicas you have
  mpd_ms_diet_table[[k]] <- data.frame(time = round(seq(from = root_ages[[k]], to = 0, by = -0.1),1), val = unlist(mpd_ms_diet_ts[[k]]))
}

dir.create("./Ecomorphological_data/competition/reach/mpd/", recursive = T)
save(mpd_ms_diet_table, file = "./Ecomorphological_data/competition/reach/mpd/reach_mpd_diet.Rdata", compress = "xz")

for (k in 1:50){ # here the "50" is the number of replicas you have
  aux <- mpd_ms_diet_table[[k]]
  names(aux) <- c("time", "reach_mpd")
  aux <- aux[order(aux$time, decreasing = F),]
  aux[,2][which(is.nan(aux[,2]))] <- 0
  write.table(x = aux, file = paste0("./Ecomorphological_data/competition/reach/mpd/reach_mpd_diet", "_", k, ".txt"), row.names = F, quote = F)
}

# site ####

load("./Spatio_temporal_analyses/Time_space/site/site_final.Rdata")

# testing if the species names are consistent
identical(colnames(mat_ms), colnames(mat.sites.nalma[[1]][[1]]))

# Here you can navigate between the preferred method

# MNND ###

# this loop will automatically multiply the morphological datasets by the preferred geographical
# filter and then apply the competition metric and save the different 10 dataset files

mat_ms_ts <- vector("list", 50) # here the "50" is the number of replicas you have, "ts" for "time_space"

for (k in 1:50){ # here the "50" is the number of replicas you have
  mat_ms_ts[[k]] <- lapply(mat.sites.nalma[[k]], function (x) x * mat_ms)
}

nnd_ms_ts <- vector("list", 50)
mnnd_ms_ts <- vector("list", 50) 

for (k in 1:50){
  nnd_ms_ts[[k]] <- vector("list", length(long.alive[[k]]))
  mnnd_ms_ts[[k]] <- vector("list", length(long.alive[[k]]))
}

# MNND loop - some workaround was required when there are 0 or only 1 species in a bin, so as to not produce false 0 distances
for (k in 1:50){ # here the "50" is the number of replicas you have
  for (i in 1:length(long.alive[[k]])){
    for (j in 1:length(long.alive[[k]][[i]])) {
      if (rlang::is_empty(long.alive[[k]][[i]])) {
        nnd_ms_ts[[k]][[i]][j] <- NaN
      }
      else if (length(long.alive[[k]][[i]]) == 1) {
        nnd_ms_ts[[k]][[i]][j] <- NaN
      }
      else {
        nnd_ms_ts[[k]][[i]][j] <- min(mat_ms_ts[[k]][[i]][long.alive[[k]][[i]],long.alive[[k]][[i]]][j,][mat_ms_ts[[k]][[i]][long.alive[[k]][[i]],long.alive[[k]][[i]]][j,] > 0], na.rm = T) # here is where we identify for each species its nearest neighbor
      }
    }
    nnd_ms_ts[[k]][[i]][nnd_ms_ts[[k]][[i]] == Inf] <- NA
    mnnd_ms_ts[[k]][[i]] <- mean(unique(nnd_ms_ts[[k]][[i]]), na.rm = T) # and then average those distances in the pool
    mnnd_ms_ts[[k]] <- unlist(mnnd_ms_ts[[k]])
  }
}

# the warning resulting from the loop above does not influence the calculation. It happens because a given species in a given  time interval does not coexist in space with any other species. Hence it is not possible to calculate an morphological distance. in those cases an "inf" is first added which is later replaced by an "NA" is added.

mnnd_ms_table <- vector("list", 50) # here the "50" is the number of replicas you have

for(k in 1:50){ # here the "50" is the number of replicas you have
  mnnd_ms_table[[k]] <- data.frame(time = round(seq(from = root_ages[[k]], to = 0, by = -0.1), 1), val = mnnd_ms_ts[[k]])
}

# save the data in .Rdata
dir.create("./Ecomorphological_data/competition/site/mnnd/", recursive = T)
save(mnnd_ms_table, file = "./Ecomorphological_data/competition/site/mnnd/site_mnnd.Rdata", compress = "xz")

# save the data in txt file
for (k in 1:50){ # here the "50" is the number of replicas you have
  aux <- mnnd_ms_table[[k]]
  names(aux) <- c("time", paste0("site_mnnd"))
  aux <- aux[order(aux$time, decreasing = F),]
  aux[,2][which(is.nan(aux[,2]))] <- 0
  write.table(x = aux, file = paste0("./Ecomorphological_data/competition/site/mnnd/site_mnnd", "_", k, ".txt"), row.names = F, quote = F)
}

# MPD ####
mat_ms_diet_ts <- vector("list", 50) # here the "50" is the number of replicas you have

for (k in 1:50){ # here the "50" is the number of replicas you have
  mat_ms_diet_ts[[k]] <- lapply(mat.sites.nalma[[k]], function (x) x * mat_ms_diet)
}

# Applying the methods

mpd_ms_diet_ts <- vector("list", 50) # here the "50" is the number of replicas you have

# you can see the MPD algorithm is a lot simpler because it averages between all possible combinations of competing species, but it adds the diet filter
for (k in 1:50){ # here the "50" is the number of replicas you have
  for (i in 1:length(long.alive[[k]])){
    mpd_ms_diet_ts[[k]][[i]] <- mean(mat_ms_diet_ts[[k]][[i]][lower.tri(mat_ms_diet_ts[[k]][[i]])][mat_ms_diet_ts[[k]][[i]][lower.tri(mat_ms_diet_ts[[k]][[i]])] > 0], na.rm = T)
  }
}

mpd_ms_diet_table <- vector("list", 50) # here the "50" is the number of replicas you have

for (k in 1:50){ # here the "50" is the number of replicas you have
  mpd_ms_diet_table[[k]] <- data.frame(time = round(seq(from = root_ages[[k]], to = 0, by = -0.1),1), val = unlist(mpd_ms_diet_ts[[k]]))
}

dir.create("./Ecomorphological_data/competition/site/mpd/", recursive = T)
save(mpd_ms_diet_table, file = "./Ecomorphological_data/competition/site/mpd/site_mpd_diet.Rdata", compress = "xz")

for (k in 1:50){ # here the "50" is the number of replicas you have
  aux <- mpd_ms_diet_table[[k]]
  names(aux) <- c("time", "site_mpd")
  aux <- aux[order(aux$time, decreasing = F),]
  aux[,2][which(is.nan(aux[,2]))] <- 0
  write.table(x = aux, file = paste0("./Ecomorphological_data/competition/site/mpd/site_mpd_diet", "_", k, ".txt"), row.names = F, quote = F)
}

# regional ####

load("./Spatio_temporal_analyses/Time_coex/mat_time_replicas.Rdata")

# testing if the species names are consistent
identical(colnames(mat_ms), colnames(mat.time.replicas[[1]][[1]]))

# Here you can navigate between the preferred method

# MNND ####

# this loop will automatically multiply the morphological datasets by the preferred geographical
# filter and then apply the competition metric and save the different 10 dataset files

mat_ms_ts <- vector("list", 50) # here the "50" is the number of replicas you have, "ts" for "time_space"


for (k in 1:50){ # here the "50" is the number of replicas you have
  mat_ms_ts[[k]] <- lapply(mat.time.replicas[[k]], function (x) x * mat_ms)
}

# here you can also change the size depending on thepoints in time scale

nnd_ms_ts <- vector("list", 50)
mnnd_ms_ts <- vector("list", 50) 

for (k in 1:50){
  nnd_ms_ts[[k]] <- vector("list", length(long.alive[[k]]))
  mnnd_ms_ts[[k]] <- vector("list", length(long.alive[[k]]))
}


# MNND loop - some workaround was required when there are 0 or only 1 species in a bin, so as to not produce false 0 distances
for (k in 1:50){ # here the "50" is the number of replicas you have
  for (i in 1:length(long.alive[[k]])){ 
    for (j in 1:length(long.alive[[k]][[i]])) {
      if (rlang::is_empty(long.alive[[k]][[i]])) {
        nnd_ms_ts[[k]][[i]][j] <- NaN
      }
      else if (length(long.alive[[k]][[i]]) == 1) {
        nnd_ms_ts[[k]][[i]][j] <- NaN
      }
      else {
        nnd_ms_ts[[k]][[i]][j] <- min(mat_ms_ts[[k]][[i]][long.alive[[k]][[i]],long.alive[[k]][[i]]][j,][mat_ms_ts[[k]][[i]][long.alive[[k]][[i]],long.alive[[k]][[i]]][j,] > 0], na.rm = T) # here is where we identify for each species its nearest neighbor
      }
    }
    nnd_ms_ts[[k]][[i]][nnd_ms_ts[[k]][[i]] == Inf] <- NA
    mnnd_ms_ts[[k]][[i]] <- mean(unique(nnd_ms_ts[[k]][[i]]), na.rm = T) # and then average those distances in the pool
    mnnd_ms_ts[[k]] <- unlist(mnnd_ms_ts[[k]])
  }
}

# the warning resulting from the loop above does not influence the calculation. It happens because a given species in a given  time interval does not coexist in space with any other species. Hence it is not possible to calculate an morphological distance. in those cases an "inf"is first added which is later replaced by an "NA" is added.

mnnd_ms_table <- vector("list", 50) # here the "50" is the number of replicas you have

for (k in 1:50){ # here the "50" is the number of replicas you have
  mnnd_ms_table[[k]] <- data.frame(time = round(seq(from = root_ages[[k]], to = 0, by = -0.1),1), val = mnnd_ms_ts[[k]])
}

# save the data in .Rdata
dir.create("./Ecomorphological_data/competition/regional/mnnd", recursive = T)
save(mnnd_ms_table, file = "./Ecomorphological_data/competition/regional/mnnd/regional_mnnd.Rdata", compress = "xz")

# save the data in txt file
for (k in 1:50){ # here the "50" is the number of replicas you have
  aux <- mnnd_ms_table[[k]]
  names(aux) <- c("time", paste0("regional_mnnd"))
  aux <- aux[order(aux$time, decreasing = F),]
  aux[,2][which(is.nan(aux[,2]))] <- 0
  write.table(x = aux, file = paste0("./Ecomorphological_data/competition/regional/mnnd/regional_mnnd", "_", k, ".txt"), row.names = F, quote = F)
}

# MPD ####

# now, for the MPD (and also when we do for site coexistence) is the exact same rationale

mat_ms_diet_ts <- vector("list", 50) # here the "50" is the number of replicas you have

for (k in 1:50){ # here the "50" is the number of replicas you have
  mat_ms_diet_ts[[k]] <- lapply(mat.time.replicas[[k]], function (x) x * mat_ms_diet)
}

# Applying the methods

mpd_ms_diet_ts <- vector("list", 50) # here the "50" is the number of replicas you have

# you can see the MPD algorithm is a lot simpler because it averages between all possible combinations of competing species, but it adds the diet filter
for (k in 1:50){ # here the "50" is the number of replicas you have
  for (i in 1:length(mat_ms_diet_ts[[k]])){ 
    mpd_ms_diet_ts[[k]][[i]] <- mean(mat_ms_diet_ts[[k]][[i]][lower.tri(mat_ms_diet_ts[[k]][[i]])][mat_ms_diet_ts[[k]][[i]][lower.tri(mat_ms_diet_ts[[k]][[i]])] > 0], na.rm = T)
  }
}

mpd_ms_diet_table <- vector("list", 50) # here the "50" is the number of replicas you have

for (k in 1:50){ # here the "50" is the number of replicas you have
  mpd_ms_diet_table[[k]] <- data.frame(time = round(seq(from = root_ages[[k]], to = 0, by = -0.1),1), val = unlist(mpd_ms_diet_ts[[k]]))
}

dir.create("./Ecomorphological_data/competition/regional/mpd", recursive = T)
save(mpd_ms_diet_table, file = "./Ecomorphological_data/competition/regional/mpd/regional_mpd_diet.Rdata", compress = "xz")

for (k in 1:50){ # here the "50" is the number of replicas you have
  aux <- mpd_ms_diet_table[[k]]
  names(aux) <- c("time", "regional_mpd")
  aux <- aux[order(aux$time, decreasing = F),]
  aux[,2][which(is.nan(aux[,2]))] <- 0
  write.table(x = aux, file = paste0("./Ecomorphological_data/competition/regional/mpd/regional_mpd_diet", "_", k, ".txt"), row.names = F, quote = F)
}


##### Data is ready to input in PyRateContinuous analyses, refer to "16_PyRate_Commands.txt"

# It is advised to create a final directory in Continuous and paste all the time series .txt files, from diversity curves to the files created here
dir.create("./Continuous/time_series/")