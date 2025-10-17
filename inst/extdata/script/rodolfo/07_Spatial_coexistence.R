# Spatial analyses of 50 replicas ####

#From the result of the previous script (gap_filling), apply the analyzes that will
#determine the coexistence between species according to the two methods (reach and site).
#The script will take the complete data, apply the distance calculations, determine coexistence
#and save it as matrices

rm(list = ls());gc()

library(geosphere)
library(plyr)
library(tidyverse)

# load("./Spatio_temporal_analyses/Time_coex/longevities.RData")
load(here::here("inst", "extdata", "script", "rodolfo", "longevities.RData"))
# load("./PBDB/df_can.Rdata")
load(here::here("inst", "extdata", "script", "rodolfo", "df_can.Rdata"))
options(scipen = 999) # this will be needed in a step that converts distances

NALMA <- c("Duchesnean", "Chadronian","Orellan","Whitneyan","Arikareean","Hemingfordian","Barstovian","Clarendonian","Hemphillian","Blancan","Irvingtonian","Rancholabrean_Present") # Barnsosky 2014
NALMA_age <- c(39.7, 37, 33.9, 31.8, 29.5, 18.5, 16.3, 12.5, 9.4, 4.7, 1.4, 0.21,0) # Barnosky 2014

# Reaching distance ####
# Setup and analyes

#setwd("./Spatio_temporal_analyses/Spatial_coex/")

# load the dataset with the new points added from the previous script
# load("./Spatio_temporal_analyses/Spatial_coex/df_can_tot.Rdata")
load(here::here("inst", "extdata", "script", "rodolfo", "df_can_tot.Rdata"))

# first, split the datasets by each species occurrence in each nalma
sp_split_nalmas_tot <- lapply(df_can_tot, function (x) split(x, list(x$species, x$NALMA), drop = T))

# creates the matrix of distances between each pair of occurrences of each NALMA, shortest distance on an ellipsoid WGS84
list.distance.matrices.nalmas.tot <- lapply(sp_split_nalmas_tot, function(x) lapply(x, function(df) distm(cbind(df$lng, df$lat),fun = distGeo)/1000))

# vectorizing
list.distance.vectors.nalmas.tot <- lapply(list.distance.matrices.nalmas.tot, function(x) lapply(x, function(x){
  a <- c(x)
  b <- unique(a)
  b
}
))

# treating the "self distance", the distance between a point and itself as NA instead
# of 0 to not interfere with following calculations
for (k in 1:50){ # here the "50" should be the number of replicas you have
  for (i in 1:length(list.distance.vectors.nalmas.tot[[k]])){
    list.distance.vectors.nalmas.tot[[k]][[i]][1] <- NA
  }
}

# applying functions over the list to recover distribution values
# you are probably not going to use all of them, but I'll keep it written down

sum.distances.nalmas.tot <- sapply(list.distance.vectors.nalmas.tot, function (x) sapply(x, function(x) summary(x)))
min.distances.nalmas.tot <- sapply(list.distance.vectors.nalmas.tot, function (x) sapply(x, function (x) min(x, na.rm = T)))
### the warning produced when running line 57 is not really a problem, it relates to "infinity distances" that are fixed below.

# replacing Inf values again to not interfere with calculations
for (k in 1:50){# here the "50" should be the number of replicas you have
  min.distances.nalmas.tot[[k]][min.distances.nalmas.tot[[k]] == Inf] <- NA
}
mean.distances.nalmas.tot <- sapply(list.distance.vectors.nalmas.tot, function (x) sapply(x, function (x) mean(x, na.rm = T)))
# same for NaN
for (k in 1:50){# here the "50" should be the number of replicas you have
  mean.distances.nalmas.tot[[k]][mean.distances.nalmas.tot[[k]] == "NaN"] <- NA
}
median.distances.nalmas.tot <- lapply(list.distance.vectors.nalmas.tot, function (x) sapply(x, function(x) median(x, na.rm = T)))
max.distances.nalmas.tot <- lapply(list.distance.vectors.nalmas.tot, function (x) sapply(x, function(x) max(x, na.rm = T)))
### the warning produced when running line 70 is not really a problem, it relates to "infinity distances" that are fixed below.

# and for -Inf
for (k in 1:50){# here the "50" should be the number of replicas you have
  max.distances.nalmas.tot[[k]][max.distances.nalmas.tot[[k]] == -Inf] <- NA
}

# grouping the main ones

distances.df.nalmas.tot <- vector("list", 50)# here the "50" should be the number of replicas you have
for (k in 1:50){# here the "50" should be the number of replicas you have
  distances.df.nalmas.tot[[k]] <- data.frame(min.distances.nalmas.tot[[k]], mean.distances.nalmas.tot[[k]], median.distances.nalmas.tot[[k]], max.distances.nalmas.tot[[k]])
}

# Estimating the pairwise distances ####

# extract lat and long values of the distributions
sp_split_coords_nalmas_tot <- lapply(sp_split_nalmas_tot, function (x) lapply(x, function(x) cbind(x$lng, x$lat)))

# applies distm function between all points of species in a pairwise fashion,
# then retrieves the minimum distance between them

# empty vector
pwdist.min.nalmas.tot <- vector("list", 50)# here the "50" should be the number of replicas you have
# warning - can take up to 1 hour to run
for (k in 1:50){# here the "50" should be the number of replicas you have
  #print(k)
  pwdist.min.nalmas.tot[[k]] <- sapply(sp_split_coords_nalmas_tot[[k]], function (x) sapply(sp_split_coords_nalmas_tot[[k]], function (y) min(distm(x,y, fun = distGeo)/1000)))
}

save(pwdist.min.nalmas.tot, file = "./Spatio_temporal_analyses/Spatial_coex/pwdist.Rdata", compress = "xz")

# in the end you'll have 50 (depends on how many replicas you have) matrices of the min distance
# between each pair of species in each NALMA

# REACH - distance threshold ####

# There is coexistence when the sum of the max.distances of two species is greater the minimal distance between them

# Make a matrix of pairwise sum of distances
# by NALMAS ####

# convert NA back to 0
for (k in 1:50){# here the "50" should be the number of replicas you have
  max.distances.nalmas.tot[[k]][is.na(max.distances.nalmas.tot[[k]])] <- 0
}

# To infer the dispersal potential of species, we estimated the max value of distance ever
# reached of a species in a NALMA

max.df.tot <- lapply(max.distances.nalmas.tot, function (x) data.frame(x, row.names = NULL, sp = names(x)))
for (k in 1:50){# here the "50" should be the number of replicas you have
  max.df.tot[[k]] <- max.df.tot[[k]][,c(2,1)]
  names(max.df.tot[[k]]) <- c("sp", "max.distances.nalmas")
}

max.df.tot <- lapply(max.df.tot, function (x) separate(x, 1, c("sp", "NALMA"), sep = "\\."))

# ab comes from absolute
ab.max.tot <- lapply(max.df.tot, function (x) tapply(X = x$max.distances.nalmas, INDEX = x$sp, FUN = max))

ab.max.tot <- lapply(ab.max.tot, function (x) x[order(factor(names(x), levels = spp))])

# now comes some workaround, I used the new vector to replace the oldest value,
# estimated from each nalma
dist.temp1.tot <- max.distances.nalmas.tot

# that lead to a conflict in naming
for (k in 1:50){# here the "50" should be the number of replicas you have
  names(dist.temp1.tot[[k]]) <- gsub(pattern = "\\..*", replacement = "", x = names(max.distances.nalmas.tot[[k]]))
}

# some reordering required
dist.temp2.tot <- lapply(dist.temp1.tot, function (x) table(names(x)))

dist.temp2.tot <- lapply(dist.temp2.tot, function(x) x[order(factor(names(x), levels = spp))])

# combine all of it and rename
# empty vectors
max.distance.ab.tot <- vector("list", 50)# here the "50" should be the number of replicas you have
for (k in 1:50){# here the "50" should be the number of replicas you have
  max.distance.ab.tot[[k]] <- rep(ab.max.tot[[k]], times = as.vector(dist.temp2.tot[[k]]))
}

for (k in 1:50){# here the "50" should be the number of replicas you have
  max.distance.ab.tot[[k]] <- max.distance.ab.tot[[k]][match(names(dist.temp1.tot[[k]]), names(max.distance.ab.tot[[k]]))]
}

for (k in 1:50){# here the "50" should be the number of replicas you have
  names(max.distance.ab.tot[[k]]) <- names(max.distances.nalmas.tot[[k]])
}

pwdist.sum.ab.tot <- vector("list", 50)# here the "50" should be the number of replicas you have
for (k in 1:50){# here the "50" should be the number of replicas you have
  #print(k)
  pwdist.sum.ab.tot[[k]] <- sapply(max.distance.ab.tot[[k]], function (x) sapply(max.distance.ab.tot[[k]], function (y) sum(x,y)))
}

# now, from the matrices of distances, estimate if species coexisted according to the method
mat.reach.ab.tot <- vector("list", 50)# here the "50" should be the number of replicas you have
for (k in 1:50){# here the "50" should be the number of replicas you have
  mat.reach.ab.tot[[k]] <- pwdist.sum.ab.tot[[k]] >= pwdist.min.nalmas.tot[[k]]
  mat.reach.ab.tot[[k]][which(mat.reach.ab.tot[[k]] == TRUE)] <- 1
}

mat.times.tot <- vector("list", 50)# here the "50" should be the number of replicas you have

for (k in 1:50){# here the "50" should be the number of replicas you have
  #print(k)
  for (i in 1:length(NALMA)){
  mat.times.tot[[k]][[i]] <- matrix(nrow = 50, ncol = 50) # both "50" are the number of replicas
  }
}

for (k in 1:50){# here the "50" should be the number of replicas you have
  #print(k)
  for (i in 1:length(NALMA)){
    #print(i)
    mat.times.tot[[k]][[i]] <- mat.reach.ab.tot[[k]][grep(pattern = NALMA[i], x = rownames(mat.reach.ab.tot[[k]])), grep(pattern = NALMA[i], x = colnames(mat.reach.ab.tot[[k]])), drop = F]
  }
}

for (k in 1:50){# here the "50" should be the number of replicas you have
  for (i in 1:length(mat.times.tot[[k]])){
    rownames(mat.times.tot[[k]][[i]]) <- gsub(pattern = "\\..*", replacement = "", x = rownames(mat.times.tot[[k]][[i]]))
    colnames(mat.times.tot[[k]][[i]]) <- gsub(pattern = "\\..*", replacement = "", x = rownames(mat.times.tot[[k]][[i]]))
  }
}

# rewrite the matrices according to our pattern
for (k in 1:50){# here the "50" should be the number of replicas you have
  names(mat.times.tot[[k]]) <- NALMA
}

mat.temp <- matrix(0, ncol = length(spp), nrow = length(spp))
rownames(mat.temp) <- spp # names the rows with sps
colnames(mat.temp) <- spp

mat.times.final.tot <- vector("list", 50)# here the "50" should be the number of replicas you have

for (k in 1:length(mat.times.final.tot)){
  mat.times.final.tot[[k]] <- rep(list(mat.temp), times = length(NALMA))
}

for (k in 1:50){# here the "50" should be the number of replicas you have
  for (i in 1:length(mat.times.final.tot[[k]])){
    mat.times.final.tot[[k]][[i]][rownames(mat.times.tot[[k]][[i]]), colnames(mat.times.tot[[k]][[i]])] <- mat.times.tot[[k]][[i]]
  }
}

for (k in 1:50){# here the "50" should be the number of replicas you have
  names(mat.times.final.tot[[k]]) <- NALMA
}

mat.times.ab.tot <- mat.times.final.tot
rm(mat.times.final.tot)

# ok, after all that recombining and renaming, you have the final matrices of coexistence
# and you can save the .Rdata

dir.create("./Spatio_temporal_analyses/Spatial_coex/reach/", recursive = T)
save(mat.times.ab.tot, file = "./Spatio_temporal_analyses/Spatial_coex/reach/mat_reach_tot.Rdata", compress = "xz")

# sites ####

# this one is actually easier to implement and to understand

sites <- unique(df_can$site)

df_can$sp.nalma <- paste(df_can$species, df_can$NALMA, sep = ".")
sp.nalma <- names(split(df_can, list(df_can$species, df_can$NALMA), T))

alive.in.site <- list(NA)
for (k in 1:length(sites)){
  alive.in.site[[k]] <- df_can[which(df_can$site == sites[k]),"sp.nalma"]
} # list with species names in coexisting sites

names(alive.in.site) <- sites

mat <- matrix(0, nrow = length(sp.nalma), ncol = length(sp.nalma))
# mat <- data.frame(mat)
rownames(mat) <- sp.nalma
colnames(mat) <- sp.nalma

mat.site <- mat

for(j in 1:length(sites)){
  mat.site[alive.in.site[[j]], alive.in.site[[j]]] <- 1
}

mat.times <- list(NA)
for (i in 1:length(NALMA)){
  mat.times[[i]] <- mat.site[grep(pattern = NALMA[i], x = rownames(mat.site)),grep(pattern = NALMA[i], x = colnames(mat.site)), drop = F]
}

for (i in 1:length(mat.times)){
  rownames(mat.times[[i]]) <- gsub(pattern = "\\..*", replacement = "", x = rownames(mat.times[[i]]))
  colnames(mat.times[[i]]) <- gsub(pattern = "\\..*", replacement = "", x = colnames(mat.times[[i]]))
}

names(mat.times) <- NALMA

mat.temp <- matrix(0, ncol = length(spp), nrow = length(spp))
rownames(mat.temp) <- spp
colnames(mat.temp) <- spp

mat.times.final <- rep(list(mat.temp), times = length(mat.times))
for (k in 1:length(mat.times.final)){
  mat.times.final[[k]][rownames(mat.times[[k]]), colnames(mat.times[[k]])] <- mat.times[[k]]
}

names(mat.times.final) <- NALMA
mat.site.tot <- mat.times.final

# you can also save that object
dir.create("./Spatio_temporal_analyses/Spatial_coex/site/", recursive = T)
save(mat.site.tot, file = here::here("inst", "extdata", "script", "rodolfo", "mat_site_tot.Rdata"), compress = "xz")
