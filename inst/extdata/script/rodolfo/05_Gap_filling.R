# Dealing with species with no geographical data for time spans which PyRate suggests are present ####

rm(list = ls()); gc()

library(sf)
library(plyr)

load("./PBDB/df_can.Rdata")
load("./Spatio_temporal_analyses/Time_coex/longevities.RData")

# in order to avoid overestimating species true geographical data, used an approximation of root age for each replica. the formula sets the root age as the next interval in 0.1 myr scale from the speciation time of the first species
root_ages <- sapply(longs, function (x) floor(max(x$TS) * 10) /10) # 

for (k in 1:length(longs)){
  longs[[k]]$TS[which.max(longs[[k]]$TS)] <- root_ages[[k]]
}

# working with spatial data - using the SPDF methods
spdf <- st_as_sf(df_can, coords = c("lng", "lat"), crs = 4326) # most common coordinate system == Geodetic CRS: WGS84

# define the name of the intervals 
levels(spdf$NALMA) <- c("Duchesnean", "Chadronian","Orellan","Whitneyan","Arikareean","Hemingfordian","Barstovian","Clarendonian","Hemphillian","Blancan","Irvingtonian","Rancholabrean_Present")
spdf$NALMA <- factor(spdf$NALMA, levels = levels(spdf$NALMA))
spdf$species <- factor(spdf$species, levels = spp)
NALMAs <- factor(levels(spdf$NALMA), levels = levels(spdf$NALMA))

# splitting dataset according to species and NALMAS

# Polygons by NALMAs ####
sp_split_nalmas <- split(spdf,f = list(spdf$species,spdf$NALMA),drop = T)

# first, find out where there is missing data
# identifying species longevities and sampling

# longs are the real longevities estimated from PyRate
# df_can holds the original longevities

# applying the functions to determine the gaps ####

# determining if species longevities cross the gaps of the NALMAS
intervals <- c(39.7, 37, 33.9, 31.8, 29.5, 18.5, 16.3, 12.5, 9.4, 4.7, 1.4, 0.21,0)
div.tot <- vector(mode = "list", length = 50) # here the 50 is the number of replicas you have
div <- list(NA)
for (k in 1:50){# here the 50 is the number of replicas you have
  for (i in 1:(length(intervals)-1)){
    FL <- (longs[[k]]$TS <= intervals[i]) & (longs[[k]]$TS >= intervals[i+1]) & (longs[[k]]$TE <= intervals[i]) & (longs[[k]]$TE >= intervals[i+1])
    bL <- (longs[[k]]$TS > intervals[i]) & (longs[[k]]$TE < intervals[i]) & (longs[[k]]$TE > intervals[i+1])
    Ft <- (longs[[k]]$TS < intervals[i]) & (longs[[k]]$TS > intervals[i+1]) & (longs[[k]]$TE < intervals[i+1])
    bt <- (longs[[k]]$TS >= intervals[i]) & (longs[[k]]$TE <= intervals[i+1])
    div[[i]] <- FL+bL+Ft+bt
    div.tot[[k]] <- div
  }
}

# estimate which species are present in each NALMA
alive.in.nalma.tot <- vector(mode = "list", length = 50) # here the 50 is the number of replicas you have
alive.in.nalma <- list(NA)
for (k in 1:50){ # here the 50 is the number of replicas you have
  for (i in 1:length(NALMAs)){
    alive.in.nalma[[i]] <- spp[c(which(div.tot[[k]][[i]] == 1))]
    alive.in.nalma.tot[[k]] <- alive.in.nalma
  }
}
alive.in.nalma.tot <- lapply(alive.in.nalma.tot, function (x) {names(x) <- NALMAs; x})

# make the matrix with the data of living species in each interval
mat.alive <- matrix(NA, ncol = length(NALMAs), nrow = length(spp))
rownames(mat.alive) <- spp
colnames(mat.alive) <- NALMAs
mat.alive.tot <- vector(mode = "list", length = 50) # here the 50 is the number of replicas you have
for (k in 1:50){ # here the 50 is the number of replicas you have
  mat.alive[,] <- NA
  for(i in 1:length(NALMAs)){# here the number 12 is the number of intervals you have
    mat.alive[which(rownames(mat.alive) %in% alive.in.nalma.tot[[k]][[i]]),i] <- 1
    mat.alive.tot[[k]] <- mat.alive
  }
}

# make a list such as the alive.in.nalma, now using only the species names which there is geographical occurrence - this one is fixed and does not change according to replica
nalma_pol_trues <- list(NA)
for (i in 1:length(NALMAs)){
  nalma_pol_trues[[i]] <- names(sp_split_nalmas)[grep(pattern = NALMAs[i], x = names(sp_split_nalmas))]
}
nalma_pol_trues_vector <- unlist(nalma_pol_trues)
names(nalma_pol_trues) <- NALMAs

# remove the names of the NALMAS from each species
nalma_pol_trues <- lapply(nalma_pol_trues, FUN = function (x) gsub(pattern = "\\..*", replacement = "", x = x))

# now, convert that list into matrix of geographic data, hence geo
mat.geo <- matrix(0, ncol = length(NALMAs), nrow = length(spp))
rownames(mat.geo) <- spp
colnames(mat.geo) <- NALMAs

for(i in 1:length(NALMAs)){
  mat.geo[which(rownames(mat.geo) %in% nalma_pol_trues[[i]]),i] <- 1
}

# combine the two matrices by matching
mat.geo.data.tot <- vector(mode = "list", length = 50) 
for(k in 1:50){ # here the 50 is the number of replicas you have
  mat.geo.data <- mat.geo == mat.alive.tot[[k]]
  mat.geo.data[which(mat.geo.data == T)] <- 1
  mat.geo.data.tot[[k]] <- mat.geo.data
}

# amount of missing data, how many gaps there are
sapply(mat.geo.data.tot, function (x) {length(which(x == 0))})

# how many SPECIES have missing data
sapply(mat.geo.data.tot,function(x) {table(rowSums(x == 0, na.rm = T) >= 1)})

# function to paste the row and column names when there is a given value
paste.row.col <- function(DF, value){
  ind <- which(DF == value, arr.ind=TRUE)
  paste(rownames(DF)[ind[,"row"]], colnames(DF)[ind[,"col"]], sep = '.')
}

# creates the final list of the points needed to add spatial data
missing.geo.nalma.tot <- vector(mode = "list", length = 50)# here the 50 is the number of replicas you have
for (k in 1:50){# here the 50 is the number of replicas you have
  missing.geo.nalma <- paste.row.col(DF = mat.geo.data.tot[[k]], value = 0)
  missing.geo.nalma.tot[[k]] <- missing.geo.nalma
}

# intermediate steps to determine which are the gaps
new.points.tot <- vector(mode = "list", length = 50)# here the 50 is the number of replicas you have
for (k in 1:50){# here the 50 is the number of replicas you have
  new.points.tot[[k]] <- rep(list(NA), length(missing.geo.nalma.tot[[k]]))
  names(new.points.tot[[k]]) <- missing.geo.nalma.tot[[k]]
}
sp_split_tot <- rep(list(sp_split_nalmas), times = 50)# here the 50 is the number of replicas you have
nalma_pol_trues_tot <- rep(list(nalma_pol_trues), times = 50)# here the 50 is the number of replicas you have
nalma_pol_trues_vector_tot <- rep(list(nalma_pol_trues_vector), times = 50)# here the 50 is the number of replicas you have
mat.geo.tot <- vector("list", 50)# here the 50 is the number of replicas you have
sp_split_tot_temp <- rep(list(vector("list", max(unlist(lapply(mat.geo.data.tot, function (x) max(rowSums(x == 0, na.rm = T))))))), times = 50)# here the 50 is the number of replicas you have
new.points.tot.temp <- rep(list(vector("list", max(unlist(lapply(mat.geo.data.tot, function (x) max(rowSums(x == 0, na.rm = T))))))), times = 50)# here the 50 is the number of replicas you have

# Now the code will determine what to do and when to do: if that point of the species distribution is a gap or if the interval does not have geographical data at all. There are some workarounds to errors if the species is the first in a list of empty names and such.

for (k in 1:50){# here the 50 is the number of replicas you have
  #print(k)
  for (j in 1:max(rowSums(mat.geo.data.tot[[k]] == 0, na.rm = T))){
    #print(j)  
    for (i in 1:nrow(mat.geo.data.tot[[k]])){
      #print(i)
      sp_split_tot_temp[[k]][[j]] <- c(sp_split_tot[[k]])
      seq1 <- which(NALMAs %in% names(mat.geo.data.tot[[k]][i,][complete.cases(mat.geo.data.tot[[k]][i,])])) # identify where the species is present by TSTE
      seq2 <- which(NALMAs %in% names(which(mat.geo.data.tot[[k]][i,][complete.cases(mat.geo.data.tot[[k]][i,])] == 0))) # identify where we don't have geographical data for that species
      
      if (length(seq2 == 0) == 0){ # no missing data, do nothing
        
      }  else if (seq2[1] < setdiff(seq1, seq2)[1]) { # logical test: is the position of the first missing nalma smaller than the first nalma with data? - just repeat the adjacent data
        new.points.tot[[k]][missing.geo.nalma.tot[[k]][grep(pattern = rownames(mat.geo.data.tot[[k]])[i], missing.geo.nalma.tot[[k]])]][which(NALMAs[match(gsub(pattern = ".*\\.", replacement = "", x = missing.geo.nalma.tot[[k]][grep(pattern = rownames(mat.geo.data.tot[[k]])[i], missing.geo.nalma.tot[[k]])]), NALMAs)] == NALMAs[setdiff(seq1,seq2)[1] - 1])] <- sp_split_tot[[k]][grep(pattern = rownames(mat.geo.data.tot[[k]])[i], nalma_pol_trues_vector_tot[[k]])[1]]
        
      } else if (seq2[1] > setdiff(seq1, seq2)[length(setdiff(seq1, seq2))]) { #logical test: is the position of the first missing nalma larger than the last nalma with data?
        new.points.tot[[k]][missing.geo.nalma.tot[[k]][grep(pattern = rownames(mat.geo.data.tot[[k]])[i], missing.geo.nalma.tot[[k]])]][which(NALMAs[match(gsub(pattern = ".*\\.", replacement = "", x = missing.geo.nalma.tot[[k]][grep(pattern = rownames(mat.geo.data.tot[[k]])[i], missing.geo.nalma.tot[[k]])]), NALMAs)] == NALMAs[setdiff(seq1,seq2)[length(setdiff(seq1,seq2))] + 1])] <- sp_split_tot[[k]][grep(pattern = rownames(mat.geo.data.tot[[k]])[i], nalma_pol_trues_vector_tot[[k]])[length(grep(pattern = rownames(mat.geo.data.tot[[k]])[i], nalma_pol_trues_vector_tot[[k]]))]]
        
      } else { # the last possible scenario is when there are gaps in the data, a species is missing from a nalma in the fossil record
        new.points.tot[[k]][missing.geo.nalma.tot[[k]][grep(pattern = rownames(mat.geo.data.tot[[k]])[i], missing.geo.nalma.tot[[k]])]][1] <- 
          st_centroid(
            st_convex_hull(
              st_union(st_union(
                sp_split_tot_temp[[k]][[1]][grep(pattern = rownames(mat.geo.data.tot[[k]])[i], names(sp_split_tot_temp[[k]][[1]]))][[grep(pattern = NALMAs[seq2[1] - 1], x = names(sp_split_tot_temp[[k]][[1]][grep(rownames(mat.geo.data.tot[[k]])[i], names(sp_split_tot_temp[[k]][[1]]))]))]]),
                st_union(sp_split_tot_temp[[k]][[1]][grep(pattern = rownames(mat.geo.data.tot[[k]])[i], names(sp_split_tot_temp[[k]][[1]]))][[grep(pattern = NALMAs[seq2[1] + 1], x = names(sp_split_tot_temp[[k]][[1]][grep(rownames(mat.geo.data.tot[[k]])[i], names(sp_split_tot_temp[[k]][[1]]))]))]]))
          )
        )
      }
    }
    
    # add the new.points
    new.points.tot.temp[[k]][[j]] <- Filter(Negate(anyNA), new.points.tot[[k]])
    sp_split_tot_temp[[k]][[j]] <- c(sp_split_tot[[k]], new.points.tot.temp[[k]][[j]])
    
    # now, when adding the nalma_pol_trues section, it will change between replicates
    
    for (l in 1:length(NALMAs)){
      nalma_pol_trues_tot[[k]][[l]] <- names(sp_split_tot_temp[[k]][[j]])[grep(pattern = NALMAs[l], x = names(sp_split_tot_temp[[k]][[j]]))]
    }
    
    nalma_pol_trues_tot[[k]] <- lapply(nalma_pol_trues_tot[[k]], FUN = function (x) gsub(pattern = "\\..*", replacement = "", x = x))
    names(nalma_pol_trues_tot[[k]]) <- NALMAs
    
    # now, convert that list into matrix of geographic data, hence geo
    mat.geo.tot[[k]] <- matrix(0, ncol = length(NALMAs), nrow = length(spp))
    rownames(mat.geo.tot[[k]]) <- spp
    colnames(mat.geo.tot[[k]]) <- NALMAs
    
    for(l in 1:length(NALMAs)){
      mat.geo.tot[[k]][which(rownames(mat.geo.tot[[k]]) %in% nalma_pol_trues_tot[[k]][[l]]),l] <- 1
    }
    
    # combine the two matrices by matching
    mat.geo.data.tot[[k]][,] <- NA
    mat.geo.data.tot[[k]] <- mat.geo.tot[[k]] == mat.alive.tot[[k]]
    mat.geo.data.tot[[k]][which(mat.geo.data.tot[[k]] == T)] <- 1
    
    missing.geo.nalma.tot[[k]] <- paste.row.col(DF = mat.geo.data.tot[[k]], value = 0)
  }
}

# vectors with the new points
new.sp.tot <- vector(mode = "list", 50)# here the 50 is the number of replicas you have
for (k in 1:50){# here the 50 is the number of replicas you have
  new.sp.tot[[k]] <- append(sp_split_nalmas, new.points.tot[[k]])
}

# detach the data points from a spdf object and reincorporate to the original dataset

new.cd.tot <- lapply(new.points.tot, function(x) lapply(x, function(y) {st_coordinates(y)}))

dt_points = vector(mode = "list", length = 50)
for (k in 1:length(new.cd.tot)){
  #print(k)
  for (i in 1:length(new.cd.tot[[k]])){
    #print(i)
    dt_points[[k]][[i]] <- data.frame(
      "lng" = new.cd.tot[[k]][[i]][,1], 
      "lat" = new.cd.tot[[k]][[i]][,2], 
      "species" = names(new.cd.tot[[k]][i]))
  }
}
    
new.cd.tot.2 <- lapply(dt_points, function (x) do.call(rbind, x))
for (k in 1:length(new.cd.tot.2)){
  rownames(new.cd.tot.2[[k]]) <- NULL
}
new.cd.tot.2 <- lapply(new.cd.tot.2, function (x) x[,c(3,1,2)])
new.cd.tot.2 <- lapply(new.cd.tot.2, function (x) separate(x, col = "species", c("species", "NALMA"), sep = "\\."))

df_can_tot <- vector("list", 50)# here the 50 is the number of replicas you have
for (i in 1:50){# here the 50 is the number of replicas you have
  df_can_tot[[i]] <- rbind.fill(df_can, new.cd.tot.2[[i]])
}

for (i in 1:50){# here the 50 is the number of replicas you have
  df_can_tot[[i]]$species <- factor(df_can_tot[[i]]$species, spp)
  df_can_tot[[i]]$NALMA <- factor(df_can_tot[[i]]$NALMA, NALMAs)
}

# Saving as object
dir.create("./Spatio_temporal_analyses/Spatial_coex/", recursive = T)
save(df_can_tot, file = "./Spatio_temporal_analyses/Spatial_coex/df_can_tot.Rdata", compress = "xz")
