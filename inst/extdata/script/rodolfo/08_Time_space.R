# Multiplying spatial and time coexistence matrices ####
# loading prep data ####

rm(list=ls());gc()

library(FSA)
library(tidyverse)
library(data.table)

load("./Spatio_temporal_analyses/Time_coex/mat_time_replicas.Rdata")
load(here::here("inst", "extdata", "script", "rodolfo", "mat_time_replicas.Rdata"))
load("./Spatio_temporal_analyses/Diversity_curves/div_curves.Rdata")
load(here::here("inst", "extdata", "script", "rodolfo", "div_curves.Rdata"))
load("./Spatio_temporal_analyses/Spatial_coex/df_can_tot.Rdata")
load(here::here("inst", "extdata", "script", "rodolfo", "df_can_tot.Rdata"))
load("./Spatio_temporal_analyses/Time_coex/longevities.RData")
load(here::here("inst", "extdata", "script", "rodolfo", "longevities.RData"))
dir.create("./Spatio_temporal_analyses/Time_space/reach/", recursive = T)
dir.create("./Spatio_temporal_analyses/Time_space/site/", recursive = T)

div.curves <- lapply(div.curves, unname)

# determine the duration of the bins and the scale
# in order to avoid overestimating species true geographical data, used an approximation of root age for each replica. the formula sets the root age as the next interval in 0.1 myr scale from the speciation time of the first species
root_ages <- sapply(longs, function (x) floor(max(x$TS) * 10) /10) #

for (k in 1:length(longs)){
  longs[[k]]$TS[which.max(longs[[k]]$TS)] <- root_ages[[k]]
}

# apply diagonal to 0
for (k in 1:50){ # here the "50" is the number of replicas you have
  #print(k)
  for (i in 1:length(mat.time.replicas[[k]])){
    diag(mat.time.replicas[[k]][[i]]) <- 0
  }
}

# defining the function that extracts the summary statistics ####
smry.coex <- function(x){
  res <- matrix(NA, ncol(x), 6)
  for(i in 1:ncol(x)){
    res[i,] <- c(summary(x[x[,i]>0, i]))
    res <- data.frame(res)
    colnames(res) <- names(summary(rnorm(50,mean = 1)))# here the "50" is the number of replicas you have
  }
  return(res)
}

# defining the function that summarizes the tables ####
smry.table <- function(x){
  for (i in 1:length(x)){
    x[[i]] <- x[[i]][x[[i]]>0]
  }
  y <- lapply(x, Summarize)
  z <- data.frame(t(do.call(cbind.data.frame, y)))
  setDT(z, keep.rownames = "id")[]
  z$id <- as.numeric(z$id)
  return(z)
}

#NALMAS
NALMA <- c("Duchesnean", "Chadronian", "Orellan", "Whitneyan", "Arikareean", "Hemingfordian", "Barstovian", "Clarendonian", "Hemphillian", "Blancan", "Irvingtonian", "Rancholabrean_Present") # Barnsosky 2014
NALMA_age <- c(39.7, 37, 33.9, 31.8, 29.5, 18.5, 16.3, 12.5, 9.4, 4.7, 1.4, 0.21,0) # Barnosky 2014

# datasets when using nalmas
points.in.time <- round(NALMA_age,1)

# subsetting the temporal matrices per nalma
duc.list <- lapply(mat.time.replicas, function (x) x[as.character(seq(from = points.in.time[1], to = points.in.time[2]+0.1, by = -0.1))])
chad.list <- lapply(mat.time.replicas, function (x) x[as.character(seq(from = points.in.time[2], to = points.in.time[3]+0.1, by = -0.1))])
ore.list <- lapply(mat.time.replicas, function (x) x[as.character(seq(from = points.in.time[3], to = points.in.time[4]+0.1, by = -0.1))])
whi.list <- lapply(mat.time.replicas, function (x) x[as.character(seq(from = points.in.time[4], to = points.in.time[5]+0.1, by = -0.1))])
ari.list <- lapply(mat.time.replicas, function (x) x[as.character(seq(from = points.in.time[5], to = points.in.time[6]+0.1, by = -0.1))])
hem.list <- lapply(mat.time.replicas, function (x) x[as.character(seq(from = points.in.time[6], to = points.in.time[7]+0.1, by = -0.1))])
bar.list <- lapply(mat.time.replicas, function (x) x[as.character(seq(from = points.in.time[7], to = points.in.time[8]+0.1, by = -0.1))])
cla.list <- lapply(mat.time.replicas, function (x) x[as.character(seq(from = points.in.time[8], to = points.in.time[9]+0.1, by = -0.1))])
hemp.list <- lapply(mat.time.replicas, function (x) x[as.character(seq(from = points.in.time[9], to = points.in.time[10]+0.1, by = -0.1))])
bla.list <- lapply(mat.time.replicas, function (x) x[as.character(seq(from = points.in.time[10], to = points.in.time[11]+0.1, by = -0.1))])
irv.list <- lapply(mat.time.replicas, function (x) x[as.character(seq(from = points.in.time[11], to = points.in.time[12]+0.1, by = -0.1))])
ran.list <- lapply(mat.time.replicas, function (x) x[as.character(seq(from = points.in.time[12], to = points.in.time[13], by = -0.1))])


# Time_space multiplications ####


# site ####

load("./Spatio_temporal_analyses/Spatial_coex/site/mat_site_tot.Rdata")
load(here::here("inst", "extdata", "script", "rodolfo", "mat_site_tot.Rdata"))

NALMAs.matrices <- mat.site.tot
NALMAs.matrices <- lapply(NALMAs.matrices, as.matrix)

# multiplying
# the multiplication is easy, please note it is NOT a PRODUCT of matrices, an algebraic multiplication;
# it is a vectorial multiplication, by each corresponding cell.

mat.site.duc <- lapply(duc.list, function (x) lapply(x, function (x) x * NALMAs.matrices[[1]]))
mat.site.chad <- lapply(chad.list, function (x) lapply(x, function (x) x * NALMAs.matrices[[2]]))
mat.site.ore <- lapply(ore.list, function (x) lapply(x, function (x) x * NALMAs.matrices[[3]]))
mat.site.whi <- lapply(whi.list, function (x) lapply(x, function (x) x * NALMAs.matrices[[4]]))
mat.site.ari <- lapply(ari.list, function (x) lapply(x, function (x) x * NALMAs.matrices[[5]]))
mat.site.hem <- lapply(hem.list, function (x) lapply(x, function (x) x * NALMAs.matrices[[6]]))
mat.site.bar <- lapply(bar.list, function (x) lapply(x, function (x) x * NALMAs.matrices[[7]]))
mat.site.cla <- lapply(cla.list, function (x) lapply(x, function (x) x * NALMAs.matrices[[8]]))
mat.site.hemp <- lapply(hemp.list, function (x) lapply(x, function (x) x * NALMAs.matrices[[9]]))
mat.site.bla <- lapply(bla.list, function (x) lapply(x, function (x) x * NALMAs.matrices[[10]]))
mat.site.irv <- lapply(irv.list, function (x) lapply(x, function (x) x * NALMAs.matrices[[11]]))
mat.site.ran <- lapply(ran.list, function (x) lapply(x, function (x) x * NALMAs.matrices[[12]]))


mat.sites.nalma <- vector("list", 50) # here the "50" is the number of replicas you have
# recombine them
for (k in 1:50){ # here the "50" is the number of replicas you have
  mat.sites.nalma[[k]] <- list(mat.site.duc[[k]],mat.site.chad[[k]], mat.site.ore[[k]], mat.site.whi[[k]], mat.site.ari[[k]], mat.site.hem[[k]], mat.site.bar[[k]], mat.site.cla[[k]], mat.site.hemp[[k]], mat.site.bla[[k]], mat.site.irv[[k]], mat.site.ran[[k]])
}
mat.sites.nalma <- lapply(mat.sites.nalma, function (x) do.call(c, x))

# after the multiplication and recombining ensures all have the same length, the next step determines the number of species coexisting in the method. It is therefore needed to remove empty matrices produced due to not all replicas including the first NALMA.

for (k in 1:50){
    mat.sites.nalma[[k]] <- mat.sites.nalma[[k]][lapply(mat.sites.nalma[[k]], length) > 0]
}

site.coex.nalma <- vector("list", 50)# here the "50" is the number of replicas you have
for (k in 1:50){# here the "50" is the number of replicas you have
  site.coex.nalma[[k]] <- rep(list(rep(NA, 133)),length(mat.sites.nalma[[k]]))
  names(site.coex.nalma[[k]]) <- names(mat.sites.nalma[[k]])
}


# number of coexisting species - reading each line (species) of each time slice of each replica and sum the number of coexistences
for (k in 1:50){# here the "50" is the number of replicas you have
  #print(k)
  for(i in 1:length(mat.sites.nalma[[k]])){
    for(j in 1:length(spp)){
      site.coex.nalma[[k]][[i]][[j]] <- sum(mat.sites.nalma[[k]][[i]][j,])
    }
  }
}

site.coex.nalma.mat <- lapply(site.coex.nalma, function (x) do.call(cbind.data.frame, x))
for (k in 1:50){# here the "50" is the number of replicas you have
  rownames(site.coex.nalma.mat[[k]]) <- spp
}

# summarizes the coexistence values
smry.site.nalma <- lapply(site.coex.nalma.mat, smry.coex)
smry.site.nalma.t <- lapply(site.coex.nalma, smry.table)

# this one we used to analyze the proportion of coexisting species to the total, hence the "prop".
# you can ignore it, but some of its by products will be used in the morphospace multiplications

site.coex.nalma.prop <- vector("list", 50)# here the "50" is the number of replicas you have
for (k in 1:50){# here the "50" is the number of replicas you have
  for (i in 1:length(site.coex.nalma[[k]])){
    site.coex.nalma.prop[[k]][[i]] <- site.coex.nalma[[k]][[i]]/div.curves[[k]][[i]]
  }
}

for (k in 1:50){# here the "50" is the number of replicas you have
  for (i in 1:length(site.coex.nalma.prop[[k]])){
    site.coex.nalma.prop[[k]][[i]][is.nan(site.coex.nalma.prop[[k]][[i]])] <- 0
  }
}

smry.site.nalma.t.prop <- lapply(site.coex.nalma.prop, smry.table)
# the warnings produced by the line above is not a concern. It happens because at some points in time there is just one species.

for (k in 1:50){# here the "50" is the number of replicas you have
  smry.site.nalma.t.prop[[k]]$id <- smry.site.nalma.t[[k]]$id
}

save(mat.sites.nalma, smry.site.nalma.t, smry.site.nalma.t.prop, file = "./Spatio_temporal_analyses/Time_space/site/site_final.Rdata", compress = "xz")
save(mat.sites.nalma, smry.site.nalma.t, smry.site.nalma.t.prop, file = here::here("inst", "extdata", "script", "rodolfo", "site_final.Rdata"), compress = "xz")

# here you get the diversity of the "site" scenario

#for (k in 1:50){ # here the "50" is the number of replicas you have
#  aux_1 <- smry.site.nalma.t[[k]][,c(1,3)]
#  names(aux_1) <- c("time", paste0("site_diversity"))
#  aux_1 <- aux_1[order(aux_1$time, decreasing = F),]
#  aux_1[,2][which(is.na(aux_1[,2]))] <- 0
#  write.table(x = aux_1, file = paste0("./Spatio_temporal_analyses/Time_space/site/site_diversity", k, ".txt"), row.names = F, quote = F)
#}

for (k in 1:50){ # here the "50" is the number of replicas you have
  aux_1 <- smry.site.nalma.t[[k]][,c(1,3)]
  names(aux_1) <- c("time", paste0("site_diversity"))
  aux_1 <- aux_1[order(aux_1$time, decreasing = F),]
  aux_1[,2][which(is.na(aux_1[,2]))] <- 0
  write.table(x = aux_1,
              file = here::here("inst", "extdata", "script", "rodolfo", "site_res", paste0("site_diversity", k, ".txt")),
              row.names = F, quote = F)
}



# reach ####

# ATTENTION
# the exact same rationale of the "site" scenario will be applied here to multiplication
# using "reach" matrices. Thus, the code that comes now is not commended in as much detail.

load("./Spatio_temporal_analyses/Spatial_coex/reach/mat_reach_tot.Rdata")

NALMAs.matrices <- mat.times.ab.tot

# multiplying

mat.reach.duc <- vector("list", 50)
mat.reach.chad <- vector("list", 50) # here the "50" is the number of replicas you have
mat.reach.ore <- vector("list", 50) # here the "50" is the number of replicas you have
mat.reach.whi <- vector("list", 50) # here the "50" is the number of replicas you have
mat.reach.ari <- vector("list", 50) # here the "50" is the number of replicas you have
mat.reach.hem <- vector("list", 50) # here the "50" is the number of replicas you have
mat.reach.bar <- vector("list", 50) # here the "50" is the number of replicas you have
mat.reach.cla <- vector("list", 50) # here the "50" is the number of replicas you have
mat.reach.hemp <- vector("list", 50) # here the "50" is the number of replicas you have
mat.reach.bla <- vector("list", 50) # here the "50" is the number of replicas you have
mat.reach.irv <- vector("list", 50) # here the "50" is the number of replicas you have
mat.reach.ran <- vector("list", 50) # here the "50" is the number of replicas you have

for (k in 1:50){# here the "50" is the number of replicas you have
  #print(k)
mat.reach.duc[[k]] <- lapply(duc.list[[k]], function (x) x * NALMAs.matrices[[k]][[1]])
mat.reach.chad[[k]] <- lapply(chad.list[[k]], function (x) x * NALMAs.matrices[[k]][[2]])
mat.reach.ore[[k]] <- lapply(ore.list[[k]], function (x) x * NALMAs.matrices[[k]][[3]])
mat.reach.whi[[k]] <- lapply(whi.list[[k]], function (x) x * NALMAs.matrices[[k]][[4]])
mat.reach.ari[[k]] <- lapply(ari.list[[k]], function (x) x * NALMAs.matrices[[k]][[5]])
mat.reach.hem[[k]] <- lapply(hem.list[[k]], function (x) x * NALMAs.matrices[[k]][[6]])
mat.reach.bar[[k]] <- lapply(bar.list[[k]], function (x) x * NALMAs.matrices[[k]][[7]])
mat.reach.cla[[k]] <- lapply(cla.list[[k]], function (x) x * NALMAs.matrices[[k]][[8]])
mat.reach.hemp[[k]] <- lapply(hemp.list[[k]], function (x) x * NALMAs.matrices[[k]][[9]])
mat.reach.bla[[k]] <- lapply(bla.list[[k]], function (x) x * NALMAs.matrices[[k]][[10]])
mat.reach.irv[[k]] <- lapply(irv.list[[k]], function (x) x * NALMAs.matrices[[k]][[11]])
mat.reach.ran[[k]] <- lapply(ran.list[[k]], function (x) x * NALMAs.matrices[[k]][[12]])
}

mat.reach.nalma <- vector("list", 50)# here the "50" is the number of replicas you have
for (k in 1:50){# here the "50" is the number of replicas you have
  mat.reach.nalma[[k]] <- list(mat.reach.duc[[k]], mat.reach.chad[[k]], mat.reach.ore[[k]], mat.reach.whi[[k]], mat.reach.ari[[k]], mat.reach.hem[[k]], mat.reach.bar[[k]], mat.reach.cla[[k]], mat.reach.hemp[[k]], mat.reach.bla[[k]], mat.reach.irv[[k]], mat.reach.ran[[k]])
}
mat.reach.nalma <- lapply(mat.reach.nalma, function (x) do.call(c, x))

# after the multiplication and recombining ensures all have the same length, the next step determines the number of species coexisting in the method. It is therefore needed to remove empty matrices produced due to not all replicas including the first NALMA.

for (k in 1:50){
  mat.reach.nalma[[k]] <- mat.reach.nalma[[k]][lapply(mat.reach.nalma[[k]], length) > 0]
}


reach.coex.nalma <- vector("list", 50)# here the "50" is the number of replicas you have
for (k in 1:50){# here the "50" is the number of replicas you have
  reach.coex.nalma[[k]] <- rep(list(rep(NA, 133)),length(mat.reach.nalma[[k]]))
  names(reach.coex.nalma[[k]]) <- names(mat.reach.nalma[[k]])
}

for (k in 1:50){# here the "50" is the number of replicas you have
  #print(k)
  for(i in 1:length(mat.reach.nalma[[k]])){
    for(j in 1:length(spp)){
      reach.coex.nalma[[k]][[i]][[j]] <- sum(mat.reach.nalma[[k]][[i]][j,])
    }
  }
}

reach.coex.nalma.mat <- lapply(reach.coex.nalma, function (x) do.call(cbind.data.frame, x))
for (k in 1:50){# here the "50" is the number of replicas you have
  rownames(reach.coex.nalma.mat[[k]]) <- spp
}

smry.reach.nalma <- lapply(reach.coex.nalma.mat, smry.coex)
smry.reach.nalma.t <- lapply(reach.coex.nalma, smry.table)

reach.coex.nalma.prop <- vector("list", 50)# here the "50" is the number of replicas you have
for (k in 1:50){# here the "50" is the number of replicas you have
  for (i in 1:length(reach.coex.nalma[[k]])){
    reach.coex.nalma.prop[[k]][[i]] <- reach.coex.nalma[[k]][[i]]/div.curves[[k]][[i]]
  }
}

for (k in 1:50){# here the "50" is the number of replicas you have
  for (i in 1:length(reach.coex.nalma.prop[[k]])){
    reach.coex.nalma.prop[[k]][[i]][is.nan(reach.coex.nalma.prop[[k]][[i]])] <- 0
  }
}

smry.reach.nalma.t.prop <- lapply(reach.coex.nalma.prop, smry.table)
# the warnings produced by the line above is not a concern. It happens because at some points in time there is just one species.

for (k in 1:50){# here the "50" is the number of replicas you have
  smry.reach.nalma.t.prop[[k]]$id <- smry.reach.nalma.t[[k]]$id
}

save(mat.reach.nalma, smry.reach.nalma.t, smry.reach.nalma.t.prop, file = "./Spatio_temporal_analyses/Time_space/reach/reach_final.Rdata", compress = "xz")

# here you get the diversity of the "reach" scenario

for (k in 1:50){ # here the "50" is the number of replicas you have
  aux_1 <- smry.reach.nalma.t[[k]][,c(1,3)]
  names(aux_1) <- c("time", paste0("reach_diversity"))
  aux_1 <- aux_1[order(aux_1$time, decreasing = F),]
  aux_1[,2][which(is.na(aux_1[,2]))] <- 0
  write.table(x = aux_1, file = paste0("./Spatio_temporal_analyses/Time_space/reach/reach_diversity", k, ".txt"), row.names = F, quote = F)
}

