# Reading true times of speciation and extinction from the previous PyRate analyses ####

#### WARNING  We are making available the .RData resulting from this script (longevities.RData), so you can skip this script entirely. ####

rm(list = ls());gc()

library(gtools)

# After step 03 described in the 16_PyRate_Commands.txt file, we compile the true times of speciation and extinction for each PyRate replica to set the basis for all analyses of coexistence in time in a .RData file.

# For the sake of completeness, if you wish to compile the results, it is advised to maintain the se_est.txt files resulting from the PyRate -ginput command in a separate directory.

path <- ("./PyRate_analyses/logs/TSTE") # path to your _se_ext.txt files

# Loading the individual TSTE tables ####
list.filenames <- mixedsort(list.files(path, full.names = T))
longevities <- list()
for (i in 1:50){
  longevities[[i]] <- read.table(list.filenames[i], header = T)
}

spp <- read.table("./PyRate_analyses/input_PyRate/Canidae_N_america_TaxonList.txt", header = T)
spp <- spp[,1]

# unifying and formatting longs as the functions use it

longs <- lapply(X = longevities, FUN = function(x) data.frame(TS = x$ts, TE = x$te, row.names = spp))
longs.order <- lapply(X = longs, function(x) x[order(x$TS, x$TE, decreasing = T),])

rm(longevities)

dir.create("~/MAIN/Spatio_temporal_analyses/Time_coex", recursive = T)
save(list = c("longs", "longs.order", "spp"), file = "~/MAIN/Spatio_temporal_analyses/Time_coex/longevities.RData", compress = "xz")
