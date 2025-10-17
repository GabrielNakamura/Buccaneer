# Preparing for PyRate analyses using pyrate_utilities.r ####

rm(list=ls()); gc()

# Refer to https://github.com/dsilvestro/PyRate/blob/master/tutorials/pyrate_tutorial_1.md on how to download and install PyRate

source("./PyRate/pyrate_utilities.r") # load from your pyrate directory
library(dplyr)

# PBDB_data ####

canidae <- read.table(file = "./PBDB/PBDB_NEW_Canidae_NA_simpler_PyRate_high_resol.csv", header = T, sep=",")

spp <- sort(unique(canidae$Species))

# Accounting for age dependence in multiple occurrences from the same site: what happens is that occurrences from the same assemblage/site are co-occurrences and their ages are probably correlated. To assign different ages per occurrence would not be realistic and induce pseudoreplicates, as they are not informative and/or redundant. 

# Site information is defined by the collection_no.

# factoring Site number to simplify
canidae$Site <- as.numeric(factor(canidae$Site))

# Extant and extinct species ####
# some species might be extant, but have no fossil data close to the present, it is best to write all the living spp manually.

extant_canidae <- c("Canis_latrans", "Canis_lupus", "Canis_rufus", "Urocyon_cinereoargenteus", "Urocyon_littoralis", "Vulpes_lagopus", "Vulpes_macrotis", "Vulpes_velox", "Vulpes_vulpes")

# now, the %in% gives the correct test and ifelse gives the correct option of identifying two possible outcomes
canidae$Status <- 
  ifelse(test = canidae$Species %in% extant_canidae, yes = "extant", no = "extinct")

# Final tables ####

canidae_input <- canidae[,c("Species", "Status", "MinT", "MaxT", "Site")]
# remove the same site occurrences
canidae_input <- unique(canidae_input)

# writing tables
# create directory for PyRate_analyses and inputs
dir.create("./PyRate_analyses/input_PyRate/", recursive = T)
filename <- paste("./PyRate_analyses/input_PyRate/Canidae_N_america.txt", sep = "")
write.table(canidae_input, filename, sep = "\t", row.names = F)

# function from silvestro that extracts the final files
extract.ages(file = filename, replicates = 50)

###### Data is ready to run PyRate analyses, refer to "16_PyRate_Commands.txt"