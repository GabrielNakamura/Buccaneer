# Preparing dataset to spatial analyses ####
# Original dataset compiled from PDBD of all occurrences
# After steps 01 to 03 from "16_PyRate_Commands.txt, it is advised to set different directories for combined and individual se_est files"

rm(list = ls()); gc()

library(tidyverse)

# loading data ####
# dataset table

df_can <- read.table(file = "./PBDB/PBDB_NEW_Canidae_NA_simpler_PyRate_high_resol.csv", header = T, sep = ",")

# rename and reorder to simplify

colnames(df_can) <- c("species", "max", "min", "lng", "lat", "site")
df_can$site <- as.numeric(factor(df_can$site))

head(df_can)

# Defining intervals for each occurrence based on their range
# NALMAs vector (or the specific time intervals that you are interested)

#NALMAS
NALMA <- c("Duchesnean","Chadronian","Orellan","Whitneyan","Arikareean","Hemingfordian","Barstovian","Clarendonian","Hemphillian","Blancan","Irvingtonian","Rancholabrean_Present") # Barnsosky 2014
NALMA_age <- c(39.7, 37, 33.9, 31.8, 29.5, 18.5, 16.3, 12.5, 9.4, 4.7, 1.4, 0.21,0) # Barnosky 2014

# operations to define NALMA
midpoints <- NALMA_age[-length(NALMA_age)] + diff(NALMA_age)/2

# defining NALMA based on range of each occurrences
for(k in 1:length(NALMA)){
 df_can$NALMA[df_can$min >= NALMA_age[k+1] & df_can$max <= NALMA_age[k]] <- NALMA[k]
}

# occurrences which do not fit NALMA
length(which(is.na(df_can$NALMA)))

# use midpoint of range to define NALMA for these occurrences
df_can$midpoint <- rowMeans(df_can[,c('max', 'min')])

missing <- df_can[which(is.na(df_can$NALMA)),]

for(k in 1:length(NALMA)){
  missing$NALMA[missing$midpoint >= NALMA_age[k+1] & missing$midpoint <= NALMA_age[k]] <- NALMA[k]
}

missing[,c(1,2,3,7,8)]

# replace
df_can[which(is.na(df_can$NALMA)),'NALMA'] <- missing$NALMA
length(which(is.na(df_can$NALMA)))

# saving ####

str(df_can)

save(df_can, file = "./PBDB/df_can.Rdata", compress = "xz")
