# Linear discriminant analysis (LDA)
# The script here is strongly based on the script used by Slater (2015)

rm(list = ls());gc()
require(MASS)

# Loading data ####
# Load training data from a CSV file
# ed <- read.csv("./Training_set.csv", header = TRUE, sep = ";", row.names = "species")
ed <- read.csv(here::here("inst",
                          "extdata",
                          "script",
                          "rodolfo",
                          "Training_set.csv"),
               header = TRUE, sep = ";", row.names = "species")

# Fossil data
# fd <- read.csv("./fossil_data.csv", header = TRUE, sep = ",",row.names = "species")
fd <- read.csv(here::here("inst",
                          "extdata",
                          "script",
                          "rodolfo",
                          "Ecomorphological_data",
                          "morphospace",
                          "fossil_data.csv"),
               header = TRUE, sep = ",", row.names = "species")
na.foo <- function(x) sum(!is.na(x))

# Perform the LDA in the training set ####

# Training set LDA
disc <- lda(diet ~ RLGA + RBL + M1BS + M2S + p4S, data = as.data.frame(ed))
dis.pred <- predict(disc)
disc

# Visualization and storage of LDA results
resultados_lda <- cbind(rownames(ed), ed[, "diet"], dis.pred$class)
plot(dis.pred$x, pch = 21, bg = dis.pred$class)
text(dis.pred$x, label = rownames(dis.pred$x), cex = 0.4)
coef(disc)

# Selection of extant species of interest we have fossil occurrences
extant <- cbind(dis.pred$class, dis.pred$posterior, dis.pred$x)
canids_extant <- extant[c("Canis_latrans", "Canis_lupus","Cuon_alpinus", "Urocyon_cinereoargenteus","Urocyon_littoralis", "Vulpes_lagopus", "Vulpes_velox", "Vulpes_vulpes"), c(1, 5, 6)]

# Classification of fossil data
fd.pred <- predict(disc, as.data.frame(fd))
fossil.classif <- cbind(fd.pred$class, fd.pred$posterior, fd.pred$x)
fossil.classif1 <- fossil.classif[-which(is.na(fossil.classif[, 1])), ]
unclass<- which(is.na(fossil.classif[, 1]))
apply(fd[unclass, ], 2, na.foo)

# Final results of discriminant classification for the fossil species
discr_result <- rbind(fossil.classif[, c(1, 5, 6)])

# Obtaining diet information
diet <- discr_result[match(rownames(fd), rownames(discr_result)), 1]

# Create the complete dataset (133 extinct and extant species) with diet information
ft <- as.data.frame(rbind(fossil.classif[, c(1, 5, 6)], canids_extant))
ft <- ft[order(row.names(ft)), ]
ft$species <- rownames(ft)
ft <- ft[, c("species", setdiff(names(ft), "species"))]
colnames(ft) <- c("species", "diet", "LD1", "LD2")

# write.table(ft, "./Ecomorphological_data/morphospace/LDA_results.csv", row.names = F, sep = ",")
write.table(ft, here::here("inst",
                           "extdata",
                           "script",
                           "rodolfo",
                           "Ecomorphological_data",
                           "morphospace",
                           "LDA_results.csv"),
            row.names = F, sep = ",")

# Body mass dataset
# dat_mass <- read.csv("./canids_body_mass.csv", header = TRUE, sep = ",")
dat_mass <- read.csv(here::here("inst",
                                "extdata",
                                "script",
                                "rodolfo",
                                "canids_body_mass.csv"),
                     header = TRUE, sep = ",")

combined_df <- as.data.frame(cbind(
  species = dat_mass$species,
  diet = ft$diet,
  LD1 = ft$LD1,
  log_mass = dat_mass$Size_Est
))

# write.table(combined_df, "./Ecomorphological_data/morphospace/final_data.csv", row.names = F, sep = ",")
write.table(combined_df, here::here("inst",
                                    "extdata",
                                    "script",
                                    "rodolfo",
                                    "Ecomorphological_data",
                                    "morphospace",
                                    "final_data.csv"),
            row.names = F, sep = ",")
