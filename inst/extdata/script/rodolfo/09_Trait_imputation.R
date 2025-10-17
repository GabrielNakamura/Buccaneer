# Data imputation procedure ####

rm(list = ls());gc()

# Load necessary packages ####
require(phytools)
require(geiger)
require(here)
require(Rphylopars)
require(purrr)
require(dplyr)

# Canidae tree and morphological data ####
Canidae_trees <- read.tree(file = "./Canidae_trees.tree")

# Read the raw data from CSV
df <- read.csv(file = "./eco_morph.csv", header = TRUE, sep = ";")
df <- df %>%
  mutate(species = case_when(
    species == "Aenocyon_dirus" ~ "Canis_dirus",
    species == "Eucyon_ferox" ~ "Canis_ferox",
    species== "Ferrucyon_avius" ~ "Cerdocyon_avius",
    TRUE ~ species  
  ))
data <- df[, !names(df) %in% c("status")] 

# Create a list to store ordered data
list_data <- list()

# Loop to order species data to match the tree tips
for (i in 1:length(Canidae_trees)) {
  data_order <- data[match(Canidae_trees[[i]]$tip.label, data$species), ]
  list_data[[i]] <- data_order
}

# Create lists to store imputation data and lambda reconstructions
imputation_data <- list()
p_lambda_ancrecon <- list()

# Loop through the trees and perform imputation
# Warning messages from model testing while loop runs, it's expected
# Warning - it can take a long time to run
for (i in 1:length(Canidae_trees)) {
  p_lambda <- phylopars(
    trait_data = list_data[[i]],
    tree = Canidae_trees[[i]],
    phylo_correlated = TRUE,
    pheno_correlated = FALSE,
    model = "lambda"
  )
  p_lambda_ancrecon[[i]] <- p_lambda$anc_recon[1:133, ][order(rownames(p_lambda$anc_recon[1:133, ])), ]
  imputation_data[[i]] <- p_lambda
}

# Create a list to store lambda reconstructions as data frames
p_lambda_dataframe <- list()

# Convert lambda reconstructions to data frames
for (i in 1:length(Canidae_trees)) {
  p_lambda_dataframe[[i]] <- as.data.frame(p_lambda_ancrecon[[i]])
}

# Calculate row means for specified columns and create a data frame
columns_to_mean <- c("RLGA", "RBL", "M2S", "M1BS", "p4S")
species_list_mean <- p_lambda_dataframe %>%
  reduce(`+`) %>%
  mutate(species = rownames(p_lambda_dataframe[[1]])) %>%
  select(species, everything()) %>%
  mutate(across(all_of(columns_to_mean), ~ . / length(p_lambda_dataframe)))
rownames(species_list_mean) <- NULL
species_list_mean <- species_list_mean[, -2]

df_list <- as.data.frame(do.call(cbind, species_list_mean))
colnames(df_list) <- columns_to_mean
row.names(df_list) <- row.names(p_lambda_dataframe[[1]])

# Add the MAT column from the original data frame
df_list <- cbind(species_list_mean)

# Write the final data frame to a CSV file for all species ####
# Creates directory for output data
dir.create("./Ecomorphological_data/morphospace/", recursive = T)
write.table(df_list, "./Ecomorphological_data/morphospace/input_values.csv", row.names = FALSE, sep = ",")

# Preparing the fossil data for the LDA analysis ####

ed <- read.csv("./Training_set.csv", header = TRUE, sep = ";")

fossil_species <- df_list %>% 
  filter(!species %in% ed$species)

write.table(fossil_species, "./Ecomorphological_data/morphospace/fossil_data.csv", row.names = FALSE, sep = ",")
saveRDS(p_lambda_ancrecon, "./Ecomorphological_data/morphospace/output_imputation.rds")

