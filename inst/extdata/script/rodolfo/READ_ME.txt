# This README provides details about the data and scripts used in the study titled "Ecological and spatial overlap indicate interspecific competition during North American Canid radiation" 
by Rodolfo P. Graciotti, Tiago B. Quental, Lucas M. V. Porto and Salatiel Gonçalves-Neto.

## Please refer to this guide for the directory structure created while running R scripts:
# MAIN: main directory for all analyses
# PBDB: occurrence data
# PyRate: your PyRate version directory from https://github.com/dsilvestro/PyRate/blob/master/ with pyratelibs and utilities.
# PyRate_analyses: inputs and results from PyRate baseline analyses
# Spatio_temporal_analyses: inputs and results from coexistence metrics
    Diversity_curves: matrix data of standing diversity    
    Time_coex: matrix data of time coexistence
    Spatial_coex: matrix data of spatial coexistence at different scales
    Time_space: time series of coexistence through time and space at different scales
# Ecomorphological data
    morphospace: body mass, diet and inputs and resulst for trait and LDA analyses
    competition: time series of competition in different scales
# Continuous: results from PyRateContinuous analyses. Suggested structure:
    time_series: paste all time series to this directory
    log_files: paste all .log files to this directory
    
# Data:
## Occurrence Data:
pbdb_data_canidae_NA_full.csv: Full Canidae occurrence dataset downloaded from Paleobiology Database in August 2024;

## Ecomorphological Data:
Canidae_trees.tree: 1000 phylogenetic trees for Canidae, pruned from Carnivora trees of Faurby et al. (2024), used for data imputation on species with missing craniodental variables;

eco_morph.csv: Includes original ecomorphological data for 6 craniodental variables from Slater (2015), Faurby et al. (2021) and Juhn et al. (2024), featuring missing data. Used for data imputation;

canids_body_mass.csv: Log transformed body mass for all 133 species from Faurby et al. (2021);

Training_set.csv: Contains ecomorphological data for 6 craniodental variables used in Linear Discriminant Analysis. Comprises 25 extant canid species, nine procyonids, and Ailurus fulgens. Used to categorize canid species into different species categories. Original data taken from Slater (2015);

nalmas.txt: text file with the intervals defined for PyRate qShift notation;

longevities.RData: R objects compiling the true times of speciation and extinction for each replica.

# R and PyRate scripts:
01_PBDB_data_manipulation.R: Curatorial work for occurrence data;

02_Generate_PyRate_input.R: Creates input for PyRate from cleaned occurrence data;

03_Dataset_spatial_prep.R: Assign occurrence data to NALMA;

04_TSTE_tables.R: Compile true times of speciation and extinction;
## We also made available the .RData file resulting from this script.

05_Gap_filling.R: Estimate spatial distribution for time intervals which species were not sampled;

06_Temporal_coexistence.R: Create the temporal coexistence matrices and time series; 

07_Spatial_coexistence.R: Create the spatial coexistence matrices and time series;

08_Time_space.R: Combine temporal and spatial coexistence matrices and creates time series;

09_Trait_imputation.R: Input missing craniodental variables;

10_Linear_Discriminant_Analysis.R: Perform LDA for craniodental variables;

11_Temperature.R: Extract temperature time series;

12_Competition_metrics.R: Create competition metrics time series;

13_Time_series_plot.R: Compile and plot time series data;

14_PyRateContinuous_results.R: Compile results from PyRateContinuous;
## This script is customizable depending on the time window used on the PyRateContinuous analyses.
## We discussed the epoch time scale in the main text and "cenozoic" scale on ESM 

15_PyRateContinuous_plots.R: Create figures of PyRateContinuous results;

16_PyRate_Commands.txt: Contains commands used in PyRate and PyRateContinuous analyses. Detailed instructions for PyRate setup and analysis can be found in the PyRate GitHub tutorials (link: https://github.com/dsilvestro/PyRate/blob/master/tutorials/README.md).
