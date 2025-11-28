# pre-processing-data

``` r
library(buccaneer)
```

## Background

In this tutorial I will present how data on longevity and fossil
occurrence can be used to investigate the imprints of species
competition in different levels. This will be done with the help of the
new R package `Buccaneer`

## Data

We will use the data from [Faurby et al
(2024)](https://royalsocietypublishing.org/doi/10.1098/rspb.2024.0473)
to demonstrate how the new Buccaneer R package can be used in processing
and data analysis.

This data contains information on fossil data of all extant and extinct
species of mammalian carnivores and related extinct groups
(Carnivoramorpha, Hyaenodonta and Oxyaenidae)retrieved from Paleobiology
Database and New and Old Worlds database.

## Loading data

``` r
library(dplyr)
df.faurby <- read.csv(here::here("data", "data_faurby.csv"), sep = ";") # dataset from https://royalsocietypublishing.org/doi/10.1098/rspb.2024.0473
df.body.mass <- read.csv(here::here("data", "Inferred_Size_Carnivoramorpha.csv"), sep = ",") # dataset from https://royalsocietypublishing.org/doi/10.1098/rspb.2024.0473
taxonomic_placement <- read.csv(here::here("data", "taxonomic_placement.csv"), sep = ";")
```

## General cleaning and processing

Here I will filter the original data keeping only Carnivora species from
suborders Caniformia and Feliformia and only records from Eurasia

``` r
taxonomic_placement2 <- 
  taxonomic_placement |> 
  filter(Suborder == "Caniformia" | Suborder == "Feliformia")


# adding a new column with genus and species to taxonomic placement and occurrence records
taxonomic_placement2$genus.species <- paste(taxonomic_placement2$Genus, taxonomic_placement2$Species, sep = "_")
df.faurby$genus.species <- paste(df.faurby$Genus, df.faurby$Species, sep = "_")


df.faurby2 <- 
  taxonomic_placement2 |> 
  left_join(df.faurby, by = c(genus.species = "genus.species")) |> 
  select(-Family.y, -Genus.y, -Species.y) |> 
  rename(Family = Family.x, Genus = Genus.x, Species = Species.x)

# keeping only Europe and NA
df.faurby3 <- 
  df.faurby2 |> 
  filter(Continent == "Europe" | Continent == "Asia")


df.faurby3 <- 
  df.faurby2 |> 
  filter(Continent == "Europe")
```

Keeping the site ID from PBDB. I will use NOW ID only when PBDB is
absent

``` r
df.faurby4 <- 
  df.faurby3 |> 
  mutate(site = ifelse(is.na(PBDB_ID), NOW_ID, PBDB_ID))
```

Processing data on body mass by calculating mean value for each species

``` r
df_body_mass <- 
  df.body.mass |> 
  group_by(Genus) |> 
  summarise(mean.size = mean(Size_Est))
```

## Processing and flagging data with `Buccaneer`

I will use the function `clean_occ_fossil` to make simple data cleaning,
processing and flagging the occurrence records. In this example I will
remove any record with NA in any of the columns and flag records with
the maximum and minimum age range greater than 5 myr, calculate midpoint
age, time of origination (TS) and time of extinction (TE)

``` r
df_occ_faurby <- 
  clean_occ_fossil(df.occ.fossil = df.faurby4, 
                 method.ages = "midpoint",
                 thresh.age.range = 5, 
                 species = "Genus",
                 Max.age = "Maximum_Age",
                 Min.age = "Minimum_Age", 
                 site = "site",
                 lat = "Latitude",
                 lng = "Longitude",
                 remove.sub.species = TRUE, 
                 comp.TS.TE = TRUE, 
                 group = "Suborder")
```

``` r
df_occ_faurby <- 
  df_occ_faurby |> 
  filter(species != "Delete")

# filtering TS and TE per species

df_TS_TE_faurby <- 
  df_occ_faurby |> 
  dplyr::distinct(species, .keep_all = T)

# joining TS and TE with body mass

df_TS_TE_mass <- 
  df_body_mass |> 
  left_join(df_TS_TE_faurby, by = c(Genus = "species")) |> 
  rename(species = "Genus") |> 
  filter(!is.na(TS) & !is.na(TE)) |> 
  filter(!is.na(site))

# creating time intervals (bins)

bins <- seq(from = max(df_occ_faurby$Max.age), to = min(df_occ_faurby$Min.age), by = -5)
```
