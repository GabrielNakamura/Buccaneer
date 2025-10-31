

# Loading Rodolfo's data

load(here::here("inst", "extdata", "script", "rodolfo", "longevities.RData"))
# load("./PBDB/df_can.Rdata")
load(here::here("inst", "extdata", "script", "rodolfo", "df_can.Rdata"))
load(here::here("inst", "extdata", "script", "rodolfo", "div_curves.Rdata"))
NALMA <- c("Duchesnean", "Chadronian","Orellan","Whitneyan","Arikareean",
           "Hemingfordian","Barstovian","Clarendonian","Hemphillian","Blancan",
           "Irvingtonian","Rancholabrean_Present") # Barnsosky 2014
NALMA_age <- c(39.7, 37, 33.9, 31.8, 29.5, 18.5, 16.3, 12.5, 9.4, 4.7, 1.4, 0.21,0) # Barnosky 2014



#####
# Using only replicate number 10 - actually are all the same
library(dplyr)

longs2 <- data.frame(species = rownames(longs[[20]]), longs[[20]])
res_rodolfo <-
  read.table(here::here("inst",
                        "extdata",
                        "script",
                        "rodolfo",
                        "site_res",
                        "site_diversity20.txt"),
             header = T)

  load(here::here("inst",
                        "extdata",
                        "script",
                        "rodolfo",
                        "div_curves.RData")
       )


df_occ_good <- readr::read_csv(here::here("inst", "extdata", "script", "rodolfo", "PBDB_NEW_Canidae_NA_simpler_PyRate_high_resol.csv"))

df_occ_good_low <-
  df_occ_good |>
  mutate(max_low_res = case_when(
    MaxT < 29.50 & MaxT > 18.50 ~ 29.5,
    MaxT < 16.30 & MaxT > 12.50 ~ 16.30,
    MaxT < 12.50 & MaxT > 9.40 ~ 12.50,
    MaxT < 9.40 & MaxT > 4.70 ~ 9.40,
    MaxT < 4.70 & MaxT > 1.40 ~ 4.70,
    MaxT < 1.40 & MaxT > 0.21 ~ 1.40,
    MaxT < 0.21 ~ 0.21,
    TRUE ~ MaxT
  )
  )


# this will be used to run the new analysis
df_occ_good_low2 <-
  df_occ_good_low |>
  mutate(min_low_res = case_when(
    MinT > 18.50 & MinT < 29.50 ~ 18.50,
    MinT > 12.50 & MinT < 16.30 ~ 12.50,
    MinT > 9.40 & MinT < 12.50 ~ 9.40,
    MinT > 4.70 & MinT < 9.40 ~ 4.70,
    MinT > 1.40 & MinT < 4.70 ~ 1.40,
    MinT > 0.21 & MinT < 1.40 ~ 0.21,
    MinT > 0.00 & MinT < 0.21 ~ 0,
    TRUE ~ MinT
  )) |>
  rename(species = Species) |>
  mutate(site.char = paste("site", Site, sep = "_")
         )

# checking the data that crosses the boundaries

crossing.nalma.occ <-
  lapply(NALMA_age,
         function(x) which(df_occ_good_low2$MaxT > x &
                             df_occ_good_low2$MinT < x)
  )

unlist(crossing.nalma.occ)

# data frame with occurrence records that cross nalma
df_crossing_nalma <- df_occ_good_low2[unlist(crossing.nalma.occ), ]

# removing the occurrences that cross nalma
df_occ_good_low_nocrossing <- df_occ_good_low2[-unlist(crossing.nalma.occ), ]

# calculating site coexistence withouth occ records crossing nalmas
res_clade_site_low_nocrossing <-
  clade_site_richness(df.TS.TE = longs2,
                      df.occ = df_occ_good_low_nocrossing,
                      time.slice = 0.1,
                      round.digits = 10,
                      species = "species",
                      TS = "TS",
                      TE = "TE",
                      Max.age = "max_low_res",
                      Min.age = "min_low_res",
                      site = "site.char",
                      group = NULL,
                      group.focal.compare = NULL,
                      type.comparison = NULL, remove.singletons = TRUE)


# calculating regional metrics with package function
res_regional_function <-
  clade_regional_richness(df.TS.TE = longs2,
                          time.slice = 0.1,
                          round.digits = 10,
                          species = "species",
                          TS = "TS",
                          TE = "TE")


# calculating site metrics using low resolution - same as used in Rodolfo
res_clade_site_low <-
  clade_site_richness(df.TS.TE = longs2,
                      df.occ = df_occ_good_low2,
                      time.slice = 0.1,
                      round.digits = 10,
                      species = "species",
                      TS = "TS",
                      TE = "TE",
                      Max.age = "max_low_res",
                      Min.age = "min_low_res",
                      site = "site.char",
                      group = NULL,
                      group.focal.compare = NULL,
                      type.comparison = NULL, remove.singletons = TRUE)


# using high resolution - different than used by rodolfo

res_clade_site_good_high <-
  clade_site_richness(df.TS.TE = longs2,
                      df.occ = df_occ_good_low2,
                      time.slice = 0.1,
                      round.digits = 1,
                      species = "species",
                      TS = "TS",
                      TE = "TE",
                      Max.age = "MaxT",
                      Min.age = "MinT",
                      site = "site.char",
                      group = NULL,
                      group.focal.compare = NULL,
                      type.comparison = NULL)


# plotting results from regional coexistence

plot(-as.numeric(names(div.curves[[50]])), div.curves[[50]], type = "l") # Rodolfo
lines(-res_regional_function$time.slice, res_regional_function$richness, type = "l", col = "red") # low resolution

# plotting results from site coexistence removing records crossing nalma
plot(-res_rodolfo$time, res_rodolfo$site_diversity, type = "l") # Rodolfo
lines(-res_clade_site_low_nocrossing$time.slice, res_clade_site_low_nocrossing$mean.coexistence, type = "l", col = "red") # low resolution


# plotting results from site coexistence
plot(-res_rodolfo$time, res_rodolfo$site_diversity, type = "l") # Rodolfo
lines(-res_clade_site_low$time.slice, res_clade_site_low$mean.coexistence, type = "l", col = "red") # low resolution
lines(-res_clade_site_low_nocrossing$time.slice, res_clade_site_low_nocrossing$mean.coexistence, type = "l", col = "green") # low resolution
abline(v = -33.9)
abline(v = -31.8)
abline(v = -20.44, col = "orange")
abline(v = -15.98, col = "orange")
abline(v = -16.3, col = "blue")
abline(v = -18.5, col = "blue")
abline(v = -1.4, col = "gray")


# plotting results from site coexistence but removing zeroes and self coex
plot(-res_rodolfo$time, res_rodolfo$site_diversity, type = "l") # Rodolfo
lines(-df_res_nozero$time.slice, df_res_nozero$mean.coexistence, type = "l", col = "red") # low resolution

# ploting with singletons
plot(-res_rodolfo$time, res_rodolfo$site_diversity, type = "l") # Rodolfo
lines(-df_res_singleton$time.slice, df_res_singleton$mean.coexistence, type = "l", col = "red") # low resolution

# plotting removing completely singleton
plot(-res_rodolfo$time, res_rodolfo$site_diversity, type = "l") # Rodolfo
lines(-df_res_removingsingleton$time.slice, df_res_removingsingleton$mean.coexistence, type = "l", col = "red") # low resolution



plot(-res_rodolfo$time, res_rodolfo$site_diversity, type = "l") # Rodolfo
lines(-df_res$time.slice, df_res$mean.coexistence, type = "l", col = "red") # low resolution

# lines(-smry.site.nalma.t[[10]]$id, smry.site.nalma.t[[10]]$mean, col = "green")
lines(-res_clade_site_low_digit10$time.slice, res_clade_site_low_digit10$mean.coexistence, col = "orange")
abline(v = -res_clade_site_low_digit10$time.slice, lwd = 0.05)
lines(-res_clade_site_good_high$time.slice, res_clade_site_good_high$mean.coexistence, type = "l", col = "blue") # high resolution
#lines(-t3_clade_site$time.slice, t3_clade_site$mean.coexistence, type = "l", col = "blue")

# only checking
unique(df_can_tot2$max)
unique(df_can_tot2$min)

sort(unique(df_can_tot2$max), decreasing = T)
sort(unique(df_can_tot2$min), decreasing = T)
NALMA_age
