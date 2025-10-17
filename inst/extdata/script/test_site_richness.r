#' Compute time series based on mean site richness
#'
#' @param df.TS.TE Data frame object containing at least four columns. Species names,
#'     origination time, extinction time and a trait value for each species.
#' @param df.occ a data frame object containing the occurrence records for each species.
#'     This must have at least a column indicating the name of species, its minimum and maximum age estimate,
#'     and its site location ID.
#' @param time.slice Scalar indicating the time interval between consecutive time slices.
#' @param round.digits Scalar indicating the precision of time slices.
#' @param species Character indicating the name of the column of the data frame
#'     containing the species name information.
#' @param TS Character indicating the name of the columns of the data frame
#'     containing the information on origination time.
#' @param TE Character indicating the name of the column of the data frame
#'     containing the information on extinction time.
#' @param group Character indicating the name of the column that contain the groups that will be used
#'     in comparison.
#' @param group.focal.compare Character vector indicating the focal (first element) and comparison (second element)
#'     groups used in the calculation. If NULL, the default, the metrics  will be calculated
#'     using all  clades.
#' @param type.comparison Character. It can be "between" to compute distances only between species/genus of two groups
#'     or "within" to calculate distance only inside the focal group. If null the distance is computed
#'     considering all clades together
#' @param Max.age Character indicating the name of the column containing the upper age limit for occurrence record.
#' @param Min.age Character indicating the name of the column containing the lower age limit for occurrence record.
#' @param site Character indicating the name of the column containing the information on site location.
#'
#' @return
#' @export
#'
#' @examples
clade_site_richness <-
  function(df.TS.TE,
           df.occ,
           time.slice,
           round.digits = 1,
           species = "species",
           TS = "TS",
           TE = "TE",
           Max.age = "Max.age",
           Min.age = "Min.age",
           site = "site",
           group = NULL,
           group.focal.compare = NULL,
           type.comparison = NULL){

    # sub-setting the dataframe to keep or remove the group information

    if(!is.null(group) == TRUE){
      df.TS.TE <- df.TS.TE[, c(species, TS, TE, group)]
      colnames(df.TS.TE) <- c("species", "TS", "TE", "group")
    } else{
      df.TS.TE <- df.TS.TE[, c(species, TS, TE)]
      colnames(df.TS.TE) <- c("species", "TS", "TE")
    }

    # sub-setting and correcting the names of occurrence dataframe
    df_occ <-
      df.occ[, c(species, Max.age, Min.age, site)]
    vars <- list(species, Max.age, Min.age, site)
    name_vars <- c("species", "Max.age", "Min.age", "site")
    names(vars) <- name_vars
    column.names <- names(unlist(vars))
    colnames(df_occ) <- column.names
    df_occ$site <- as.factor(df_occ$site)

    # Generating time intervals used to compute temporal coexistence
    seq_interval <- seq(from = max(df.TS.TE[, "TS"]), to = min(df.TS.TE[, "TE"]), by = -time.slice)
    seq_interval <- c(round(seq_interval, digits = round.digits))

    # Time coexistence matrix for all species
    matrix_coex <-
      aux_matrix_regional_coex(df.TS.TE, time.slice, round.digits = 1,
                               species = "species",
                               TS = "TS",
                               TE = "TE")

    # modified coexistence matrix containing group comparison
    if(!is.null(group.focal.compare) == TRUE){
      focal <- group.focal.compare[1]
      compare <- group.focal.compare[2]
      spp_focal <- df.TS.TE[which(df.TS.TE$group == focal), "species"]$species
      spp_compare <- df.TS.TE[which(df.TS.TE$group == compare), "species"]$species

      if(type.comparison == "between"){# comparison between species of two groups
        matrix_coex <- lapply(matrix_coex, function(x) x[spp_focal, spp_compare]) # focal species in lines and comparison in columns
      }
      if(type.comparison == "within"){ # comparison only within the focal group
        matrix_coex <- lapply(matrix_coex, function(x) x[spp_focal, spp_focal])
      }
    } else{
      matrix_coex <- matrix_coex
    }

    # species composition at each timeslice
    spp_slice <-
      lapply(matrix_coex, function(x){
        names(which(rowSums(x) >= 1))
      })

    names(spp_slice) <- seq_interval

    # calculating matrix of species cooccurrence for site

    list_matrix_cooccur_site <-
      comp_site_cooccurr(spp_slice = spp_slice, df.occ = df_occ)

    # only presence-absence
    list_matrix_cooccur_site2 <- lapply(list_matrix_cooccur_site, function(x) ifelse(x >= 1, 1, 0))

    # calculating mean coexistence and variance in coexistence for each species in each timeslice

    mean_species_site_coex <-
      lapply(list_matrix_cooccur_site2,
             function(x){
               mean((rowSums(x) - 1))
             })

    var_species_site_coex <-
      lapply(list_matrix_cooccur_site2,
             function(x){
               var((rowSums(x) - 1))
             })

    # function output

    df_res <-
      data.frame(mean.coexistence = unlist(mean_species_site_coex),
                 var.distance = unlist(var_species_site_coex),
                 time.slice = seq_interval)

    return(df_res)


  }




#####
# Using only replicate number 10 - actually are all the same

longs2 <- data.frame(species = rownames(longs[[50]]), longs[[50]])
res_rodolfo <-
  read.table(here::here("inst",
                        "extdata",
                        "script",
                        "rodolfo",
                        "site_res",
                        "site_diversity50.txt"),
             header = T)

  load(here::here("inst",
                        "extdata",
                        "script",
                        "rodolfo",
                        "div_curves.RData")
       )


df_occ_good <- read_csv(here::here("inst", "extdata", "script", "rodolfo", "PBDB_NEW_Canidae_NA_simpler_PyRate_high_resol.csv"))

df_occ_good_low <-
  df_occ_good %>%
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

# using low resolution

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
                      type.comparison = NULL)

res_clade_site_low_digit10 <-
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
                      type.comparison = NULL)

# using high resolution

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

# without zeros

res_clade_site_good_low_nozero <-
  clade_site_richness2(df.TS.TE = longs2,
                      df.occ = df_occ_good_low2,
                      time.slice = 0.1,
                      round.digits = 1,
                      species = "species",
                      TS = "TS",
                      TE = "TE",
                      Max.age = "max_low_res",
                      Min.age = "min_low_res",
                      site = "site.char",
                      group = NULL,
                      group.focal.compare = NULL,
                      type.comparison = NULL)

res_regional_function <-
  clade_regional_richness(df.TS.TE = longs2,
                          time.slice = 0.1,
                          round.digits = 10,
                          species = "species",
                          TS = "TS",
                          TE = "TE")


quartz()
plot(-div.curves.frame[[50]]$time, div.curves.frame[[50]]$div, type = "l") # Rodolfo
lines(-res_regional_function$time.slice, res_regional_function$richness, type = "l", col = "red") # low resolution


quartz()
plot(-res_rodolfo$time, res_rodolfo$site_diversity, type = "l") # Rodolfo
lines(-res_clade_site_low$time.slice, res_clade_site_low$mean.coexistence, type = "l", col = "red") # low resolution

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
