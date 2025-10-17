df_toy <- data.frame(
  species = paste0("Species", 1:6),
  trait1 = c(1, 1, 1, 5, 5, 4),
  trait2 = c(2, 1, 3, 2, 1, 3),
  TS = c(10, 8, 12, 5, 9, 7),
  TE = c(2,  1,  6, 0, 4, 3),
  site = c(1, 2, 1, 2, 3, 3))
df_ts_te <- df_toy[, c("species", "TS", "TE")]
df_occ <- df_toy[, c("species", "TS", "TE", "site")]
names(df_occ) <- c("species", "Max.age", "Min.age", "site")

devtools::load_all()

res_clade_toy <-
  clade_site_richness(df.TS.TE = df_ts_te,
                      df.occ = df_occ,
                      time.slice = 1,
                      round.digits = 1,
                      species = "species",
                      TS = "TS",
                      TE = "TE",
                      Max.age = "Max.age",
                      Min.age = "Min.age",
                      site = "site")


res_clade_faurby <-
  clade_site_richness(df.TS.TE = df_TS_TE_faurby,
                      df.occ = df_occ_faurby,
                      time.slice = 0.1,
                      round.digits = 1,
                      species = "species",
                      TS = "TS",
                      TE = "TE",
                      Max.age = "Max.age",
                      Min.age = "Min.age",
                      site = "site")

