test_that("IndivSpecies_reach_distance calculates individual species reach distance", {
  # Load test data
  data("df_longevities_canidae")
  data("traits_canidae")
  data("df_occurrence_canidae")

  # Prepare data
  df_ts_te <- merge(df_longevities_canidae, traits_canidae, by = "species")
  df_occ <- df_occurrence_canidae

  # Rename columns
  colnames(df_occ)[colnames(df_occ) == "MaxT"] <- "Max.age"
  colnames(df_occ)[colnames(df_occ) == "MinT"] <- "Min.age"
  colnames(df_occ)[colnames(df_occ) == "lng"] <- "lng"
  colnames(df_occ)[colnames(df_occ) == "lat"] <- "lat"

  # Create a distance matrix
  dist_mat <- dist(df_ts_te$log_mass)

  # Run the function
  result <- IndivSpec_reach_distance(
    df.TS.TE = df_ts_te,
    df.occ = df_occ,
    time.slice = 5,
    dist.trait = dist_mat,
    crs = 4326, lon = "lng",
    lat = "lat",
    nearest.taxon = 1
  )

  # Test that result is a data frame
  expect_s3_class(result, "data.frame")

  # Test that result has expected columns
  expect_true("species" %in% colnames(result) || "time.slice" %in% colnames(result))
})
