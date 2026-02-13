test_that("IndivSpecies_regional_distance calculates individual species regional distance", {
  # Load test data
  data("df_longevities_canidae")
  data("traits_canidae")

  # Prepare data
  df_ts_te <- merge(df_longevities_canidae, traits_canidae, by = "species")

  # Create a distance matrix
  dist_mat <- dist(df_ts_te$log_mass)

  # Run the function
  result <- IndivSpec_regional_distance(
    df.TS.TE = df_ts_te,
    time.slice = 5,
    dist.trait = dist_mat,
    nearest.taxon = 1
  )

  # Test that result is a data frame
  expect_s3_class(result, "data.frame")

  # Test that result has expected columns
  expect_true("species" %in% colnames(result) || "time.slice" %in% colnames(result))
})
