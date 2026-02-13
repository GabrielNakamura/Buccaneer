test_that("IndivSpecies_reach_coexistence calculates individual species reach coexistence", {
  # Load test data
  data("df_longevities_canidae")
  data("df_occurrence_canidae")

  # Prepare data
  df_ts_te <- df_longevities_canidae
  df_occ <- df_occurrence_canidae

  # Rename columns
  colnames(df_occ)[colnames(df_occ) == "MaxT"] <- "Max.age"
  colnames(df_occ)[colnames(df_occ) == "MinT"] <- "Min.age"
  colnames(df_occ)[colnames(df_occ) == "lng"] <- "lng"
  colnames(df_occ)[colnames(df_occ) == "lat"] <- "lat"

  # Run the function
  result <- IndivSpecies_reach_coexistence( species = "species",
    df.TS.TE = df_ts_te,
    TS = "TS",
    TE = "TE",
    df.occ = df_occ,
    time.slice = 5,
    crs = 4326
  )

  # Test that result is a data frame
  expect_s3_class(result, "data.frame")

  # Test that result has expected columns
  expect_true("species" %in% colnames(result))
  expect_true("n.coexistence" %in% colnames(result) || "time.slice" %in% colnames(result))
})
