test_that("IndivSpecies_regional_coexistence calculates individual species regional coexistence", {
  # Load test data
  data("df_longevities_canidae")
  data("df_occurrence_canidae")

  # Prepare data
  df_ts_te <- df_longevities_canidae
  df_occ <- df_occurrence_canidae

  # Rename columns
  colnames(df_occ)[colnames(df_occ) == "MaxT"] <- "Max.age"
  colnames(df_occ)[colnames(df_occ) == "MinT"] <- "Min.age"
  colnames(df_occ)[colnames(df_occ) == "Site"] <- "site"

  # Run the function
  result <- IndivSpec_regional_coex(
    df.TS.TE = df_ts_te,
    time.slice = 5
  )

  # Test that result is a data frame
  expect_s3_class(result, "data.frame")

  # Test that result has expected columns
  expect_true("species" %in% colnames(result) || "n.coexistence" %in% colnames(result))
})
