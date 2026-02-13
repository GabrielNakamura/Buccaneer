test_that("Assemblage_regional_richness calculates regional richness", {
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
  result <- Assemblage_regional_richness(
    df.TS.TE = df_ts_te,
    df.occ = df_occ,
    time.slice = 5,
    grid.size = 10
  )
  
  # Test that result is a list
  expect_type(result, "list")
  
  # Test that result has expected elements
  expect_true("grid_mean_age" %in% names(result) || "time_series_rich" %in% names(result))
})
