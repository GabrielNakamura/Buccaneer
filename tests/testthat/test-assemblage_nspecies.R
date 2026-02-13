test_that("assemblage_nspecies calculates species count", {
  # Load test data
  data("df_longevities_canidae")
  data("df_occurrence_canidae")
  
  # Prepare data with required column names
  df_ts_te <- df_longevities_canidae
  df_occ <- df_occurrence_canidae
  
  # Rename columns to match function expectations
  colnames(df_occ)[colnames(df_occ) == "MaxT"] <- "Max.age"
  colnames(df_occ)[colnames(df_occ) == "MinT"] <- "Min.age"
  colnames(df_occ)[colnames(df_occ) == "Site"] <- "site"
  
  # Run the function
  result <- assemblage_nspecies(
    df.TS.TE = df_ts_te,
    df.occ = df_occ,
    time.slice = 5
  )
  
  # Test that result is a data frame
  expect_s3_class(result, "data.frame")
  
  # Test that result has expected columns
  expect_true("time.slice" %in% colnames(result))
  expect_true("sites" %in% colnames(result))
  expect_true("n.species" %in% colnames(result))
  
  # Test that n.species contains numeric values
  expect_true(is.numeric(result$n.species) | all(is.na(result$n.species)))
})
