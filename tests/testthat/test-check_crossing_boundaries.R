test_that("check_crossing_boundaries analyzes boundary crossing", {
  # Create test data
  df_test <- data.frame(
    species = c("sp1", "sp2", "sp1", "sp2"),
    Maximum_Age = c(100, 95, 90, 85),
    Minimum_Age = c(80, 75, 70, 65)
  )
  
  # Define time intervals
  interval <- c(100, 80, 60)
  
  # Run the function
  result <- check_crossing_boundaries(
    df.occ.fossil = df_test,
    interval = interval,
    species = "species",
    Max.age = "Maximum_Age",
    Min.age = "Minimum_Age"
  )
  
  # Test that result is a list
  expect_type(result, "list")
  
  # Test that result has expected elements
  expect_true("summary_table" %in% names(result) || "df_records_interval" %in% names(result))
})
