test_that("get_summary_interval summarizes occurrences by interval", {
  # Create test data
  df_test <- data.frame(
    species = c("sp1", "sp2", "sp1", "sp2"),
    Maximum_Age = c(100, 95, 90, 85),
    Minimum_Age = c(80, 75, 70, 65),
    midpoint = c(90, 85, 80, 75)
  )
  
  # Define intervals
  interval <- c(100, 80, 60)
  
  # Run the function
  result <- get_summary_interval(
    df.occ.fossil = df_test,
    interval = interval,
    species = "species",
    Max.age = "Maximum_Age",
    Min.age = "Minimum_Age"
  )
  
  # Test that result is a list
  expect_type(result, "list")
})
