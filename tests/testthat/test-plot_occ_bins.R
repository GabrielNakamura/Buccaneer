test_that("plot_occ_bins creates a ggplot object", {
  # Create test data
  df_test <- data.frame(
    accepted_name = rep(c("sp1", "sp2", "sp3"), each = 3),
    max_ma = c(100, 95, 90, 85, 80, 75, 70, 65, 60),
    min_ma = c(95, 90, 85, 80, 75, 70, 65, 60, 55)
  )
  
  # Run the function
  result <- plot_occ_bins(
    df.occ.fossil = df_test,
    species = "accepted_name",
    max.age = "max_ma",
    min.age = "min_ma",
    age.scheme = FALSE
  )
  
  # Test that result is a ggplot object
  expect_s3_class(result, "ggplot")
})
