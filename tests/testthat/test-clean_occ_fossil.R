test_that("clean_occ_fossil cleans and flags occurrence data", {
  # Create test occurrence data
  df_test <- data.frame(
    species = c("sp1", "sp2", "sp1", "sp2"),
    Maximum_Age = c(100, 95, 90, 85),
    Minimum_Age = c(80, 75, 70, 65)
  )
  
  # Run the function
  result <- clean_occ_fossil(
    df.occ.fossil = df_test,
    thresh.age.range = 20,
    species = "species",
    Max.age = "Maximum_Age",
    Min.age = "Minimum_Age",
    comp.TS.TE = TRUE
  )
  
  # Test that result is a data frame
  expect_s3_class(result, "data.frame")
  
  # Test that result has added columns
  expect_true("midpoint" %in% colnames(result))
  expect_true("age.range" %in% colnames(result))
  expect_true("flag.age.range" %in% colnames(result))
})
