test_that("Assemblage_site_trait_distance calculates trait distance", {
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
  colnames(df_occ)[colnames(df_occ) == "Site"] <- "site"

  # Create a simple distance matrix for trait
  n_spp <- nrow(df_ts_te)
  dist_mat <- dist(df_ts_te$log_mass)

  # Run the function
  result <- assemblage_site_trait_distance(
    df.TS.TE = df_ts_te,
    df.occ = df_occ,
    time.slice = 5,
    dist.trait = dist_mat,
    trait = "log_mass"
  )

  # Test that result is a data frame
  expect_s3_class(result, "data.frame")

  # Test that result has expected columns
  expect_true("mean_dist_to_cooccur" %in% colnames(result) || "sites" %in% colnames(result) || "time.slice" %in% colnames(result))

  # Test that result has any numeric value
  expect_false(all(is.na(result$mean_dist_to_cooccur)))
})
